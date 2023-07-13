import csv
import itertools
import logging
import multiprocessing
import os
import random
import statistics
import sys

import matplotlib.pyplot as plt
import pyfastaq

from sprog import bifrost, fastani, mash, sample, tree, utils

csv.field_size_limit(sys.maxsize)


def _match_unitigs(species, ref_files, query_file, outprefix, k, kmer_ratio, debug):
    array_length, got_tsv = bifrost.build_and_query(
        ref_files,
        query_file,
        outprefix,
        k=k,
        kmer_ratio=kmer_ratio,
        debug=debug,
    )
    return (species, got_tsv, array_length)


class SampleSet:
    def __init__(self, root_dir, debug=False):
        self.root_dir = os.path.abspath(root_dir)
        self.samples_root = os.path.join(self.root_dir, "Samples")
        self.all_unitigs_dir = os.path.join(self.root_dir, "All_unitigs")
        self.all_unitigs_fa_prefix = os.path.join(self.all_unitigs_dir, "all_unitigs")
        self.all_unitigs_fa = f"{self.all_unitigs_fa_prefix}.fasta.gz"
        self.all_unitigs_fofn = os.path.join(self.all_unitigs_dir, "genome_list.txt")
        self.all_query_dir = os.path.join(self.root_dir, "All_queries")
        self.all_query_json = os.path.join(self.all_query_dir, "files.json")
        self.mash_dir = os.path.join(self.root_dir, "Mash")
        self.mash_distances = None
        self.fastani_dir = os.path.join(self.root_dir, "FastANI")
        self.fastani_json = os.path.join(self.root_dir, "FastANI", "fastANI.json")
        self.mash_fastani_dir = os.path.join(self.root_dir, "Mash_v_fastANI")
        self.probes_dir = os.path.join(self.root_dir, "Probes")
        self.node2unitig_tsv = os.path.join(self.probes_dir, "node_to_unitig.tsv")
        self.unitig_counts_tree = os.path.join(
            self.probes_dir, "tree_with_unitig_counts.txt"
        )
        self.node2unitigs = None
        self.probes_fa = os.path.join(self.probes_dir, "probes.fa.gz")
        self.myk_hierarchy_json = os.path.join(
            self.probes_dir, "mykrobe_hierarchy.json"
        )

        if not os.path.exists(self.root_dir):
            os.mkdir(self.root_dir)
            os.mkdir(self.samples_root)
        self.metadata_json = os.path.join(root_dir, "metadata.json")
        self.debug = debug
        if os.path.exists(self.metadata_json):
            self.metadata = utils.load_json(self.metadata_json)
        else:
            self.metadata = {}

        self.samples = {}
        for sample_name in self.metadata:
            sample_dir = os.path.join(self.samples_root, sample_name)
            self.samples[sample_name] = sample.Sample(sample_dir, debug=self.debug)

    def add_samples_from_tsv(self, tsv_file):
        with open(tsv_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for d in reader:
                if d["sample"] in self.metadata:
                    logging.info(f"Add samples. Already got {d['sample']}. Skipping")
                    continue

                logging.info(f"Add sample {d['sample']}")
                new_sample = sample.Sample(
                    os.path.join(self.samples_root, d["sample"]),
                    name=d["sample"],
                    data_source=d["data_source"],
                    accession=None if d["accession"] == "." else d["accession"],
                    species=d["species"],
                    original_file=None if d["file"] in [None, "."] else d["file"],
                    debug=self.debug,
                )

                self.metadata[d["sample"]] = {
                    "species": d["species"],
                    "dir": os.path.relpath(new_sample.root_dir, self.root_dir),
                }

                self.samples[d["sample"]] = new_sample

        utils.write_json(self.metadata, self.metadata_json)
        logging.info(f"Finished loading samples. Total samples: {len(self)}")

    def __len__(self):
        return len(self.metadata)

    def download_missing_genomes(self):
        for s in self.samples.values():
            s.obtain_genome()

    def get_species2sample(self):
        species2sample = {}
        for name, d in self.metadata.items():
            if d["species"] not in species2sample:
                species2sample[d["species"]] = []
            species2sample[d["species"]].append(self.samples[name].genome_fasta)
        return species2sample

    def fastani_all_v_all(self, force=False, threads=1, fraglen=None, kmer=None):
        genomes = {s.genome_fasta: n for n, s in self.samples.items()}
        name2species = {n: s.metadata["species"] for n, s in self.samples.items()}
        distances_tsv = fastani.all_v_all_fastani(
            genomes,
            name2species,
            self.fastani_dir,
            force=force,
            threads=threads,
            fraglen=fraglen,
            kmer=kmer,
        )
        distances = fastani.load_distances_tsv(distances_tsv)
        utils.write_json(distances, self.fastani_json)
        fastani.analyse_ani_data(distances, os.path.join(self.fastani_dir, "analysis"))
        return distances

    def mash_all_v_all(self, force=False):
        genomes = {s.genome_fasta: n for n, s in self.samples.items()}
        name2species = {n: s.metadata["species"] for n, s in self.samples.items()}
        distances_tsv = mash.all_v_all_mash(
            genomes, name2species, self.mash_dir, force=force
        )
        self.mash_distances = mash.load_distances_tsv(distances_tsv)
        mash.analyse_distance_data(
            self.mash_distances, os.path.join(self.mash_dir, "analysis")
        )

    def fastani_mash_scatter(self):
        mash_tsv = os.path.join(self.mash_dir, "distances.tsv")
        if not (os.path.exists(mash_tsv) and os.path.exists(self.fastani_json)):
            return

        if os.path.exists(self.mash_fastani_dir):
            utils.syscall(f"rm -rf {self.mash_fastani_dir}")

        fastani_data = utils.load_json(self.fastani_json)
        os.mkdir(self.mash_fastani_dir)
        x_points = []
        y_points = []

        with open(mash_tsv) as f:
            for line in f:
                genome1, genome2, distance, *_ = line.split("\t")
                key = ";".join(sorted([genome1, genome2]))
                if key in fastani_data:
                    identities = [
                        v
                        for k, v in fastani_data[key].items()
                        if k.startswith("identity_") and v is not None
                    ]
                    if len(identities) == 0:
                        identity = 0
                    else:
                        identity = statistics.mean(identities)
                else:
                    identity = 0
                x_points.append(identity)
                y_points.append(float(distance))

        fig, ax = plt.subplots()
        hb = ax.hexbin(
            x_points, y_points, bins="log", mincnt=1, cmap="rainbow", gridsize=50
        )
        plt.colorbar(hb)
        ax.set_xlabel("ANI")
        ax.set_ylabel("Mash distance")
        plt.savefig(os.path.join(self.mash_fastani_dir, "scatter.pdf"))
        plt.close()

        fig, ax = plt.subplots()
        hb = ax.hexbin(
            [100 - x for x in x_points],
            y_points,
            bins="log",
            mincnt=1,
            cmap="rainbow",
            gridsize=50,
        )
        plt.colorbar(hb)
        ax.set_xlabel("100 - ANI")
        ax.set_ylabel("Mash distance")
        plt.savefig(os.path.join(self.mash_fastani_dir, "scatter.2.pdf"))
        plt.close()

    def make_unitigs_one_graph_all_samples(self, force=False, threads=1, k=None):
        if os.path.exists(self.all_unitigs_dir):
            if force:
                logging.info(
                    "--force used: deleting existing unitigs directory if it exists"
                )
                utils.syscall(f"rm -r {self.all_unitigs_dir}")
            else:
                raise Exception(
                    f"Cannot build graph from all samples. Unitigs directory already exists ({self.all_unitigs_dir})\nYou could use the --force option to replace existing unitigs directory"
                )

        logging.info("Gathering sample info to build unitigs")
        os.mkdir(self.all_unitigs_dir)
        with open(self.all_unitigs_fofn, "w") as f:
            for sample_name, this_sample in self.samples.items():
                if not this_sample.has_genome():
                    raise Exception(
                        f"Sample {sample_name} has not had genome downloaded. Cannot continue"
                    )
                print(this_sample.genome_fasta, file=f)

        logging.info(
            f"Start building graph of unitigs from {len(self)} samples using {threads} threads"
        )
        got_fasta = bifrost.build(
            self.all_unitigs_fofn,
            self.all_unitigs_fa_prefix,
            "fasta",
            threads=threads,
            k=k,
        )
        assert got_fasta == self.all_unitigs_fa
        logging.info("Finished building graph of unitigs")

    def query_all_unitigs_to_each_species(
        self, threads=1, force=False, ratio_kmers=None, k=None
    ):
        if os.path.exists(self.all_query_dir):
            if force:
                logging.info(
                    "--force used: deleting existing queries directory if it exists"
                )
                utils.syscall(f"rm -r {self.all_query_dir}")
            else:
                raise Exception(
                    f"Cannot query unitigs. Directory already exists ({self.all_query_dir})\nYou could use the --force option to replace existing directory"
                )

        if not os.path.exists(self.all_unitigs_fa):
            raise Exception(f"Unitigs of all samples not found: {self.all_unitigs_fa}")
        if not os.path.exists(self.all_unitigs_fofn):
            raise Exception(
                f"File of genomes used to make all unitigs not found: {self.all_unitigs_fofn}"
            )

        expect_fastas = len([1 for k, v in self.samples.items() if v.has_genome()])

        with open(self.all_unitigs_fofn) as f:
            got_fastas = len([x for x in f])

        if expect_fastas != got_fastas:
            raise Exception(
                f"Different number of samples found compared to the number used to make all unitigs fasta {self.all_unitigs_fa}"
            )

        species2sample = self.get_species2sample()
        species, ref_files = zip(*species2sample.items())
        os.mkdir(self.all_query_dir)
        outprefixes = [
            os.path.join(self.all_query_dir, str(i)) for i in range(len(species))
        ]

        with multiprocessing.Pool(processes=threads) as pool:
            results = pool.starmap(
                _match_unitigs,
                zip(
                    species,
                    ref_files,
                    itertools.repeat(self.all_unitigs_fa),
                    outprefixes,
                    itertools.repeat(k),
                    itertools.repeat(ratio_kmers),
                    itertools.repeat(self.debug),
                ),
            )
        all_lengths = set(t[2] for t in results)
        assert len(all_lengths) == 1
        number_of_unitigs = all_lengths.pop()
        results = {t[0]: os.path.basename(t[1]) for t in results}
        to_write = {
            "number_of_unitigs": number_of_unitigs,
            "filenames": results,
        }
        utils.write_json(to_write, self.all_query_json)
        logging.info("Finished processing all samples")

    def load_species_query_files(self):
        json_data = utils.load_json(self.all_query_json)
        number_of_unitigs = json_data["number_of_unitigs"]
        unitig_presence = {}

        for species, presence_file in json_data["filenames"].items():
            filename = os.path.join(self.all_query_dir, presence_file)
            unitig_presence[species] = bifrost.load_bitarray_file(
                filename, number_of_unitigs
            )
            if len(unitig_presence) % 25 == 0:
                logging.info(
                    f"Loaded {len(unitig_presence)} of {len(json_data['filenames'])} species"
                )

        return unitig_presence

    def match_probes_to_tree_nodes(self, species_tree, unitig_presence=None):
        if unitig_presence is None:
            unitig_presence = self.load_species_query_files()
        species_names = list(unitig_presence.keys())
        leaf_combos = species_tree.leaf_combinations_to_nodes(species_names)

        total_unitigs = len(unitig_presence[species_names[0]])
        logging_number = int(total_unitigs / 10)
        node2unitigs = {}
        logging.info(f"Processing {total_unitigs} unitigs")
        for i in range(total_unitigs):
            if i % logging_number == 0 and i > 0:
                pc = round(100 * i / total_unitigs, 1)
                logging.info(f"Processed {pc}% ({i}/{total_unitigs}) unitigs")
            key = tuple(
                j for j, name in enumerate(species_names) if unitig_presence[name][i]
            )
            if key in leaf_combos:
                node = leaf_combos[key]
                if node not in node2unitigs:
                    node2unitigs[node] = set()
                node2unitigs[node].add(i)

        return node2unitigs

    def match_probes_to_tree_nodes_2(
        self,
        tree_file,
        min_prop_contain=1,
        max_prop_outside=0.0,
        force=False,
    ):
        if os.path.exists(self.probes_dir):
            if force:
                logging.info(
                    "--force used: deleting existing unitigs directory if it exists"
                )
                utils.syscall(f"rm -r {self.probes_dir}")
            else:
                raise Exception(
                    f"Cannot match unitigs to tree. Directory already exists ({self.probes_dir})\nYou could use the --force option to replace existing directory"
                )

        if not os.path.exists(self.all_query_json):
            raise Exception(
                "Did not find results of running `sprog per_sample_presence`"
            )

        os.mkdir(self.probes_dir)
        unitig_matrix = self.load_species_query_files()
        species_names = list(unitig_matrix.keys())
        species_tree = tree.Tree(tree_file)
        species_tree.init_leaf_combination_data(species_names)

        total_unitigs = len(unitig_matrix[species_names[0]])
        logging_number = max(1, int(total_unitigs / 20))
        self.node2unitigs = {}
        logging.info(f"Processing {total_unitigs} unitigs")
        for i in range(total_unitigs):
            if i % logging_number == 0 and i > 0:
                pc = round(100 * i / total_unitigs, 1)
                logging.info(f"Processed {pc}% ({i}/{total_unitigs}) unitigs")
            combo = {
                j for j, name in enumerate(species_names) if unitig_matrix[name][i]
            }
            node = species_tree.leaf_combo_parent_nodes(
                combo,
                min_prop_contain=min_prop_contain,
                max_prop_outside=max_prop_outside,
            )
            if node is None:
                continue
            if node not in self.node2unitigs:
                self.node2unitigs[node] = set()
            self.node2unitigs[node].add(i)

        species_tree.write_to_file(
            self.unitig_counts_tree, node_unitig_combos=self.node2unitigs
        )

        self.write_node2unitig_tsv()

    def write_node2unitig_tsv(self):
        with open(self.node2unitig_tsv, "w") as f:
            print("Node", "Unitigs", sep="\t", file=f)
            for node, unitigs in self.node2unitigs.items():
                print(
                    node.name,
                    ",".join(sorted(str(x) for x in unitigs)),
                    sep="\t",
                    file=f,
                )

    def load_node2untig_tsv(self, species_tree):
        self.node2unitigs = {}
        logging.info(f"Loading file {self.node2unitig_tsv}")
        with open(self.node2unitig_tsv) as f:
            for d in csv.DictReader(f, delimiter="\t"):
                node = species_tree.find_node_by_name(d["Node"])
                assert node is not None
                self.node2unitigs[node] = set(int(x) for x in d["Unitigs"].split(","))

    def write_probes(
        self,
        tree_file,
        max_probes_per_node=2000,
        mykrobe_tb=False,
    ):
        random.seed(42)

        species_tree = tree.Tree(tree_file)
        if mykrobe_tb:
            species_tree.write_mykrobe_hierarchy_json(self.myk_hierarchy_json)

        if self.node2unitigs is None:
            self.load_node2untig_tsv(species_tree)

        logging.info("Gathering unitigs/probes for each node in the tree")
        wanted_unitigs = set()
        wanted_node2unitig = {}
        for node, unitigs in self.node2unitigs.items():
            if len(unitigs) <= max_probes_per_node:
                wanted_node2unitig[node] = unitigs
            else:
                wanted_node2unitig[node] = random.sample(unitigs, max_probes_per_node)

            wanted_unitigs.update(wanted_node2unitig[node])

        logging.info(
            f"Loading {len(wanted_unitigs)} required unitigs from file {self.all_unitigs_fa}"
        )
        reader = pyfastaq.sequences.file_reader(self.all_unitigs_fa)
        unitig_seqs = {int(s.id): s.seq for s in reader if int(s.id) in wanted_unitigs}
        logging.info(f"Total unitigs loaded: {len(unitig_seqs)}")

        if mykrobe_tb:
            edges_to_panel_type = [
                None,  # root, don't expect any
                "sub-complex",
                "species",
                "sub-species",
            ]

        logging.info(f"Writing probes to file {self.probes_fa}")
        f = pyfastaq.utils.open_file_write(self.probes_fa)

        for node, unitigs in wanted_node2unitig.items():
            if mykrobe_tb:
                panel_type = edges_to_panel_type[species_tree.edges_to_root(node)]
                assert panel_type is not None
            for i in unitigs:
                if mykrobe_tb:
                    # example name from tb-species-170421.fasta.gz:
                    # Mycobacterium_abscessus?name=Mycobacterium_abscessus&panel_type=species&length=40
                    name = "_".join(node.name.split())
                    length = len(unitig_seqs[i])
                    print(
                        f">{name}?name={name}&panel_type={panel_type}&length={length}",
                        unitig_seqs[i],
                        sep="\n",
                        file=f,
                    )
                else:
                    print(
                        f">{'_'.join(node.name.split())}_{i}",
                        unitig_seqs[i],
                        sep="\n",
                        file=f,
                    )
        f.close()
