import copy
import json
import logging
import os

import matplotlib.pyplot as plt

from sprog import utils


def all_v_all_fastani(
    genome2name,
    name2species,
    outdir,
    fraglen=None,
    kmer=None,
    threads=1,
    force=False,
):
    if os.path.exists(outdir):
        if force:
            utils.syscall(f"rm -rf {outdir}")
        else:
            raise Exception(f"Output directory already exists {outdir}")

    os.mkdir(outdir)

    genomes_fofn = os.path.join(outdir, "genomes.fofn")
    fastani_out = os.path.join(outdir, "fastANI.out.tsv")

    with open(genomes_fofn, "w") as f:
        print(*genome2name.keys(), sep="\n", file=f)

    logging.info("Running all vs all fastANI")
    fraglen = "" if fraglen is None else f"--fragLen {fraglen}"
    kmer = "" if kmer is None else f"--kmer {kmer}"
    utils.syscall(
        f"fastANI {fraglen} {kmer} --rl {genomes_fofn} --ql {genomes_fofn}  -t {threads} -o {fastani_out}"
    )

    # the output file uses the filenames for the genome names. Replace
    ## with the actual names instead
    with open(fastani_out) as f:
        ani_lines = [x.rstrip() for x in f]

    ani_lines.sort()

    with open(fastani_out, "w") as f:
        for line in ani_lines:
            try:
                genome1, genome2, identity, frags_aln, frags = line.split("\t")
                identity = float(identity)
                frags_aln = int(frags_aln)
                frags = int(frags)
            except:
                raise Exception(f"Error parsing fastANI line:\n{line}")

            genome1_name = genome2name[genome1]
            genome2_name = genome2name[genome2]
            print(
                genome1_name,
                genome2_name,
                identity,
                frags_aln,
                frags,
                name2species[genome1_name],
                name2species[genome2_name],
                sep="\t",
                file=f,
            )

    return fastani_out


def load_distances_tsv(infile):
    distances = {}
    dist_keys = [
        "identity_1",
        "identity_2",
        "aligned_frags_1",
        "frags_1",
        "aligned_frags_2",
        "frags_2",
        "species_1",
        "species_2",
    ]

    with open(infile) as f:
        for line in map(str.rstrip, f):
            fields = line.split("\t")
            # fields are: genome1, genome2, identity, aligned frags, total frags, species1, species2
            if fields[0] == fields[1]:
                continue
            key = tuple(sorted((fields[0], fields[1])))
            key_str = ";".join(key)
            if key_str not in distances:
                distances[key_str] = {x: None for x in dist_keys}
            suffix = "1" if fields[0] == key[0] else "2"
            distances[key_str][f"aligned_frags_{suffix}"] = int(fields[3])
            distances[key_str][f"frags_{suffix}"] = int(fields[4])
            if fields[0] == key[0]:
                distances[key_str]["identity_1"] = float(fields[2])
                distances[key_str]["species_1"] = fields[5]
                distances[key_str]["species_2"] = fields[6]
            else:
                distances[key_str]["identity_2"] = float(fields[2])
                distances[key_str]["species_1"] = fields[6]
                distances[key_str]["species_2"] = fields[5]

    return distances


def analyse_ani_data(ani_data, outprefix):
    results = {}
    empty_results = {
        "within": None,
        "between": {
            "best_match": None,
            "identities": [],
        },
    }

    for genomes_str, ani_dict in ani_data.items():
        genome1, genome2 = genomes_str.split(";")
        assert genome1 != genome2
        identities = [
            v for k, v in ani_dict.items() if k.startswith("identity") and v is not None
        ]
        ani_dict["identity_max"] = max(identities + [0])
        ani_dict["genome1"] = genome1
        ani_dict["genome2"] = genome2
        s1 = ani_dict["species_1"]
        s2 = ani_dict["species_2"]
        for i in "1", "2":
            if ani_dict[f"frags_{i}"] is None:
                ani_dict[f"aligned_frag_%_{i}"] = None
            else:
                ani_dict[f"aligned_frag_%_{i}"] = round(
                    100 * (ani_dict[f"aligned_frags_{i}"] / ani_dict[f"frags_{i}"]), 2
                )

        for s in s1, s2:
            if s not in results:
                results[s] = copy.deepcopy(empty_results)

        if s1 == s2:
            assert results[s1]["within"] is None
            results[s1]["within"] = copy.deepcopy(ani_dict)
            continue

        for s in s1, s2:
            best_match = results[s]["between"]["best_match"]
            if (
                best_match is None
                or best_match["identity_max"] < ani_dict["identity_max"]
            ):
                results[s]["between"]["best_match"] = copy.deepcopy(ani_dict)
            results[s]["between"]["identities"].append(ani_dict["identity_max"])

    scatter_x = []
    scatter_y = []
    hist_within = []
    check_species_lines = []

    for s in results:
        within = results[s]["within"]
        between = results[s]["between"]
        max_within = None if within is None else within["identity_max"]

        if max_within is None or max_within == 0:
            hist_within.append(between["best_match"]["identity_max"])
        else:
            scatter_x.append(max_within)
            scatter_y.append(between["best_match"]["identity_max"])
            if max_within <= between["best_match"]["identity_max"]:
                check_species_lines.extend(
                    [
                        f"CHECK SPECIES {s}, intra species match:",
                        json.dumps(within, indent=2),
                        f"CHECK SPECIES {s}, inter species best match:",
                        json.dumps(between["best_match"], indent=2),
                    ]
                )

    if len(check_species_lines) > 0:
        outfile = f"{outprefix}.info.txt"
        logging.warning(
            f"At least one species where inter-species percent identity is higher than intra-species percent identity. Check the file {outfile}"
        )
        with open(outfile, "w") as f:
            print(*check_species_lines, sep="\n", file=f)

    plot_min = min(scatter_x + scatter_y)
    fig, ax = plt.subplots()
    plt.plot([plot_min, 100], [plot_min, 100], color="grey", linestyle=":")
    ax.scatter(scatter_x, scatter_y, alpha=0.4, s=30)
    ax.set_xlabel("Intra-species percent identity")
    ax.set_ylabel("Inter-species max percent identity")
    ax.set_aspect("equal")
    plt.xlim(plot_min, 100)
    plt.ylim(plot_min, 100)
    plt.savefig(f"{outprefix}.scatter.pdf")
    plt.close()

    fig, ax = plt.subplots()
    plt.hist(hist_within, bins=range(int(plot_min), 101, 1))
    ax.set_title("Inter-species ANI for single-genome species")
    ax.set_ylabel("Number of species")
    ax.set_xlabel("Inter-species max percent identity")
    plt.savefig(f"{outprefix}.hist.pdf")
    plt.close()
