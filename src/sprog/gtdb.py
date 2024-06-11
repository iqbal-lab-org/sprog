import csv
from collections import Counter
import gzip
import logging
from operator import itemgetter
import os
import random
import re

from sprog import utils


MYCO_SP_REGEX = re.compile("Mycobacterium sp[0-9]+")
GTDB_VERSION = "220"
GTDB_TSV_GZ = f"bac120_metadata_r{GTDB_VERSION}.tsv.gz"
GTDB_TSV_URL = f"https://data.gtdb.ecogenomic.org/releases/release{GTDB_VERSION}/{GTDB_VERSION}.0/{GTDB_TSV_GZ}"
MTB_VARIANTS = {"africanum", "caprae", "microti", "pinnipedii"}


def get_gtdb_tsv(outdir):
    logging.info(f"Download GTDB TSV file release version {GTDB_VERSION}")
    utils.syscall(f"wget -q {GTDB_TSV_URL}", cwd=outdir)
    logging.info(f"Finished getting GTDB TSV file {GTDB_TSV_GZ}")
    return os.path.join(outdir, GTDB_TSV_GZ)


def taxon_to_species(taxon):
    fields = [x.strip() for x in taxon.split(";")]
    taxons = {}
    for f in fields:
        assert "__" in f
        level, value = f.split("__")
        assert level not in taxons
        taxons[level] = value

    if taxons["s"] == "":
        return None
    else:
        return taxons.get("s", None)


def load_gtdb_tsv(gtdb_tsv):
    species = {}
    logging.info(f"Getting Mycobaterium species from GTDB file {gtdb_tsv}")
    wanted_keys = [
        "accession",
        "gtdb_representative",
        "gtdb_type_designation_ncbi_taxa",
        "ncbi_genbank_assembly_accession",
        "ncbi_organism_name",
    ]
    fopen = gzip.open if gtdb_tsv.endswith(".gz") else open
    with fopen(gtdb_tsv, "rt") as f:
        for count, d in enumerate(csv.DictReader(f, delimiter="\t")):
            if count > 0 and count % 100_000 == 0:
                logging.info(f" ... parsed {count} lines of TSV file")
            this_species = taxon_to_species(d["gtdb_taxonomy"])
            if not this_species.startswith("Mycobacterium"):
                continue
            if MYCO_SP_REGEX.match(this_species) is not None:
                continue
            if this_species not in species:
                species[this_species] = {
                    "gtdb_taxon": d["gtdb_taxonomy"],
                    "samples": [],
                }
            assert species[this_species]["gtdb_taxon"] == d["gtdb_taxonomy"]
            species[this_species]["samples"].append({k: d[k] for k in wanted_keys})

    logging.info(f"Got {len(species)} species from GTDB TSV file")
    return species


def strip_ncbi_name(name):
    name = name.replace("[", "").replace("]", "")
    fields = name.split()
    if len(fields) == 2:
        return name
    elif "subsp." in fields or "variant" in fields:
        return " ".join(fields[:4])
    else:
        return " ".join(fields[:2])


def gtdb_dict_to_manifest_line(d, species):
    return "\t".join(
        [
            d["ncbi_genbank_assembly_accession"],
            species,
            "download",
            d["ncbi_genbank_assembly_accession"],
            ".",
        ]
    )


def gtdb_dict_list_to_sampled_manifest_lines(gtdb_dicts, species, to_sample):
    if len(gtdb_dicts) == 0:
        return None
    lines = []
    random.seed(42)
    for d in random.sample(gtdb_dicts, min(len(gtdb_dicts), to_sample)):
        lines.append(gtdb_dict_to_manifest_line(d, species))
    return lines


def species_dict_to_manifest_lines(gtdb_species, gtdb_d, genomes_per_species=2):
    lines_out = []
    species_aliases = {}
    samples = gtdb_d["samples"]
    gtdb_reps = [x for x in samples if x["gtdb_representative"] == "t"]
    assert len(gtdb_reps) == 1
    gtdb_rep = gtdb_reps[0]
    gb_accessions = {}

    type_strains = [
        x
        for x in samples
        if x["gtdb_type_designation_ncbi_taxa"].startswith("type strain")
    ]
    type_strain_names = {strip_ncbi_name(x["ncbi_organism_name"]) for x in type_strains}
    is_mtb = gtdb_species == "Mycobacterium tuberculosis"
    has_subspecies = is_mtb or any("subsp" in x for x in type_strain_names)
    tree_lines = set()
    sub_ntm_complex = "subNon_tuberculosis_mycobacterium_complex"
    sub_mtb_complex = "subMycobacterium_tuberculosis_complex"

    if has_subspecies:
        subspecies = {}
        for d in samples:
            if (
                d["ncbi_organism_name"]
                == "Mycobacterium tuberculosis subsp. tuberculosis"
            ):
                continue
            if "subsp" in d["ncbi_organism_name"] or is_mtb:
                name = strip_ncbi_name(d["ncbi_organism_name"])
                if name not in subspecies:
                    subspecies[name] = {"type": [], "non_type": [], "ncbi_names": {}}
                if d["gtdb_type_designation_ncbi_taxa"].startswith("type strain"):
                    subspecies[name]["type"].append(d)
                else:
                    subspecies[name]["non_type"].append(d)

                ncbi_name = strip_ncbi_name(d["ncbi_organism_name"]).replace(" ", "_")
                if not ncbi_name.endswith(" sp."):
                    subspecies[name]["ncbi_names"][ncbi_name] = (
                        subspecies[name]["ncbi_names"].get(ncbi_name, 0) + 1
                    )

        for name, d in sorted(subspecies.items()):
            matching_gtdb_rep = [
                x for x in d["type"] if x["accession"] == gtdb_rep["accession"]
            ]
            if len(matching_gtdb_rep) == 1:
                record_to_use = gtdb_rep
            elif len(d["type"]) == 1:
                record_to_use = d["type"][0]
            elif len(d["type"]) > 1:
                d["type"].sort(key=itemgetter("ncbi_genbank_assembly_accession"))
                record_to_use = d["type"][0]
            else:
                assert len(d["non_type"]) > 0
                d["non_type"].sort(key=itemgetter("ncbi_genbank_assembly_accession"))
                record_to_use = d["non_type"][0]

            subspecies_name = name.split()[-1]

            if gtdb_species.endswith(" tuberculosis"):
                if subspecies_name in MTB_VARIANTS:
                    species = f"Mycobacterium tuberculosis variant {subspecies_name}"
                else:
                    species = f"Mycobacterium {subspecies_name}"
                tree_lines.add(f"{sub_mtb_complex}\t{species}")
            else:
                species = f"{gtdb_species} subsp. {subspecies_name}"
                tree_lines.add(f"{sub_ntm_complex}\t{gtdb_species}\t{species}")

            gb_accessions[species] = {
                "for_probes": [record_to_use["ncbi_genbank_assembly_accession"]],
                "non_probes": [],
            }
            lines_out.append(gtdb_dict_to_manifest_line(record_to_use, species))
            species_aliases[name.replace(" ", "_")] = subspecies[name]["ncbi_names"]

            if genomes_per_species > 1:
                other_records = [
                    x for x in d["type"] + d["non_type"] if x != record_to_use
                ]
                assert len(other_records) == len(d["type"]) + len(d["non_type"]) - 1
                new_lines = gtdb_dict_list_to_sampled_manifest_lines(
                    other_records, species, genomes_per_species - 1
                )
                if new_lines is not None:
                    lines_out.extend(new_lines)
                    gb_accessions[species]["for_probes"].extend(
                        [x.split("\t")[0] for x in new_lines]
                    )

            gb_accessions[species]["non_probes"] = [
                x["ncbi_genbank_assembly_accession"]
                for x in d["type"] + d["non_type"]
                if x["ncbi_genbank_assembly_accession"]
                not in gb_accessions[species]["for_probes"]
            ]

    else:
        ncbi_names = [
            strip_ncbi_name(x["ncbi_organism_name"]).replace(" ", "_") for x in samples
        ]
        species_aliases = {
            gtdb_species.replace(" ", "_"): Counter(
                x.replace(" ", "_") for x in ncbi_names if not x.endswith(" sp.")
            )
        }
        lines_out.append(gtdb_dict_to_manifest_line(gtdb_rep, gtdb_species))
        gb_accessions[gtdb_species] = {"for_probes": [], "non_probes": []}
        gb_accessions[gtdb_species]["for_probes"].append(
            gtdb_rep["ncbi_genbank_assembly_accession"]
        )
        if genomes_per_species > 1:
            non_reps = [x for x in samples if x["gtdb_representative"] != "t"]
            new_lines = gtdb_dict_list_to_sampled_manifest_lines(
                non_reps, gtdb_species, genomes_per_species - 1
            )
            if new_lines is not None:
                lines_out.extend(new_lines)
                gb_accessions[gtdb_species]["for_probes"].extend(
                    [x.split("\t")[0] for x in new_lines]
                )
        tree_lines.add(f"{sub_ntm_complex}\t{gtdb_species}")
        gb_accessions[gtdb_species]["non_probes"] = [
            x["ncbi_genbank_assembly_accession"]
            for x in samples
            if x["ncbi_genbank_assembly_accession"]
            not in gb_accessions[gtdb_species]["for_probes"]
        ]

    return lines_out, tree_lines, species_aliases, gb_accessions


def write_manifest_and_tree_tsv_and_ncbi_json(
    species, manifest_tsv, tree_tsv, ncbi_json, gb_json, genomes_per_species=2
):
    tree_lines = set()
    species_aliases = {}
    genbank_accessions = {}

    with open(manifest_tsv, "w") as f:
        print("sample", "species", "data_source", "accession", "file", sep="\t", file=f)
        for gtdb_species, gtdb_d in sorted(species.items()):
            (
                tsv_lines,
                new_tree_lines,
                aliases,
                new_gb_acc,
            ) = species_dict_to_manifest_lines(
                gtdb_species, gtdb_d, genomes_per_species=genomes_per_species
            )
            if len(tsv_lines) > 0:
                print(*tsv_lines, sep="\n", file=f)
                species_aliases.update(aliases)
                genbank_accessions.update(new_gb_acc)

            tree_lines.update(new_tree_lines)

    with open(tree_tsv, "w") as f:
        tree_lines = sorted(list(tree_lines))
        print(*tree_lines, sep="\n", file=f)

    utils.write_json(species_aliases, ncbi_json)
    utils.write_json(genbank_accessions, gb_json)


def make_mykrobe_files(outdir, gtdb_tsv=None, genomes_per_species=2):
    os.mkdir(outdir)
    #gtdb_tsv = f"{outdir}/{GTDB_TSV_GZ}"
    if gtdb_tsv is None:
        gtdb_tsv = get_gtdb_tsv(outdir)
    species = load_gtdb_tsv(gtdb_tsv)
    logging.info("Generating manifest and tree TSV files from GTDB data")
    manifest_tsv = os.path.join(outdir, "manifest.tsv")
    tree_tsv = os.path.join(outdir, "tree.tsv")
    ncbi_json = os.path.join(outdir, "ncbi_names.json")
    gb_json = os.path.join(outdir, "genbank_accessions.json")
    write_manifest_and_tree_tsv_and_ncbi_json(
        species,
        manifest_tsv,
        tree_tsv,
        ncbi_json,
        gb_json,
        genomes_per_species=genomes_per_species,
    )
    logging.info(f"Manifest file for input to get_data: {manifest_tsv}")
    logging.info(f"Tree file for input to make_probes: {tree_tsv}")
    logging.info(f"Species lookup to get NCBI names for mykrobe panel: {ncbi_json}")
    logging.info("All done")
