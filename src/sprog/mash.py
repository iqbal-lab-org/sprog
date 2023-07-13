import copy
import json
import logging
import os

import matplotlib.pyplot as plt

from sprog import utils


def all_v_all_mash(genome2name, name2species, outdir, force=False):
    if os.path.exists(outdir):
        if force:
            utils.syscall(f"rm -rf {outdir}")
        else:
            raise Exception(f"Output directory already exists {outdir}")

    os.mkdir(outdir)

    genomes_fofn = os.path.join(outdir, "genomes.fofn")
    sketch_file = os.path.join(outdir, "mash_sketch")
    distance_tsv = os.path.join(outdir, "distances.tsv")

    with open(genomes_fofn, "w") as f:
        print(*genome2name.keys(), sep="\n", file=f)

    logging.info("Sketching all genomes using mash sketch")
    utils.syscall(
        f"mash sketch -l {genomes_fofn} -o {sketch_file} &> {sketch_file}.stdouterr"
    )
    logging.info("Getting all v all distance matrix using mash dist")
    dist_result = utils.syscall(
        f"mash dist -l {sketch_file}.msh {genomes_fofn}", quiet=True
    )

    # the mash distance tsv uses the filenames for the genome names. Replace
    # with the actual names instead. Add two columns on the end for
    # species name of each genome
    dist_lines = dist_result.stdout.strip().split("\n")
    with open(distance_tsv, "w") as f:
        for line in dist_lines:
            try:
                genome1, genome2, dist, pval, matches = line.split()
            except:
                raise Exception(f"Error parsing mash line:\n{line}")

            genome1_name = genome2name[genome1]
            genome2_name = genome2name[genome2]
            print(
                genome2name[genome1],
                genome2name[genome2],
                dist,
                pval,
                matches,
                name2species[genome1_name],
                name2species[genome2_name],
                sep="\t",
                file=f,
            )

    return distance_tsv


def load_distances_tsv(infile):
    distances = {}
    with open(infile) as f:
        for line in map(str.rstrip, f):
            fields = line.split("\t")
            # fields: ref, qry, distance, p value, matching hashes, spec1, spec2
            if fields[0] == fields[1]:
                continue
            key = tuple(sorted((fields[0], fields[1])))
            key_str = ";".join(key)
            # Sanity check mash results are symmetric
            if key in distances:
                assert distances[key]["distance"] == float(fields[2])
                assert distances[key]["p-value"] == float(fields[3])
                assert distances[key]["matching_hashes"] == fields[4]
            else:
                distances[key_str] = {
                    "distance": float(fields[2]),
                    "p-value": float(fields[3]),
                    "matching_hashes": fields[4],
                    "species1": fields[5],
                    "species2": fields[6],
                }
    return distances


def analyse_distance_data(mash_data, outprefix):
    results = {}
    empty_results = {
        "within": None,
        "between": {
            "best_match": None,
            "distances": [],
        },
    }

    for genomes_str, mash_dict in mash_data.items():
        genome1, genome2 = genomes_str.split(";")
        assert genome1 != genome2
        mash_dict["genome1"] = genome1
        mash_dict["genome2"] = genome2
        s1 = mash_dict["species1"]
        s2 = mash_dict["species2"]

        for s in s1, s2:
            if s not in results:
                results[s] = copy.deepcopy(empty_results)

        if s1 == s2:
            assert results[s1]["within"] is None
            results[s1]["within"] = copy.deepcopy(mash_dict)
            continue

        for s in s1, s2:
            best_match = results[s]["between"]["best_match"]
            if best_match is None or best_match["distance"] > mash_dict["distance"]:
                results[s]["between"]["best_match"] = copy.deepcopy(mash_dict)
            results[s]["between"]["distances"].append(mash_dict["distance"])

    scatter_x = []
    scatter_y = []
    hist_within = []
    check_species_lines = []

    for s in results:
        within = results[s]["within"]
        between = results[s]["between"]
        min_within = None if within is None else within["distance"]

        if min_within is None or min_within == 0:
            hist_within.append(between["best_match"]["distance"])
        else:
            scatter_x.append(min_within)
            scatter_y.append(between["best_match"]["distance"])
            if min_within >= between["best_match"]["distance"]:
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
            f"At least one species where inter-species mash distance less than intra-species mash distance. Check the file {outfile}"
        )
        with open(outfile, "w") as f:
            print(*check_species_lines, sep="\n", file=f)

    plot_max = max(scatter_x + scatter_y)
    fig, ax = plt.subplots()
    plt.plot([0, plot_max], [0, plot_max], color="grey", linestyle=":")
    ax.scatter(scatter_x, scatter_y, alpha=0.4, s=30)
    ax.set_xlabel("Intra-species mash distance")
    ax.set_ylabel("Inter-species min mash distance")
    ax.set_aspect("equal")
    plt.xlim(0, plot_max)
    plt.ylim(0, plot_max)
    plt.savefig(f"{outprefix}.scatter.pdf")
    plt.close()

    fig, ax = plt.subplots()
    plt.hist(hist_within)
    ax.set_title("Inter-species mash distances for single-genome species")
    ax.set_ylabel("Number of species")
    ax.set_xlabel("Inter-species min mash distance")
    plt.savefig(f"{outprefix}.hist.pdf")
    plt.close()
