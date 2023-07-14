#!/usr/bin/env python3

import argparse
import logging

import sprog


def main(args=None):
    parser = argparse.ArgumentParser(
        prog="sprog",
        usage="sprog <command> <options>",
        description="sprog: species probe generator",
    )
    parser.add_argument("--version", action="version", version=sprog.__version__)

    subparsers = parser.add_subparsers(title="Available commands", help="", metavar="")

    # ----------- general options common to all tasks ------------------------
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument(
        "--debug",
        action="store_true",
        help="More verbose logging, and less file cleaning",
    )
    samples_dir_parser = argparse.ArgumentParser(add_help=False)
    samples_dir_parser.add_argument(
        "samples_dir",
        help="Root directory for data",
    )

    force_parser = argparse.ArgumentParser(add_help=False)
    force_parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing files from previous run (if they exist)",
    )
    common_other = argparse.ArgumentParser(add_help=False)
    common_other.add_argument(
        "--threads",
        type=int,
        help="Number of threads [%(default)s]",
        default=1,
        metavar="INT",
    )
    kmer_parser = argparse.ArgumentParser(add_help=False)
    kmer_parser.add_argument(
        "--kmer",
        "-k",
        help="kmer length used by Bifrost [%(default)s]",
        type=int,
        metavar="INT",
        default=31,
    )

    # ----------------------------- get_data ----------------------------------
    subparser_get_data = subparsers.add_parser(
        "get_data",
        parents=[common_parser, samples_dir_parser],
        help="Get genomes etc",
        usage="sprog get_data [options] <samples_dir> <manifest_tsv>",
        description="Get genomes. Create new dir of samples, or add to existing dir",
    )
    subparser_get_data.add_argument(
        "manifest_tsv",
        help="TSV file of sample information",
    )
    subparser_get_data.set_defaults(func=sprog.tasks.get_data.run)

    # ------------------------------ mash_all ---------------------------------
    subparser_mash_all = subparsers.add_parser(
        "mash_all",
        parents=[common_parser, samples_dir_parser, force_parser],
        help="Run all vs all genomes mash",
        usage="sprog mash_all [options] <samples_dir>",
        description="Run mash of all vs all genomes",
    )
    subparser_mash_all.set_defaults(func=sprog.tasks.mash_all.run)

    # ---------------------------- fastani_all --------------------------------
    subparser_fastani_all = subparsers.add_parser(
        "fastani_all",
        parents=[common_parser, samples_dir_parser, common_other, force_parser],
        help="Run all vs all genomes fastANI",
        usage="sprog fastani_all [options] <samples_dir>",
        description="Run fastANI of all vs all genomes",
    )
    subparser_fastani_all.set_defaults(func=sprog.tasks.fastani_all.run)

    # ----------------------- build_all_samples_unitigs -----------------------
    subparser_build_all_samples_unitigs = subparsers.add_parser(
        "build_all_samples_unitigs",
        parents=[
            common_parser,
            samples_dir_parser,
            force_parser,
            kmer_parser,
            common_other,
        ],
        help="Build unitigs from all samples",
        usage="sprog get_data [options] <samples_dir>",
        description="Build unitigs from all samples",
    )
    subparser_build_all_samples_unitigs.set_defaults(
        func=sprog.tasks.build_all_samples_unitigs.run
    )

    # --------------------- per_sample_presence -------------------------------
    subparser_per_sample_presence = subparsers.add_parser(
        "per_sample_presence",
        parents=[
            common_parser,
            samples_dir_parser,
            force_parser,
            kmer_parser,
            common_other,
        ],
        help="Match all unitigs to each sample",
        usage="sprog get_data [options] <samples_dir>",
        description="Match all unitigs to each sample. Assumes `sprog build_all_samples_unitigs` has been run",
    )
    subparser_per_sample_presence.set_defaults(func=sprog.tasks.per_sample_presence.run)

    # ---------------------- make_probes --------------------------------------
    subparser_make_probes = subparsers.add_parser(
        "make_probes",
        parents=[common_parser, samples_dir_parser, force_parser],
        help="Make probes",
        usage="sprog get_data [options] <samples_dir>",
        description="Make probes. Assumes per_sample_presence has been run",
    )
    subparser_make_probes.add_argument(
        "--mykrobe_tb",
        action="store_true",
        help="Write files for mykrobe TB speciation",
    )
    subparser_make_probes.add_argument(
        "--probes_per_node",
        type=int,
        help="Target number of probes per species (or per node in the tree) [%(default)s]",
        default=2000,
        metavar="INT",
    )
    subparser_make_probes.add_argument(
        "tree_file",
        help="Filename of tree",
    )
    subparser_make_probes.set_defaults(func=sprog.tasks.make_probes.run)

    # ---------------------- mykrobe_from_gtdb --------------------------------
    subparser_mykrobe_from_gtdb = subparsers.add_parser(
        "mykrobe_from_gtdb",
        parents=[common_parser],
        help="Generate sprog input files to make mykrobe probes, using data from GTDB",
        usage="sprog mykrobe_from_gtdb [options] <outdir>",
        description="Generate sprog input files to make mykrobe probes, using data from GTDB",
    )
    subparser_mykrobe_from_gtdb.add_argument("outdir", help="Output directory")
    subparser_mykrobe_from_gtdb.add_argument(
        "--genomes_per_species",
        type=int,
        default=2,
        help="Genoems to use per species [%(default)s]",
        metavar="INT",
    )
    subparser_mykrobe_from_gtdb.set_defaults(func=sprog.tasks.mykrobe_from_gtdb.run)

    # ---------------------- sim_perfect_reads --------------------------------
    subparser_sim_perfect_reads = subparsers.add_parser(
        "sim_perfect_reads",
        parents=[common_parser, samples_dir_parser, force_parser],
        help="Simulate Illumina reads with zero errors from each genome",
        usage="sprog sim_perfect_reads [options] <outdir>",
        description="Simulate Illumina reads with zero errors from each genome",
    )
    subparser_sim_perfect_reads.add_argument("outdir", help="Output directory")
    subparser_sim_perfect_reads.add_argument(
        "--read_length",
        type=int,
        help="Read length [%(default)s]",
        default=150,
        metavar="INT",
    )
    subparser_sim_perfect_reads.add_argument(
        "--frag_length",
        type=int,
        help="Fragment length [%(default)s]",
        default=400,
        metavar="INT",
    )
    subparser_sim_perfect_reads.add_argument(
        "--depth",
        type=int,
        help="Mean read depth [%(default)s]",
        default=10,
        metavar="INT",
    )
    subparser_sim_perfect_reads.set_defaults(func=sprog.tasks.sim_perfect_reads.run)

    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        return

    logging.basicConfig(
        format="[%(asctime)s sprog %(levelname)s] %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S%z",
    )
    log = logging.getLogger()
    if args.debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
