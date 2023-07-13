import os

from sprog import sample_set


def run(options):
    if not os.path.exists(options.samples_dir):
        raise FileNotFoundError(f"Directory not found: {options.samples_dir}")
    samples = sample_set.SampleSet(options.samples_dir, debug=options.debug)
    samples.query_all_unitigs_to_each_species(
        threads=options.threads,
        force=options.force,
        k=options.kmer,
    )
