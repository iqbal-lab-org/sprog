import os

from sprog import sample_set


def run(options):
    if not os.path.exists(options.samples_dir):
        raise FileNotFoundError(f"Directory not found: {options.samples_dir}")
    samples = sample_set.SampleSet(options.samples_dir, debug=options.debug)
    samples.make_unitigs_one_graph_all_samples(
        force=options.force,
        threads=options.threads,
        k=options.kmer,
    )
