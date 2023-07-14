import os

from sprog import sample_set


def run(options):
    if not os.path.exists(options.samples_dir):
        raise FileNotFoundError(f"Directory not found: {options.samples_dir}")
    samples = sample_set.SampleSet(options.samples_dir, debug=options.debug)
    samples.simulate_perfect_reads(
        options.outdir,
        length=options.read_length,
        frag_length=options.frag_length,
        depth=options.depth,
        force=options.force,
    )
