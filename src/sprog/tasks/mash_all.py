import os

from sprog import sample_set


def run(options):
    if not os.path.exists(options.samples_dir):
        raise FileNotFoundError(f"Directory not found: {options.samples_dir}")
    samples = sample_set.SampleSet(options.samples_dir, debug=options.debug)
    samples.mash_all_v_all(force=options.force)
    samples.fastani_mash_scatter()
