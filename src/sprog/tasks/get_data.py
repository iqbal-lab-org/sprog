import logging

from sprog import sample_set


def run(options):
    samples = sample_set.SampleSet(options.samples_dir, debug=options.debug)
    samples.add_samples_from_tsv(options.manifest_tsv)
    logging.info("Copying/downloading any missing genomes")
    samples.download_missing_genomes()
    logging.info(f"Finished processing. Total samples: {len(samples)}")
