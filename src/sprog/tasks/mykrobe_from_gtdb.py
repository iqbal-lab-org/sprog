from sprog import gtdb


def run(options):
    gtdb.make_mykrobe_files(
        options.outdir, genomes_per_species=options.genomes_per_species
    )
