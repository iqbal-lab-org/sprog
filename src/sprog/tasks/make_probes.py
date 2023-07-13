from sprog import sample_set


def run(options):
    samples = sample_set.SampleSet(options.samples_dir)
    samples.match_probes_to_tree_nodes_2(options.tree_file, force=options.force)
    samples.write_probes(
        options.tree_file,
        mykrobe_tb=options.mykrobe_tb,
        max_probes_per_node=options.probes_per_node,
    )
