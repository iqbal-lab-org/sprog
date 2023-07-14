import csv
import itertools
import os
import pytest
import random

import bitarray
import pyfastaq

from sprog import sample_set, tree, utils


def test_sample_set():
    # --------------- make sample set with toy genomes -------------------------
    outdir = "tmp.sample_set"
    utils.syscall(f"rm -rf {outdir}")
    os.mkdir(outdir)
    random.seed(42)
    fastas = []
    seqs = ["".join(random.choices(["A", "C", "G", "T"], k=200)) for _ in range(8)]
    genomes = [
        seqs[0] + seqs[1] + seqs[2] + seqs[6],
        seqs[3] + seqs[1] + seqs[4] + seqs[6],
        seqs[5] + seqs[6],
        seqs[0] + seqs[7],
    ]

    for i, genome in enumerate(genomes):
        fa = os.path.join(outdir, f"ref.{i}.fa")
        fastas.append(fa)
        with open(fa, "w") as f:
            print(f">genome.{i}", file=f)
            print(genome, file=f)

    meta_fields = [
        ["sample", "species", "data_source", "accession", "file"],
        ["sample1", "species1", "copy_file", ".", fastas[0]],
        ["sample2", "species1", "symlink", ".", fastas[1]],
        ["sample3", "species2", "copy_file", ".", fastas[2]],
        ["sample4", "species3", "copy_file", ".", fastas[3]],
    ]

    samples_tsv = os.path.join(outdir, "samples.tsv")
    with open(samples_tsv, "w") as f:
        for fields in meta_fields:
            print(*fields, sep="\t", file=f)

    sample_dir = os.path.join(outdir, "Sample_set")
    samples = sample_set.SampleSet(sample_dir, debug=True)
    samples.add_samples_from_tsv(samples_tsv)
    samples.download_missing_genomes()

    # --------------- mash all v all genomes -----------------------------------
    samples.mash_all_v_all()
    # basic check of the resulting distances. Not sure fow reproducible exact
    # mash distances and p-values are, so basic change they are a wide acceptable
    # range
    assert len(samples.mash_distances) == 6
    sample_names = ["sample1", "sample2", "sample3", "sample4"]
    for (s1, s2) in itertools.combinations(sample_names, 2):
        key = ";".join(sorted([s1, s2]))
        assert key in samples.mash_distances
        assert 0 <= samples.mash_distances[key]["p-value"] <= 1
        assert 0 <= samples.mash_distances[key]["distance"] <= 1
    assert os.path.exists(os.path.join(samples.mash_dir, "analysis.hist.pdf"))
    assert os.path.exists(os.path.join(samples.mash_dir, "analysis.scatter.pdf"))
    # should do nothing because fastani not run yet
    samples.fastani_mash_scatter()
    assert not os.path.exists(samples.mash_fastani_dir)

    # --------------- fastani all v all genomes -----------------------------------
    fastani_distances = samples.fastani_all_v_all(fraglen=100, kmer=7)
    # basic check of the resulting distances. Not sure fow reproducible exact
    # fastani distances and p-values are, so basic change they are a wide acceptable
    # range. Note samples 2 and 4 share nothing in common, and also
    # 3 and 4 share nothing in common. So those pairs will not be in the
    # output. All other combinations of samples 1, 2, 3, 4 will be, so expect
    # 4 results not 6.
    assert len(fastani_distances) == 4
    sample_names = ["sample1", "sample2", "sample3", "sample4"]
    for key_tuple in itertools.combinations(sample_names, 2):
        key = ";".join(sorted(key_tuple))
        if key_tuple == ("sample2", "sample4") or key_tuple == ("sample3", "sample4"):
            assert key not in fastani_distances
        else:
            assert key in fastani_distances
            d = fastani_distances[key]
            assert 99.0 <= d["identity_1"] <= 100.0
            assert 99.0 <= d["identity_2"] <= 100.0
            assert 2 <= d["aligned_frags_1"] <= d["frags_1"] <= 8
            assert 2 <= d["aligned_frags_2"] <= d["frags_2"] <= 8
    assert os.path.exists(samples.fastani_dir)
    assert os.path.exists(samples.fastani_json)
    assert os.path.exists(os.path.join(samples.fastani_dir, "analysis.hist.pdf"))
    assert os.path.exists(os.path.join(samples.fastani_dir, "analysis.scatter.pdf"))

    # --------------- scatter fastani vs mash ----------------------------------
    assert not os.path.exists(samples.mash_fastani_dir)
    samples.fastani_mash_scatter()
    assert os.path.exists(samples.mash_fastani_dir)
    assert os.path.exists(os.path.join(samples.mash_fastani_dir, "scatter.pdf"))
    assert os.path.exists(os.path.join(samples.mash_fastani_dir, "scatter.2.pdf"))

    # --------------- make unitigs from all the genomes ------------------------
    assert not os.path.exists(samples.all_unitigs_fa)
    samples.make_unitigs_one_graph_all_samples()
    assert os.path.exists(samples.all_unitigs_fa)
    reader = pyfastaq.sequences.file_reader(samples.all_unitigs_fa)
    got_seq_count = len([_ for _ in reader])
    assert got_seq_count == 10

    # --------- make and then query each sample unitigs vs all unitigs ---------
    samples.query_all_unitigs_to_each_species(threads=2)
    expect_json = samples.all_query_json
    assert os.path.exists(expect_json)
    got_json_data = utils.load_json(expect_json)
    assert len(got_json_data) == 2
    assert got_json_data["number_of_unitigs"] == 10
    assert len(got_json_data["filenames"]) == 3
    assert "species1" in got_json_data["filenames"]
    assert "species2" in got_json_data["filenames"]
    for species, filename in got_json_data["filenames"].items():
        assert os.path.exists(os.path.join(samples.all_query_dir, filename))

    # ------------- presence absence matrix ------------------------------------
    got_presence = samples.load_species_query_files()
    # we don't know the order of the unitigs, which means don't know exact
    # contents of the bitarrays. The order of 0s and 1s depend on unitig order.
    # Something like this:
    unitig_presence = {
        "species1": bitarray.bitarray("1111111001"),
        "species2": bitarray.bitarray("0000100101"),
        "species3": bitarray.bitarray("1000000011"),
    }
    # ... which we'll use later to test matching unitigs to nodes
    assert len(got_presence) == 3
    assert list(got_presence.keys()) == ["species1", "species2", "species3"]

    # ------------- match unitigs to tree nodes --------------------------------
    tmp_newick = os.path.join(outdir, "tree.newick")
    with open(tmp_newick, "w") as f:
        print("((species1,species2)nodeX,species3);", file=f)
    t = tree.Tree(tmp_newick)
    got_probe_matches = samples.match_probes_to_tree_nodes(
        t, unitig_presence=unitig_presence
    )
    expect = {
        "species1": {1, 2, 3, 5, 6},
        "species2": {7},
        "species3": {8},
        "nodeX": {4},
        "root": {9},
    }
    got_by_name = {k.name: v for k, v in got_probe_matches.items()}
    assert got_by_name == expect

    tmp_tree_out = os.path.join(outdir, "tree.txt")
    t.write_to_file(tmp_tree_out, node_unitig_combos=got_probe_matches)
    assert os.path.exists(tmp_tree_out)

    # We've got two different methods to match probes. Test the second method
    samples.match_probes_to_tree_nodes_2(tmp_newick)
    assert os.path.exists(samples.unitig_counts_tree)
    expect = {
        "species1": {2},
        "species2": {7},
        "species3": {0, 8},
        "nodeX": {4},
    }
    with open(samples.node2unitig_tsv) as f:
        for d in csv.DictReader(f, delimiter="\t"):
            matches = set(int(x) for x in d["Unitigs"].split(","))
            # the unitig numbers are not deterministic, so check we get
            # the right number of them
            assert len(matches) == len(expect[d["Node"]])

    # ---------------- simulate perfect reads ----------------------------------
    reads_outdir = os.path.join(outdir, "sim_reads")
    samples.simulate_perfect_reads(reads_outdir, length=30, frag_length=100, depth=1)
    for sample_name in samples.samples:
        p = os.path.join(reads_outdir, sample_name)
        assert os.path.exists(f"{p}_1.fq.gz")
        assert os.path.exists(f"{p}_2.fq.gz")
    utils.syscall(f"rm -r {outdir}")
