import filecmp
import pytest
import os

from sprog import gtdb, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data")


def test_taxon_to_species():
    assert gtdb.taxon_to_species("g__genus;s__genus species") == "genus species"
    assert gtdb.taxon_to_species("g__genus;s__") is None


def test_strip_ncbi_name():
    assert gtdb.strip_ncbi_name("gen spec") == "gen spec"
    assert gtdb.strip_ncbi_name("gen spec blah") == "gen spec"
    assert gtdb.strip_ncbi_name("gen spec subsp. foo") == "gen spec subsp. foo"
    assert gtdb.strip_ncbi_name("gen spec subsp. foo blah") == "gen spec subsp. foo"


def test_make_mykrobe_files():
    gtdb_tsv = os.path.join(data_dir, "gtdb.make_mykrobe_files.tsv")
    outdir = "tmp.gtdb.make_mykrobe_files.out"
    utils.syscall(f"rm -rf {outdir}")
    gtdb.make_mykrobe_files(outdir, gtdb_tsv=gtdb_tsv, genomes_per_species=2)
    expect_manifest = os.path.join(
        data_dir, "gtdb.make_mykrobe_files.expect_manifest.tsv"
    )
    expect_tree = os.path.join(data_dir, "gtdb.make_mykrobe_files.expect_tree.tsv")
    got_manifest = os.path.join(outdir, "manifest.tsv")
    got_tree = os.path.join(outdir, "tree.tsv")
    assert filecmp.cmp(got_manifest, expect_manifest, shallow=False)
    assert filecmp.cmp(got_tree, expect_tree, shallow=False)

    expect_ncbi_names = {
        "Mycobacterium_chelonae_subsp._chelonae": {
            "Mycobacterium_chelonae_subsp._chelonae": 1
        },
        "Mycobacterium_chelonae_subsp._gwanakae": {
            "Mycobacterium_chelonae_subsp._gwanakae": 1
        },
        "Mycobacterium_tuberculosis": {"Mycobacterium_tuberculosis": 3},
        "Mycobacterium_tuberculosis_variant_africanum": {
            "Mycobacterium_tuberculosis_variant_africanum": 1
        },
        "Mycobacterium_tuberculosis_variant_bovis": {
            "Mycobacterium_tuberculosis_variant_bovis": 1
        },
        "Mycobacterium_xenopi": {"Mycobacterium_xenopi": 3},
    }
    got_ncbi_names = utils.load_json(os.path.join(outdir, "ncbi_names.json"))
    assert got_ncbi_names == expect_ncbi_names
    utils.syscall(f"rm -r {outdir}")
