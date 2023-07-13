import os
import pytest
import random

from sprog import mash, utils


def test_mash():
    outdir = "tmp.mash_test"
    utils.syscall(f"rm -rf {outdir}")
    os.mkdir(outdir)
    random.seed(42)
    fastas = []
    seqs = [random.choices(["A", "C", "G", "T"], k=500) for _ in range(2)]
    for i in range(0, 500, 50):
        seqs[0][i] = "A"
    genome1 = "".join(seqs[0])
    for i in range(0, 500, 50):
        seqs[0][i] = "T"
    genome2 = "".join(seqs[0])
    fastas = [os.path.join(outdir, f"g{i}.fa") for i in range(1, 4, 1)]
    genomes = {fastas[i]: f"genome{i}" for i in range(len(fastas))}
    with open(fastas[0], "w") as f:
        print(">g1", genome1, sep="\n", file=f)
    with open(fastas[1], "w") as f:
        print(">g2", genome2, sep="\n", file=f)
    with open(fastas[2], "w") as f:
        print(">g3", "".join(seqs[1]), sep="\n", file=f)
    genome2name = {x: x.replace("genome", "species") for x in genomes.values()}

    mash_outdir = os.path.join(outdir, "Mash")
    got_tsv = mash.all_v_all_mash(genomes, genome2name, mash_outdir)
    assert os.path.exists(got_tsv)
    got_distances = mash.load_distances_tsv(got_tsv)
    assert got_distances == {
        "genome0;genome1": {
            "distance": 0.0239955,
            "p-value": 0.0,
            "matching_hashes": "290/670",
            "species1": "species0",
            "species2": "species1",
        },
        "genome0;genome2": {
            "distance": 1.0,
            "p-value": 1.0,
            "matching_hashes": "0/960",
            "species1": "species0",
            "species2": "species2",
        },
        "genome1;genome2": {
            "distance": 1.0,
            "p-value": 1.0,
            "matching_hashes": "0/960",
            "species1": "species1",
            "species2": "species2",
        },
    }
    utils.syscall(f"rm -r {outdir}")
