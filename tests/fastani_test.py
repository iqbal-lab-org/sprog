import os
import pytest
import random

from sprog import fastani, utils


def test_fastani():
    outdir = "tmp.fastani_test"
    utils.syscall(f"rm -rf {outdir}")
    os.mkdir(outdir)
    random.seed(42)
    fastas = []
    seqs = [random.choices(["A", "C", "G", "T"], k=10000) for _ in range(2)]
    for i in range(0, 10000, 50):
        seqs[0][i] = "A"
    genome1 = "".join(seqs[0]) + "".join(random.choices(["A", "C", "G", "T"], k=3000))
    for i in range(0, 10000, 50):
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

    fastani_outdir = os.path.join(outdir, "FastANI")
    genome2name = {x: x.replace("genome", "species") for x in genomes.values()}
    got_tsv = fastani.all_v_all_fastani(genomes, genome2name, fastani_outdir)
    assert os.path.exists(got_tsv)
    got_distances = fastani.load_distances_tsv(got_tsv)
    assert got_distances == {
        "genome0;genome1": {
            "identity_1": 97.1662,
            "identity_2": 97.1871,
            "frags_1": 4,
            "frags_2": 3,
            "aligned_frags_1": 3,
            "aligned_frags_2": 3,
            "species_1": "species0",
            "species_2": "species1",
        },
    }
    utils.syscall(f"rm -r {outdir}")
