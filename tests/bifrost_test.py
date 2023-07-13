import gzip
import pytest
import os
import random

from sprog import bifrost, utils


def looks_like_fa(filename):
    with gzip.open(filename, "rt") as f:
        line = f.readline()
        return line.startswith(">")


def test_build_and_query():
    # make test files
    outprefix = "tmp.test_bifrost_build"
    utils.syscall(f"rm -f {outprefix}*")
    random.seed(42)
    test_fa1 = "tmp.test_bifrost_build.1.fa"
    test_fa2 = "tmp.test_bifrost_build.2.fa"
    seq1 = "".join(random.choices(["A", "C", "G", "T"], k=100))
    seq2 = "".join(random.choices(["A", "C", "G", "T"], k=100))
    with open(test_fa1, "w") as f:
        print(">seq1", seq1, sep="\n", file=f)
    with open(test_fa2, "w") as f:
        print(">seq2", seq2, sep="\n", file=f)

    query_fa = f"{outprefix}.query.fa"
    with open(query_fa, "w") as f:
        print(">0", "A" * 100, sep="\n", file=f)
        print(">1", seq1, sep="\n", file=f)
        print(">2", seq2, sep="\n", file=f)

    # test build bfg
    got_bfg = bifrost.build(test_fa1, outprefix, "bfg")
    expect_fa = f"{outprefix}.fasta.gz"
    expect_gfa = f"{outprefix}.gfa.gz"
    expect_bfg = f"{outprefix}.bfg"
    expect_bfi = f"{outprefix}.bfi"
    assert got_bfg == expect_bfg
    assert not os.path.exists(expect_fa)
    assert not os.path.exists(expect_gfa)
    assert os.path.exists(expect_bfg)
    assert os.path.exists(expect_bfi)

    # test querying the bfg we just made
    qry_out = bifrost.query(expect_bfg, query_fa, f"{outprefix}.query")
    assert os.path.exists(qry_out)
    with open(qry_out) as f:
        lines = [x.rstrip() for x in f]
    assert lines == [
        "query_name\tpresence_query",
        "0\t0",
        "1\t1",
        "2\t0",
    ]
    os.unlink(qry_out)
    os.unlink(expect_bfg)
    os.unlink(expect_bfi)

    # test build fasta unitigs
    filenames_file = f"{outprefix}.filenames"
    with open(filenames_file, "w") as f:
        print(test_fa1, test_fa2, sep="\n", file=f)
    got = bifrost.build(filenames_file, outprefix, "fasta")
    assert got == expect_fa
    assert os.path.exists(expect_fa)
    assert not os.path.exists(expect_gfa)
    assert not os.path.exists(expect_bfg)
    assert not os.path.exists(expect_bfi)
    assert looks_like_fa(expect_fa)
    os.unlink(expect_fa)

    b_and_q_out = f"{outprefix}.b_and_q"
    got_array_length, qry_out = bifrost.build_and_query(
        [test_fa1, test_fa2],
        query_fa,
        b_and_q_out,
        combine_method="union",
    )
    assert got_array_length == 3
    assert os.path.exists(qry_out)
    got_array = bifrost.load_bitarray_file(qry_out, got_array_length)
    assert list(got_array) == [0, 1, 1]

    got_array_length, qry_out = bifrost.build_and_query(
        [test_fa1, test_fa2],
        query_fa,
        b_and_q_out,
        combine_method="intersection",
    )
    assert got_array_length == 3
    assert os.path.exists(qry_out)
    got_array = bifrost.load_bitarray_file(qry_out, got_array_length)
    assert list(got_array) == [0, 0, 0]
    utils.syscall(f"rm {outprefix}*")


def test_query_tsv_to_bitarray_file_and_then_load():
    infile = "tmp.bifrost.query_tsv_to_bitarray_file.tsv"
    outfile = "tmp.bifrost.query_tsv_to_bitarrayt_file.txt"
    utils.syscall(f"rm -f {infile} {outfile}")
    expect = [0, 1, 0, 1, 1]
    with open(infile, "w") as f:
        print("query_name\tpresence_query", file=f)
        for i, bit in enumerate(expect):
            print(i, bit, sep="\t", file=f)
    got_length = bifrost.query_tsv_to_bitarray_file(infile, outfile)
    assert got_length == len(expect)
    got_array = bifrost.load_bitarray_file(outfile, len(expect))
    assert list(got_array) == expect
    os.unlink(infile)
    os.unlink(outfile)
