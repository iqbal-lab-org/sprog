import bitarray
import csv
import os

from sprog import utils


def build(ref_file, outprefix, build_type, threads=1, k=None, debug=False):
    if k is None:
        k = 31

    if build_type == "fasta":
        extra_opts = "--fasta-out --no-index-out"
        outfile = f"{outprefix}.fasta.gz"
    elif build_type == "bfg":
        extra_opts = "--bfg-out"
        outfile = f"{outprefix}.bfg"
    elif build_type == "color":
        extra_opts = "--colors --no-index-out"
        outfile = f"{outprefix}.gfa.gz"
    else:
        raise NotImplementedError(
            f"build_type must be one of 'fasta', 'bfg', 'color', not '{build_type}'"
        )

    command = f"Bifrost build {extra_opts} -v -r {ref_file} -k {k} -t {threads} -o {outprefix}"
    utils.syscall(command, quiet=not debug)
    return outfile


def query(
    ref_file, query_file, outprefix, color_bfg=None, kmer_ratio=None, debug=False
):
    if kmer_ratio is None:
        extra_opts = ""
    else:
        extra_opts = f"--kmer-ratio {kmer_ratio}"

    if color_bfg is not None:
        extra_opts += f"--input-color-file {color_bfg}"

    command = (
        f"Bifrost query {extra_opts} -v -g {ref_file} -q {query_file} -o {outprefix}"
    )
    utils.syscall(command, quiet=not debug)
    return f"{outprefix}.tsv"


def query_tsv_to_bitarray_file(infile, outfile, combine_method="intersection"):
    if combine_method == "intersection":
        combine_function = all
    elif combine_method == "union":
        combine_function = any
    else:
        raise NotImplementedError(f"unknown method: {combine_method}")

    presence = bitarray.bitarray()
    with open(infile) as f:
        for d in csv.DictReader(f, delimiter="\t"):
            if int(d["query_name"]) != len(presence):
                raise Exception(
                    f"Expect query names to be numbered 0, 1, ..., but got this: {d['query_name']} in file {infile}"
                )
            if combine_function(v == "1" for k, v in d.items() if k != "query_name"):
                presence.append(1)
            else:
                presence.append(0)

    # Note: when writing to file, fills any empty bits in the last byte with zeros
    with open(outfile, "wb") as f:
        presence.tofile(f)

    return len(presence)


def load_bitarray_file(filename, expect_length):
    a = bitarray.bitarray()
    with open(filename, "rb") as f:
        a.fromfile(f)
    assert len(a) >= expect_length > 0
    while len(a) > expect_length:
        a.pop()
    return a


def build_and_query(
    ref_files,
    query_file,
    outprefix,
    k=None,
    kmer_ratio=None,
    debug=False,
    combine_method="intersection",
):
    ref_fofn = f"{outprefix}.fofn"
    with open(ref_fofn, "w") as f:
        print(*ref_files, sep="\n", file=f)
    build_prefix = f"{outprefix}.build"
    graph_file = build(ref_fofn, build_prefix, "color", k=k, debug=debug)
    color_file = graph_file.replace("gfa.gz", "color.bfg")
    got_tsv = query(
        graph_file,
        query_file,
        outprefix,
        color_bfg=color_file,
        kmer_ratio=kmer_ratio,
        debug=debug,
    )
    outfile = f"{outprefix}.bitarray"
    array_length = query_tsv_to_bitarray_file(
        got_tsv, outfile, combine_method=combine_method
    )
    debug = True
    if not debug:
        os.unlink(ref_fofn)
        os.unlink(graph_file)
        os.unlink(color_file)
        bfi = graph_file.rsplit(".", maxsplit=1)[0] + ".bfi"
        if os.path.exists(bfi):
            os.unlink(bfi)
        os.unlink(got_tsv)

    return array_length, outfile
