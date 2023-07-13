import os
import pytest

import pyfastaq

from sprog import genome_download, utils


def test_download_genome_from_genbank():
    # We just want a small genome to test, saving internet bandwidth. This is
    # an HIV-1 genome
    genbank_id = "NC_001802.1"
    fasta_out = "tmp.download_genome_from_genbank.fa"
    utils.syscall(f"rm -f {fasta_out}")
    genome_download.download_genome(genbank_id, fasta_out)
    assert os.path.exists(fasta_out)
    got_seqs = {}
    pyfastaq.tasks.file_to_dict(fasta_out, got_seqs)
    assert len(got_seqs) == 1
    assert genbank_id in got_seqs
    assert len(got_seqs[genbank_id]) == 9181
    os.unlink(fasta_out)


def test_download_ncbi_assembly():
    fasta_out = "tmp.download_ncbi_assembly.fa"
    utils.syscall(f"rm -rf {fasta_out}*")
    genome_download.download_genome("GCA_000069185.1", fasta_out)
    assert os.path.exists(fasta_out)
    got_seqs = {}
    pyfastaq.tasks.file_to_dict(fasta_out, got_seqs)
    assert len(got_seqs) == 2
    os.unlink(fasta_out)
