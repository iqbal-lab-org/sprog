import logging
import os
import zipfile

import pyfastaq

from sprog import utils


def download_genome_from_genbank(genbank_id, fasta_out):
    tmp_fasta = f"{fasta_out}.tmp"
    command = f"wget -O {tmp_fasta} 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id={genbank_id}'"
    # The file always has an empty line at the end. Remove it, and while we're
    # here remove anything after whitespace in header lines
    utils.syscall(command)
    reader = pyfastaq.sequences.file_reader(tmp_fasta)
    with open(fasta_out, "w") as f:
        for seq in reader:
            seq.id = seq.id.split()[0]
            print(seq, file=f)
    os.unlink(tmp_fasta)


def download_ncbi_assembly(accession, fasta_out):
    tmp_dir = f"{fasta_out}.tmp"
    utils.syscall(f"rm -rf {tmp_dir}")
    os.mkdir(tmp_dir)
    command = f"datasets download genome accession {accession} --no-progressbar"
    utils.syscall(command, cwd=tmp_dir)
    expect_zip = os.path.join(tmp_dir, "ncbi_dataset.zip")
    if not os.path.exists(expect_zip):
        raise Exception(f"Expected downloaded file not found: {expect_zip}")
    with zipfile.ZipFile(expect_zip, "r") as z:
        z.extractall(tmp_dir)
    fasta_dir = os.path.join(tmp_dir, "ncbi_dataset", "data", accession)
    if not os.path.exists(fasta_dir):
        raise Exception(
            f"Did not find expected directory of FASTA file for {accession}: {fasta_dir}"
        )
    fasta_files = list(os.listdir(fasta_dir))
    if len(fasta_files) != 1:
        raise Exception(
            f"Did not find exactly 1 FASTA file for {accession}: "
            + ";".join(fasta_files)
        )
    logging.info(f"Got FASTA for {accession}: {fasta_files[0]}")
    os.rename(os.path.join(fasta_dir, fasta_files[0]), fasta_out)
    logging.info(f"Cleaning up temp files made during downloading {accession}")
    utils.syscall(f"rm -r {tmp_dir}")


def download_genome(accession, fasta_out):
    if "_" not in accession:
        raise Exception(f"Accession format not recognised: {accession}")

    prefix = accession.split("_")[0]

    if len(prefix) == 2:
        download_genome_from_genbank(accession, fasta_out)
    elif len(prefix) == 3:
        download_ncbi_assembly(accession, fasta_out)
    else:
        raise Exception(f"Accession format not recognised: {accession}")
