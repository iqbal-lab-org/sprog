import os
import pytest

from sprog import genome_download, sample, utils


@pytest.fixture(scope="module")
def hiv_genome_fasta():
    accession = "NC_001802.1"
    hiv_genome_fasta = os.path.abspath("tmp.NC_001802.1.fa")
    genome_download.download_genome(accession, hiv_genome_fasta)
    yield hiv_genome_fasta
    os.unlink(hiv_genome_fasta)


def test_sample_from_genbank():
    # This is small HIV-1 genome so won't use up bandwidth downloading.
    sample_dir = "tmp.sample.from_genbank"
    utils.syscall(f"rm -rf {sample_dir}")
    genbank_id = "NC_001802.1"
    smpl = sample.Sample(
        sample_dir,
        name="test_sample",
        data_source="download",
        accession=genbank_id,
        species="hiv",
        debug=True,
    )

    expect_metadata = {
        "name": "test_sample",
        "data_source": "download",
        "accession": "NC_001802.1",
        "species": "hiv",
        "original_file": None,
        "genome_obtained": False,
    }

    assert os.path.exists(sample_dir)
    assert smpl.metadata == expect_metadata
    assert not os.path.exists(smpl.genome_fasta)
    smpl.obtain_genome()
    expect_metadata["genome_obtained"] = True
    assert smpl.metadata == expect_metadata
    assert os.path.exists(smpl.genome_fasta)
    assert smpl.has_genome()

    # create the sample object again to check metadata loaded correctly
    smpl = sample.Sample(sample_dir)
    assert smpl.metadata == expect_metadata
    assert smpl.has_genome()
    utils.syscall(f"rm -r {sample_dir}")


def test_sample_copy_genome(hiv_genome_fasta):
    sample_dir = "tmp.sample.copy_genome"
    utils.syscall(f"rm -rf {sample_dir}")
    smpl = sample.Sample(
        sample_dir,
        name="test_sample",
        data_source="copy_file",
        accession=None,
        original_file=hiv_genome_fasta,
        species="hiv",
        debug=True,
    )
    expect_metadata = {
        "name": "test_sample",
        "data_source": "copy_file",
        "accession": None,
        "species": "hiv",
        "original_file": hiv_genome_fasta,
        "genome_obtained": False,
    }
    assert smpl.metadata == expect_metadata
    assert not smpl.has_genome()
    assert not os.path.exists(smpl.genome_fasta)
    smpl.obtain_genome()
    assert os.path.exists(smpl.genome_fasta)
    utils.syscall(f"rm -r {sample_dir}")


def test_sample_symlink_genome(hiv_genome_fasta):
    sample_dir = "tmp.sample.symlink_genome"
    utils.syscall(f"rm -rf {sample_dir}")
    smpl = sample.Sample(
        sample_dir,
        name="test_sample",
        data_source="symlink",
        accession=None,
        original_file=hiv_genome_fasta,
        species="hiv",
        debug=True,
    )
    expect_metadata = {
        "name": "test_sample",
        "data_source": "symlink",
        "accession": None,
        "species": "hiv",
        "original_file": hiv_genome_fasta,
        "genome_obtained": True,
    }
    assert smpl.metadata == expect_metadata
    assert smpl.has_genome()
    assert os.path.exists(smpl.genome_fasta)
    smpl.obtain_genome()
    utils.syscall(f"rm -r {sample_dir}")
