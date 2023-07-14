import logging
import os

from sprog import genome_download, utils


class Sample:
    def __init__(
        self,
        directory,
        name=None,
        data_source=None,
        accession=None,
        species=None,
        original_file=None,
        debug=False,
    ):
        self.root_dir = os.path.abspath(directory)
        if not os.path.exists(self.root_dir):
            os.mkdir(self.root_dir)
        self.metadata_json = os.path.join(self.root_dir, "metadata.json")
        self.genome_fasta = os.path.join(self.root_dir, "genome.fasta")
        self.debug = debug

        if not os.path.exists(self.metadata_json):
            assert name is not None
            assert species is not None
            if data_source == "download":
                assert accession is not None
                assert original_file is None
            elif data_source in ["symlink", "copy_file"]:
                assert accession is None
                assert original_file is not None
                original_file = os.path.abspath(original_file)
            else:
                raise NotImplementedError(f"data source '{data_source}' not recognised")

            self.metadata = {
                "name": name,
                "data_source": data_source,
                "accession": accession,
                "species": species,
                "original_file": original_file,
                "genome_obtained": False,
            }
            if self.metadata["data_source"] == "symlink":
                self.obtain_genome()
            else:
                self.write_metadata_json()
        else:
            self.metadata = utils.load_json(self.metadata_json)

    def write_metadata_json(self):
        utils.write_json(self.metadata, self.metadata_json)

    def obtain_genome(self):
        if self.metadata["genome_obtained"]:
            return

        if self.metadata["data_source"] == "download":
            logging.info(
                f"Downloading sample {self.metadata['name']}, {self.metadata['accession']}"
            )
            genome_download.download_genome(
                self.metadata["accession"], self.genome_fasta
            )
        else:
            assert self.metadata["data_source"] in ["copy_file", "symlink"]
            symlink = self.metadata["data_source"] == "symlink"
            method = (
                "Copying"
                if self.metadata["data_source"] == "copy_file"
                else "Symlinking"
            )
            logging.info(
                f"{method} sample {self.metadata['name']}, {self.metadata['original_file']}"
            )
            utils.copy_or_symlink(
                self.metadata["original_file"], self.genome_fasta, symlink=symlink
            )

        self.metadata["genome_obtained"] = True
        self.write_metadata_json()

    def has_genome(self):
        return self.metadata["genome_obtained"]

    def simulate_perfect_reads(self, outprefix, length=150, frag_length=400, depth=10):
        assert self.has_genome()
        command = f"fastaq to_perfect_reads {self.genome_fasta} - {frag_length} 1 {depth} {length} | fastaq deinterleave - {outprefix}_1.fq.gz {outprefix}_2.fq.gz"
        utils.syscall(command)
