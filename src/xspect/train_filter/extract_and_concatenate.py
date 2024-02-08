"""

"""

__author__ = "Berger, Phillip"

import os
import shutil
from pathlib import Path
from Bio import SeqIO
from loguru import logger
import xspect.file_io as file_io


class ExtractConcatenate:
    _fasta_endings = ["fasta", "fna", "fa", "ffn", "frn"]

    def __init__(self, dir_name: str, delete: bool):
        self._path = Path(os.getcwd()) / "genus_metadata" / dir_name
        self._dir_name = dir_name
        self._delete = delete
        self._all_assemblies: list[str] = list()

    def bf(self):
        zip_path = self._path / "zip_files"
        unzipped_path = self._path / "zip_files_extracted"
        concatenate_path = self._path / "concatenate"
        logger.info("Extracting assemblies")
        file_io.extract_zip(zip_path, unzipped_path)
        logger.info("Concatenating assemblies")
        self._concatenate_bf(unzipped_path, concatenate_path)
        self._save_all_assemblies()
        if self._delete:
            logger.info("Deleting copies")
            file_io.delete_zip_files(zip_path)
            shutil.rmtree(zip_path, ignore_errors=False)
            shutil.rmtree(unzipped_path, ignore_errors=False)

    def svm(self, species_accessions: dict):
        zip_path = self._path / "training_data_zipped"
        unzipped_path = self._path / "training_data_unzipped"
        assemblies_path = self._path / "training_data"
        logger.info("Extracting assemblies")
        file_io.extract_zip(zip_path, unzipped_path)
        logger.info("Copying assemblies into one folder")
        self._copy_assemblies(unzipped_path, assemblies_path)
        logger.info("Changing headers")
        self._change_header(assemblies_path, species_accessions)
        if self._delete:
            logger.info("Deleting copies")
            file_io.delete_zip_files(zip_path)
            shutil.rmtree(zip_path, ignore_errors=False)
            shutil.rmtree(unzipped_path, ignore_errors=False)

    def _copy_assemblies(self, unzipped_path, assemblies_path):
        os.mkdir(assemblies_path)
        for folder in os.listdir(unzipped_path):
            for root, dirs, files in os.walk(unzipped_path / str(folder)):
                for file in files:
                    file_ending = file.split(".")[-1]
                    if file_ending in self._fasta_endings:
                        file_path = Path(root) / file
                        shutil.copy(file_path, (assemblies_path / file))

    @staticmethod
    def _change_header(assemblies_path, species_accessions: dict):
        files = os.listdir(assemblies_path)
        # Iterate through all species.
        for name, accessions in species_accessions.items():
            # Iterate through all accessions of the current species.
            for accession in accessions:
                # Iterate through all file names.
                for file in files:
                    if accession in file:
                        file_path = assemblies_path / str(file)
                        # Change the header.
                        with open(file_path, "r") as f:
                            sequence = ""
                            for line in f.readlines():
                                if line[0] != ">":
                                    sequence += line
                            new_header = f">{name}\n"
                        with open(file_path, "w") as f:
                            f.write(new_header)
                            f.write(sequence)

    def _concatenate_bf(self, unzipped_path, concatenate_path):
        # Make new directory
        os.mkdir(concatenate_path)

        genus = self._dir_name.split("_")[0]
        meta_path = self._path / (genus + ".fasta")

        # Open the fasta files for each species.
        for folder in os.listdir(unzipped_path):
            species_files = list()
            # Walk through dirs to get all fasta files.
            for root, dirs, files in os.walk(unzipped_path / folder):
                for file in files:
                    file_ending = file.split(".")[-1]
                    if file_ending in self._fasta_endings:
                        species_files.append(Path(root) / file)
                        self._all_assemblies.append(".".join(str(file).split(".")[:-1]))

            # Gather all sequences and headers.
            sequences = list()
            headers = list()
            for file in species_files:
                records = SeqIO.parse(file, "fasta")
                for record in records:
                    headers.append(record.id)
                    sequences.append(str(record.seq))

            # Concatenate sequences
            species_sequence = "".join(sequences)
            species_header = ">" + " § ".join(headers) + "\n"

            # Save concatenated sequences and headers
            species_path = concatenate_path / (folder + ".fasta")
            with open(species_path, "w") as species_file:
                species_file.write(species_header)
                species_file.write(species_sequence)

    def _save_all_assemblies(self):
        path = self._path / "all_bf_assemblies.txt"
        with open(path, "w") as file:
            for assembly in self._all_assemblies:
                file.write(f"{assembly}\n")
