""" Module for extracting and concatenating assemblies. """

__author__ = "Berger, Phillip"

import os
import shutil
from pathlib import Path
from Bio import SeqIO
from xspect import file_io
from xspect.definitions import get_xspect_tmp_path, fasta_endings


def change_header(assemblies_path, species_accessions: dict):
    """Change the header of the assemblies to the species name."""
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
                    with open(file_path, "r", encoding="utf-8") as f:
                        sequence = ""
                        for line in f.readlines():
                            if line[0] != ">":
                                sequence += line
                        new_header = f">{name}\n"
                    with open(file_path, "w", encoding="utf-8") as f:
                        f.write(new_header)
                        f.write(sequence)


def copy_assemblies(unzipped_path, assemblies_path):
    """Copy all assemblies to a new directory."""
    os.mkdir(assemblies_path)
    for folder in os.listdir(unzipped_path):
        for root, _, files in os.walk(unzipped_path / str(folder)):
            for file in files:
                file_ending = file.split(".")[-1]
                if file_ending in fasta_endings:
                    file_path = Path(root) / file
                    shutil.copy(file_path, (assemblies_path / file))


def concatenate_bf(unzipped_path, concatenate_path):
    """Concatenate all assemblies for Bloom filter training."""

    all_assemblies = []

    # Make new directory
    os.mkdir(concatenate_path)

    # Open the fasta files for each species.
    for folder in os.listdir(unzipped_path):
        species_files = []
        # Walk through dirs to get all fasta files.
        for root, _, files in os.walk(unzipped_path / folder):
            for file in files:
                file_ending = file.split(".")[-1]
                if file_ending in fasta_endings:
                    species_files.append(Path(root) / file)
                    all_assemblies.append(".".join(str(file).split(".")[:-1]))

        # Gather all sequences and headers.
        sequences = []
        headers = []
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
        with open(species_path, "w", encoding="utf-8") as species_file:
            species_file.write(species_header)
            species_file.write(species_sequence)

    return all_assemblies


def save_all_assemblies(dir_path: Path, all_assemblies: list[str]):
    """Save all assemblies to a file."""
    path = dir_path / "all_bf_assemblies.txt"
    with open(path, "w", encoding="utf-8") as file:
        for assembly in all_assemblies:
            file.write(f"{assembly}\n")


def bf(dir_name: str, delete: bool):
    """Extract and concatenate assemblies for Bloom filter training."""
    dir_path = get_xspect_tmp_path() / dir_name
    zip_path = dir_path / "zip_files"
    unzipped_path = dir_path / "zip_files_extracted"
    concatenate_path = dir_path / "concatenate"
    file_io.extract_zip(zip_path, unzipped_path)
    all_assemblies = concatenate_bf(unzipped_path, concatenate_path)
    save_all_assemblies(dir_path, all_assemblies)
    if delete:
        file_io.delete_zip_files(zip_path)
        shutil.rmtree(zip_path, ignore_errors=False)
        shutil.rmtree(unzipped_path, ignore_errors=False)


def svm(species_accessions: dict, dir_name: str, delete: bool):
    """Extract and concatenate assemblies for generating SVM training data."""
    dir_path = get_xspect_tmp_path() / dir_name
    zip_path = dir_path / "training_data_zipped"
    unzipped_path = dir_path / "training_data_unzipped"
    assemblies_path = dir_path / "training_data"
    file_io.extract_zip(zip_path, unzipped_path)
    copy_assemblies(unzipped_path, assemblies_path)
    change_header(assemblies_path, species_accessions)
    if delete:
        file_io.delete_zip_files(zip_path)
        shutil.rmtree(zip_path, ignore_errors=False)
        shutil.rmtree(unzipped_path, ignore_errors=False)
