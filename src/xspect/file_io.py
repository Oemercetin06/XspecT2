"""
File IO module.
"""

from linecache import getline
import os
from pathlib import Path
import zipfile

from loguru import logger


def check_folder_structure():
    """Checks the folder structure and creates new folders if needed."""
    # Create list of all folder paths.
    root_path = Path(os.getcwd())
    filter_path = root_path / "filter"
    meta_path = root_path / "genus_metadata"
    filter_folder_names = [
        "array_sizes",
        "Metagenomes",
        "species_names",
        "translation_dicts",
    ]
    folder_paths = [filter_path, meta_path]
    for filter_folder_name in filter_folder_names:
        filter_folder_path = filter_path / filter_folder_name
        folder_paths.append(filter_folder_path)

    # Check if folders exist. If not create them.
    for folder_path in folder_paths:
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)


def delete_non_fasta(files):
    """Delete all non fasta files from the list and return the list without those file names.

    :param files: List of file names.
    :type files: list[str]
    :return: List with only fasta files.
    """
    # All possible fasta file endings.
    fasta_endings = ["fasta", "fna", "fa", "ffn", "frn"]

    # Iterate through file list backwards and delete all non fasta files.
    for i in range(len(files) - 1, -1, -1):
        file = files[i].split(".")
        if file[-1] in fasta_endings:
            continue
        else:
            del files[i]

    return files


def get_accessions(file_names: list[str]) -> list[str]:
    """Extract accessions from file names.

    :param files: List of file names.
    :type files: list[str]
    :return: List of all accessions.
    :rtype: list[str]
    """
    accessions = []
    for idx, file in enumerate(file_names):
        accessions.append(file.split("_"))
        accessions[idx] = accessions[idx][0] + "_" + accessions[idx][1]

    return accessions


def get_file_paths(base_path: Path, file_names: list[str]) -> list[Path]:
    """Make a list with the paths to the files.

    :param base_path: Path of the parent directory.
    :type base_path: Path
    :param files: List of file names.
    :type files: list[str]
    :return: A list with all file paths.
    :rtype: list[Path]
    """
    return [base_path / file for file in file_names]


def get_species_names(file_paths: list[Path]):
    """Extracts the species names.

    :param file_paths: List with the file paths.
    :type file_paths: list[Path]
    :return: List with all species names.
    """
    names = list()
    for path in file_paths:
        header = getline(str(path), 1)
        name = header.replace("\n", "").replace(">", "")
        if not name.isdigit():
            logger.error(
                "The header of file: {path} does not contain a correct ID: {name}. The ID needs to be "
                "just numbers"
            )
            logger.error("Aborting")
            exit()
        names.append(name)
    return names


def delete_zip_files(dir_path):
    """Delete all zip files in the given directory."""
    files = os.listdir(dir_path)
    for file in files:
        if zipfile.is_zipfile(file):
            file_path = dir_path / str(file)
            os.remove(file_path)


def extract_zip(zip_path, unzipped_path):
    """Extracts all files from a directory with zip files."""
    # Make new directory.
    os.mkdir(unzipped_path)

    file_names = os.listdir(zip_path)
    for file in file_names:
        file_path = zip_path / file
        if zipfile.is_zipfile(file_path):
            with zipfile.ZipFile(file_path) as item:
                directory = unzipped_path / file.replace(".zip", "")
                item.extractall(directory)


def concatenate_meta(path: Path, genus: str):
    """Concatenates all concatenated fasta files that are used to train bloomfilters to one fasta file.

    :param path: Path to the directory with the concatenated fasta files.
    :type path: Path
    :param genus: Genus name.
    :type genus: str
    """
    files_path = path / "concatenate"
    fasta_endings = ["fasta", "fna", "fa", "ffn", "frn"]
    meta_path = path / (genus + ".fasta")
    files = os.listdir(files_path)

    with open(meta_path, "w") as meta_file:
        # Write the header.
        meta_header = f">{genus} metagenome\n"
        meta_file.write(meta_header)

        # Open each concatenated species file and write the sequence in the meta file.
        for file in files:
            file_ending = str(file).split(".")[-1]
            if file_ending in fasta_endings:
                with open((files_path / str(file)), "r") as species_file:
                    for line in species_file:
                        if line[0] != ">":
                            meta_file.write(line.replace("\n", ""))
