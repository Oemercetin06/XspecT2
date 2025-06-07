"""
File IO module.
"""

from json import loads
import os
from pathlib import Path
import zipfile
from typing import Callable, Iterator
from Bio import SeqIO
from xspect.definitions import fasta_endings, fastq_endings


def delete_zip_files(dir_path) -> None:
    """
    Delete all zip files in the given directory.

    This function checks each file in the specified directory and removes it if it is a zip file.

    Args:
        dir_path (Path): Path to the directory where zip files should be deleted.
    """
    files = os.listdir(dir_path)
    for file in files:
        if zipfile.is_zipfile(file):
            file_path = dir_path / str(file)
            os.remove(file_path)


def extract_zip(zip_path: Path, unzipped_path: Path) -> None:
    """
    Extracts all files from a zip file.

    Extracts the contents of the specified zip file to the given directory.

    Args:
        zip_path (Path): Path to the zip file to be extracted.
        unzipped_path (Path): Path to the directory where the contents will be extracted.
    """
    unzipped_path.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(zip_path) as item:
        item.extractall(unzipped_path)


def get_record_iterator(file_path: Path) -> Iterator:
    """
    Returns a record iterator for a fasta or fastq file.

    This function checks the file extension to determine if the file is in fasta or fastq format
    and returns an iterator over the records in the file using Biopython's SeqIO module.

    Args:
        file_path (Path): Path to the fasta or fastq file.

    Returns:
        Iterator: An iterator over the records in the file.

    Raises:
        ValueError: If the file path is not a Path object, does not exist, is not a file,
                    or has an invalid file format.
    """
    if not isinstance(file_path, Path):
        raise ValueError("Path must be a Path object")

    if not file_path.exists():
        raise ValueError("File does not exist")

    if not file_path.is_file():
        raise ValueError("Path must be a file")

    if file_path.suffix[1:] in fasta_endings:
        return SeqIO.parse(file_path, "fasta")

    if file_path.suffix[1:] in fastq_endings:
        return SeqIO.parse(file_path, "fastq")

    raise ValueError("Invalid file format, must be a fasta or fastq file")


def concatenate_species_fasta_files(
    input_folders: list[Path], output_directory: Path
) -> None:
    """
    Concatenate fasta files from different species into one file per species.

    This function iterates through each species folder within the given input folder,
    collects all fasta files, and concatenates their contents into a single fasta file
    named after the species.

    Args:
        input_folders (list[Path]): List of paths to species folders.
        output_directory (Path): Path to the output directory.
    """
    for species_folder in input_folders:
        species_name = species_folder.name
        fasta_files = [
            f for ending in fasta_endings for f in species_folder.glob(f"*.{ending}")
        ]
        if len(fasta_files) == 0:
            raise ValueError(f"no fasta files found in {species_folder}")

        # concatenate fasta files
        concatenated_fasta = output_directory / f"{species_name}.fasta"
        with open(concatenated_fasta, "w", encoding="utf-8") as f:
            for fasta_file in fasta_files:
                with open(fasta_file, "r", encoding="utf-8") as f_in:
                    f.write(f_in.read())


def concatenate_metagenome(fasta_dir: Path, meta_path: Path) -> None:
    """
    Concatenate all fasta files in a directory into one file.

    This function searches for all fasta files in the specified directory and writes their contents
    into a single output file. The output file will contain the concatenated sequences from all fasta files.

    Args:
        fasta_dir (Path): Path to the directory with the fasta files.
        meta_path (Path): Path to the output file.
    """
    fasta_files = [
        file for ending in fasta_endings for file in fasta_dir.glob(f"*.{ending}")
    ]
    with open(meta_path, "w", encoding="utf-8") as meta_file:
        for fasta_file in fasta_files:
            with open(fasta_file, "r", encoding="utf-8") as f_in:
                meta_file.write(f_in.read())


def get_ncbi_dataset_accession_paths(
    ncbi_dataset_path: Path,
) -> dict[str, Path]:
    """
    Get the paths of the NCBI dataset accessions.

    This function reads the dataset catalog from the NCBI dataset directory and returns a dictionary
    mapping each accession to its corresponding file path. The first item in the dataset catalog is
    assumed to be a data report, and is skipped.

    Args:
        ncbi_dataset_path (Path): Path to the NCBI dataset directory.

    Returns:
        dict[str, Path]: Dictionary with the accession as key and the path as value.

    Raises:
        ValueError: If the dataset path does not exist or is invalid.
    """
    data_path = ncbi_dataset_path / "ncbi_dataset" / "data"
    if not data_path.exists():
        raise ValueError(f"Path {data_path} does not exist.")

    accession_paths = {}
    with open(data_path / "dataset_catalog.json", "r", encoding="utf-8") as f:
        res = loads(f.read())
        for assembly in res["assemblies"][1:]:  # the first item is the data report
            accession = assembly["accession"]
            assembly_path = data_path / assembly["files"][0]["filePath"]
            accession_paths[accession] = assembly_path
    return accession_paths


def filter_sequences(
    input_file: Path,
    output_file: Path,
    included_ids: list[str],
) -> None:
    """
    Filter sequences by IDs from an input file and save them to an output file.

    This function reads a fasta or fastq file, filters the sequences based on the provided IDs,
    and writes the matching sequences to an output file. If no IDs are provided, no output file
    is created.

    Args:
        input_file (Path): Path to the input file.
        output_file (Path): Path to the output file.
        included_ids (list[str], optional): List of IDs to include. If None, no output file
            is created.
    """
    if not included_ids:
        print("No IDs provided, no output file will be created.")
        return

    with open(output_file, "w", encoding="utf-8") as out_f:
        for record in get_record_iterator(input_file):
            if record.id in included_ids:
                SeqIO.write(record, out_f, "fasta")


def prepare_input_output_paths(
    input_path: Path,
) -> tuple[list[Path], Callable[[int, Path], Path]]:
    """
    Processes the input path into a list of input paths and a function generating output paths.

    This function checks if the input path is a directory or a file. If it is a directory,
    it collects all files with specified fasta and fastq endings. If it is a file, it uses that file
    as the input path. It then returns a list of input file paths and a function that generates
    output paths based on the index of the input file and a specified output path.

    Args:
        input_path (Path): Path to the directory or file.

    Returns:
        tuple[list[Path], Callable[[int, Path], Path]]: A tuple containing:
            - A list of input file paths
            - A function that takes an index and the output path,
              and returns the processed output path.

    Raises:
        ValueError: If the input path is invalid.
    """
    input_is_dir = input_path.is_dir()
    ending_wildcards = [f"*.{ending}" for ending in fasta_endings + fastq_endings]

    if input_is_dir:
        input_paths = [p for e in ending_wildcards for p in input_path.glob(e)]
    elif input_path.is_file():
        input_paths = [input_path]
    else:
        raise ValueError("Invalid input path")

    def get_output_path(idx: int, output_path: Path) -> Path:
        return (
            output_path.parent / f"{output_path.stem}_{idx+1}{output_path.suffix}"
            if input_is_dir
            else output_path
        )

    return input_paths, get_output_path
