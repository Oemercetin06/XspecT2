"""
File IO module.
"""

from json import loads
import os
from pathlib import Path
import zipfile
from Bio import SeqIO
from xspect.definitions import fasta_endings, fastq_endings


def delete_zip_files(dir_path):
    """Delete all zip files in the given directory."""
    files = os.listdir(dir_path)
    for file in files:
        if zipfile.is_zipfile(file):
            file_path = dir_path / str(file)
            os.remove(file_path)


def extract_zip(zip_path: Path, unzipped_path: Path):
    """Extracts all files from a directory with zip files."""
    # Make new directory.
    unzipped_path.mkdir(parents=True, exist_ok=True)

    file_names = os.listdir(zip_path)
    for file in file_names:
        file_path = zip_path / file
        if zipfile.is_zipfile(file_path):
            with zipfile.ZipFile(file_path) as item:
                directory = unzipped_path / file.replace(".zip", "")
                item.extractall(directory)


def concatenate_meta(path: Path, genus: str):
    """Concatenates all species files to one fasta file.

    :param path: Path to the directory with the concatenated fasta files.
    :type path: Path
    :param genus: Genus name.
    :type genus: str
    """
    files_path = path / "concatenate"
    meta_path = path / (genus + ".fasta")
    files = os.listdir(files_path)

    with open(meta_path, "w", encoding="utf-8") as meta_file:
        # Write the header.
        meta_header = f">{genus} metagenome\n"
        meta_file.write(meta_header)

        # Open each concatenated species file and write the sequence in the meta file.
        for file in files:
            file_ending = str(file).rsplit(".", maxsplit=1)[-1]
            if file_ending in fasta_endings:
                with open(
                    (files_path / str(file)), "r", encoding="utf-8"
                ) as species_file:
                    for line in species_file:
                        if line[0] != ">":
                            meta_file.write(line.replace("\n", ""))


def get_record_iterator(file_path: Path):
    """Returns a record iterator for a fasta or fastq file."""
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


def get_records_by_id(file: Path, ids: list[str]):
    """Return records with the specified ids."""
    records = get_record_iterator(file)
    return [record for record in records if record.id in ids]


def concatenate_species_fasta_files(input_folders: list[Path], output_directory: Path):
    """Concatenate fasta files from different species into one file per species.

    Args:
        input_species_folders (list[Path]): List of paths to species folders.
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


def concatenate_metagenome(fasta_dir: Path, meta_path: Path):
    """Concatenate all fasta files in a directory into one file.

    Args:
        fasta_dir (Path): Path to the directory with the fasta files.
        meta_path (Path): Path to the output file.
    """
    with open(meta_path, "w", encoding="utf-8") as meta_file:
        for fasta_file in fasta_dir.glob("*.fasta"):
            with open(fasta_file, "r", encoding="utf-8") as f_in:
                meta_file.write(f_in.read())


def get_ncbi_dataset_accession_paths(
    ncbi_dataset_path: Path,
) -> dict[str, Path]:
    """Get the paths of the NCBI dataset accessions.

    Args:
        ncbi_dataset_path (Path): Path to the NCBI dataset directory.

    Returns:
        dict[str, Path]: Dictionary with the accession as key and the path as value.
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
):
    """Filter sequences by IDs from an input file and save them to an output file.

    Args:
        input_file (Path): Path to the input file.
        output_file (Path): Path to the output file.
        included_ids (list[str], optional): List of IDs to include. If None, no output file is created.
    """
    if not included_ids:
        print("No IDs provided, no output file will be created.")
        return

    with open(output_file, "w", encoding="utf-8") as out_f:
        for record in get_record_iterator(input_file):
            if record.id in included_ids:
                SeqIO.write(record, out_f, "fasta")
