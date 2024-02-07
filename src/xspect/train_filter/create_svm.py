import csv
import os
import pickle
from pathlib import Path
from time import sleep

from Bio import SeqIO
from loguru import logger

import xspect.BF_v2 as BF_v2
from xspect.file_io import (
    delete_non_fasta,
    get_accessions,
    get_file_paths,
    get_species_names,
)
from xspect.train_filter.ncbi_api import download_assemblies


def select_assemblies(accessions):
    """Selects up to 4 assemblies, ideally assemblies that were not used for training the filters.

    :param accessions: All selected assembly accessions for every species.
    :type accessions: dict
    :return: Dict with species name as key and selected accessions as value.
    """
    all_accessions = {}

    for sci_name, current_accessions in accessions.items():
        selected_accessions = []
        # Select 4 assemblies beginning from the last one.
        for i in range(len(current_accessions) - 1, -1, -1):
            selected_accessions.append(current_accessions[i])
            if len(selected_accessions) == 4:
                break

        all_accessions[sci_name] = selected_accessions

    return all_accessions


def get_svm_assemblies(all_accessions, dir_name):
    """Download assemblies for svm creation.

    :param all_accessions: Contains lists with all previously selected assemblies for every species.
    :type all_accessions: dict
    :param dir_name: Name of the parent directory.
    :type dir_name: str
    """
    # Select accessions for download.
    selected_accessions = select_assemblies(all_accessions)

    # Download assemblies.
    for sci_name, accessions in selected_accessions.items():
        sleep(5)
        logger.info("Downloading {name}", name=sci_name)
        file_name = sci_name + ".zip"
        download_assemblies.download_assemblies(
            accessions=accessions,
            dir_name=dir_name,
            target_folder="training_data_zipped",
            zip_file_name=file_name,
        )
    logger.info("Downloads finished")


def init_bf(genus, array_size, hashes=7, k=21):
    """Initializes bloomfilter.

    :param genus: Name of the genus.
    :type genus: str
    :param array_size: Size of the bloomfilter.
    :type array_size: int
    :param hashes: The number of hash functions the bf uses.
    :type hashes: int
    :param k: Length of k-mers.
    :type k: int
    :return: The bloomfilter object.
    """
    path = Path(os.getcwd()) / "filter"

    # Initialize bloomfilter for genus.
    BF = BF_v2.AbaumanniiBloomfilter(array_size)
    BF.set_arraysize(array_size)
    BF.set_hashes(hashes)
    BF.set_k(k)

    # Get all species names.
    names_path = path / "species_names" / ("Filter" + genus + ".txt")
    with open(names_path, "rb") as f:
        clonetypes = pickle.load(f)

    # Get bloomfilter paths.
    bf_path = path / genus
    paths = sorted(os.listdir(bf_path))
    for i in range(len(paths)):
        paths[i] = str(bf_path / str(paths[i]))
    # Setup bloomfilters.
    BF.read_clonetypes(paths, clonetypes)

    return BF


def perform_lookup(bloomfilter, files, file_paths, accessions, names, spacing):
    """Performs a lookup on a bloomfilter object and gives the scores as a list.

    :param bloomfilter: The bloomfilter object on which the lookup is performed.
    :param files: List of file names.
    :type files: list[str]
    :param file_paths: List with the file paths.
    :type file_paths: list[str]
    :param accessions: List of all accessions.
    :type accessions: list[str]
    :param names: List with all species names.
    :type names: list[str]
    :return: List with all scores of the lookup.
    """
    scores = list()
    BF = bloomfilter

    # Lookup.
    for i in range(len(files)):
        BF.number_of_kmeres = 0
        BF.hits_per_filter = [0] * BF.clonetypes

        for sequence in SeqIO.parse(file_paths[i], "fasta"):
            # Dominik: changed sample size to var
            for j in range(0, len(sequence.seq) - BF.k, spacing):
                BF.number_of_kmeres += 1
                BF.lookup(str(sequence.seq[j : j + BF.k]))

        score = BF.get_score()
        score = [str(x) for x in score]
        score = ",".join(score)
        scores.append(accessions[i] + "," + score + "," + names[i])

    return scores


# https://stackoverflow.com/questions/21431052/sort-list-of-strings-by-a-part-of-the-string
def sort_list(scores, names):
    """Sorts the scores list by species name.

    :param scores: The scores gathered by a lookup of a bloomfilter.
    :type scores: list
    :param names: List with all species names.
    :type names: list[str]
    :return: The sorted scores list.
    """
    scores.sort(key=lambda x: x.split(",")[-1][:2])
    names = [x for x in names if x != "none"]
    names = list(dict.fromkeys(names))
    scores.insert(0, sorted(names))
    scores[0] = ["File"] + scores[0] + ["Label"]

    for i in range(1, len(scores)):
        line = scores[i].split(",")
        scores[i] = line

    return scores


def save_csv(genus, scores):
    """Saves the scores as csv file.

    :param genus: Name of the genus.
    :type genus: str
    :param scores: The scores gathered by a lookup of a bloomfilter.
    :type scores: list
    """
    training_data_path = Path(os.getcwd()) / "Training_data"
    if not os.path.exists(training_data_path):
        os.mkdir(training_data_path)

    path = training_data_path / (genus + "_Training_data_spec.csv")
    with open(path, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerows(scores)


# Dominik: added spacing
def new_helper(spacing, genus, dir_name, array_size, k=21):
    """Create support vector machine for bloomfilters of a genus.

    :param spacing:
    :param genus: Name of the genus.
    :type genus: str
    :param dir_name: Name of the parent directory.
    :type dir_name: str
    :param array_size: Size for the byte array which is the bloomfilter.
    :type array_size: int
    :param k: Length of the k-mers.
    :type k: int
    """
    # Get all files.
    base_path = Path(os.getcwd()) / "genus_metadata" / dir_name / "training_data"
    files = os.listdir(base_path)

    # Delete all non fasta files.
    files = delete_non_fasta(files)

    # Get accessions from file names.
    accessions = get_accessions(files)

    # Get all complete file paths.
    file_paths = get_file_paths(base_path, files)

    # Get all species names from the header in the fasta files.
    names = get_species_names(file_paths)

    # Initialize bloomfilter.
    bf = init_bf(genus, array_size)

    # Perform lookup on bloomfilter.
    # Dominik: added spacing
    scores = perform_lookup(bf, files, file_paths, accessions, names, spacing)

    # Sort score list by species names.
    scores = sort_list(scores, names)

    # Save results in csv file.
    save_csv(genus, scores)
