"""This module contains functions to select and download assemblies for SVM creation."""

from time import sleep
from loguru import logger
from xspect.train_filter.ncbi_api import download_assemblies


def select_assemblies(accessions):
    """Selects up to 4 assemblies, ideally assemblies that were not used for training the filters.

    :param accessions: All selected assembly accessions for every species.
    :type accessions: dict
    :return: Dict with species name as key and selected accessions as value.
    """

    all_accessions = {
        sci_name: curr_accessions[-4:]
        for sci_name, curr_accessions in accessions.items()
    }

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
