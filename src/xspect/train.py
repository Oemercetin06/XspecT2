"""
This module contains the main functions for training the models.
"""

import shutil
from pathlib import Path
from tempfile import TemporaryDirectory
from loguru import logger
from xspect.definitions import get_xspect_model_path
from xspect.file_io import (
    concatenate_species_fasta_files,
    concatenate_metagenome,
    extract_zip,
    get_ncbi_dataset_accession_paths,
)
from xspect.models.probabilistic_filter_model import ProbabilisticFilterModel
from xspect.models.probabilistic_filter_svm_model import ProbabilisticFilterSVMModel
from xspect.models.probabilistic_single_filter_model import (
    ProbabilisticSingleFilterModel,
)
from xspect.ncbi import AssemblySource, NCBIHandler


def train_from_directory(
    display_name: str,
    dir_path: Path,
    meta: bool = False,
    training_accessions: dict[str, list[str]] | None = None,
    svm_accessions: dict[str, list[str]] | None = None,
    svm_step: int = 1,
    translation_dict: dict[str, str] | None = None,
    author: str | None = None,
    author_email: str | None = None,
):
    """
    Train a model from a directory containing training data.

    This function trains a probabilistic filter model using the data in the specified directory.
    The training data should be organized in the following way:
    - dir_path
        - cobs
            - <species_name_1>
                - <fasta_file_1>
                - <fasta_file_2>
            - <species_name_2>
                - <fasta_file_1>
                - <fasta_file_2>
        - svm (optional)
            - <species_name_1>
                - <svm_file_1>
                - <svm_file_2>
            - <species_name_2>
                - <svm_file_1>
                - <svm_file_2>
    If no SVM directory is found, the model will be trained without SVM.
    The training data should be in FASTA format. The model is saved to the xspect_data directory.

    Args:
        display_name (str): Name of the model to be trained.
        dir_path (Path): Path to the directory containing training data.
        meta (bool, optional): Whether to train a metagenome model. Defaults to False.
        training_accessions (list[str], optional): List of training accessions. Defaults to None.
        svm_accessions (list[str], optional): List of SVM accession identifiers. Defaults to None.
        svm_step (int, optional): Step size for SVM training. Defaults to 1.
        translation_dict (dict[str, str], optional): Dictionary for display names. Defaults to None.
        author (str, optional): Author of the model. Defaults to None.
        author_email (str, optional): Author's email. Defaults to None.

    Raises:
        TypeError: If `display_name` is not a string.
        TypeError: If `dir_path` is not a Path object to a valid directory.
        ValueError: If the "cobs" directory is not found in `dir_path`.
        ValueError: If no folders are found in the "cobs" directory.
        ValueError: If the number of SVM folders does not match the number of COBS folders.
        ValueError: If the names of COBS folders and SVM folders do not match.
        ValueError: If no FASTA files are found in a COBS folder.

    Notes:
        - If the "svm" directory is not found, the model will be trained without SVM.
        - Temporary directories are used for intermediate processing.
    """

    if not isinstance(display_name, str):
        raise TypeError("display_name must be a string")

    if not isinstance(dir_path, Path) and dir_path.exists() and dir_path.is_dir():
        raise TypeError("dir must be Path object to a valid directory")

    cobs_training_path = dir_path / "cobs"
    if not cobs_training_path.exists():
        raise ValueError("cobs directory not found")

    cobs_folders = [f for f in cobs_training_path.iterdir() if f.is_dir()]
    if len(cobs_folders) == 0:
        raise ValueError("no folders found in cobs directory")

    svm_path = dir_path / "svm"
    if svm_path.exists():
        svm_folders = [f for f in svm_path.iterdir() if f.is_dir()]
        if len(svm_folders) != len(cobs_folders):
            raise ValueError(
                "number of svm folders does not match number of cobs folders"
            )

        for cobs_folder, svm_folder in zip(cobs_folders, svm_folders):
            if cobs_folder.name != svm_folder.name:
                raise ValueError("cobs folder and svm folder names do not match")
    else:
        print("SVM directory not found. Model will be trained without SVM.")

    with TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)
        species_dir = tmp_dir / "species"
        species_dir.mkdir(parents=True, exist_ok=True)

        logger.info("Concatenating genomes for species training...")
        concatenate_species_fasta_files(cobs_folders, species_dir)

        if svm_path.exists():
            logger.info("Training species SVM model...")
            species_model = ProbabilisticFilterSVMModel(
                k=21,
                model_display_name=display_name,
                author=author,
                author_email=author_email,
                model_type="Species",
                base_path=get_xspect_model_path(),
                kernel="rbf",
                c=1.0,
            )
            species_model.fit(
                species_dir,
                svm_path,
                display_names=translation_dict,
                svm_step=svm_step,
                training_accessions=training_accessions,
                svm_accessions=svm_accessions,
            )
        else:
            logger.info("Training species model...")
            species_model = ProbabilisticFilterModel(
                k=21,
                model_display_name=display_name,
                author=author,
                author_email=author_email,
                model_type="Species",
                base_path=get_xspect_model_path(),
            )
            species_model.fit(
                species_dir,
                display_names=translation_dict,
                training_accessions=training_accessions,
            )

        species_model.save()

        if meta:
            logger.info("Concatenating genomes for metagenome training...")
            meta_fasta = tmp_dir / f"{display_name}.fasta"
            concatenate_metagenome(species_dir, meta_fasta)

            logger.info("Training metagenome model...")
            genus_model = ProbabilisticSingleFilterModel(
                k=21,
                model_display_name=display_name,
                author=author,
                author_email=author_email,
                model_type="Genus",
                base_path=get_xspect_model_path(),
            )
            genus_model.fit(
                meta_fasta,
                display_name,
                training_accessions=(
                    sum(training_accessions.values(), [])
                    if training_accessions
                    else None
                ),
            )
            genus_model.save()


def train_from_ncbi(
    genus: str,
    svm_step: int = 1,
    author: str | None = None,
    author_email: str | None = None,
    ncbi_api_key: str | None = None,
):
    """
    Train a model using NCBI assembly data for a given genus.

    This function trains a probabilistic filter model using the assembly data from NCBI.
    The training data is downloaded and processed, and the model is saved to the
    xspect_data directory.

    Args:
        genus (str): Genus name for which the model will be trained.
        svm_step (int, optional): Step size for SVM training. Defaults to 1.
        author (str, optional): Author of the model. Defaults to None.
        author_email (str, optional): Author's email. Defaults to None.
        ncbi_api_key (str, optional): NCBI API key for accessing NCBI resources. Defaults to None.

    Raises:
        TypeError: If `genus` is not a string.
        ValueError: If no species with accessions are found.

    Notes:
        - The function uses NCBI API to fetch assembly metadata.
        - Temporary directories are used for intermediate processing.
    """
    if not isinstance(genus, str):
        raise TypeError("genus must be a string")

    logger.info("Getting NCBI metadata...")
    ncbi_handler = NCBIHandler(api_key=ncbi_api_key)
    genus_tax_id = ncbi_handler.get_genus_taxon_id(genus)
    species_ids = ncbi_handler.get_species(genus_tax_id)
    species_names = ncbi_handler.get_taxon_names(species_ids)

    filtered_species_ids = [
        tax_id
        for tax_id in species_ids
        if "candidatus" not in species_names[tax_id].lower()
        and " sp." not in species_names[tax_id].lower()
    ]
    filtered_species_names = {
        str(tax_id): species_names[tax_id] for tax_id in filtered_species_ids
    }

    accessions = {}
    for tax_id in filtered_species_ids:
        taxon_accessions = ncbi_handler.get_highest_quality_accessions(
            tax_id, AssemblySource.REFSEQ, 8
        )
        if not taxon_accessions:
            logger.warning(f"No assemblies found for tax_id {tax_id}. Skipping.")
            filtered_species_names.pop(str(tax_id), None)
            continue
        accessions[tax_id] = taxon_accessions

    if not accessions:
        raise ValueError(
            "No species with accessions found. Please check the genus name."
        )

    with TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)
        cobs_dir = tmp_dir / "cobs"
        svm_dir = tmp_dir / "svm"
        cobs_dir.mkdir(parents=True, exist_ok=True)
        svm_dir.mkdir(parents=True, exist_ok=True)

        logger.info("Downloading genomes from NCBI...")
        all_accessions = sum(accessions.values(), [])
        batch_size = 100
        accession_paths = {}
        for i in range(0, len(all_accessions), batch_size):
            batch = all_accessions[i : i + batch_size]
            ncbi_handler.download_assemblies(accessions=batch, output_dir=tmp_dir)
            extract_zip(
                tmp_dir / "ncbi_dataset.zip", tmp_dir / f"batch-{i}-{i+batch_size}"
            )
            accession_paths.update(
                get_ncbi_dataset_accession_paths(tmp_dir / f"batch-{i}-{i+batch_size}")
            )

        # select accessions
        cobs_accessions = {}
        svm_accessions = {}
        for tax_id, accession_list in accessions.items():
            cobs_accessions[tax_id] = accession_list[:4]
            svm_accessions[tax_id] = accession_list[-4:]

        # move files
        for tax_id, accession_list in cobs_accessions.items():
            tax_id_dir = cobs_dir / str(tax_id)
            tax_id_dir.mkdir(parents=True, exist_ok=True)
            for accession in accession_list:
                accession_path = accession_paths[accession]
                shutil.copy(accession_path, tax_id_dir / f"{accession}.fasta")
        for tax_id, accession_list in svm_accessions.items():
            tax_id_dir = svm_dir / str(tax_id)
            tax_id_dir.mkdir(parents=True, exist_ok=True)
            for accession in accession_list:
                accession_path = accession_paths[accession]
                shutil.copy(accession_path, tax_id_dir / f"{accession}.fasta")

        train_from_directory(
            display_name=genus,
            dir_path=tmp_dir,
            meta=True,
            training_accessions=cobs_accessions,
            svm_accessions=svm_accessions,
            svm_step=svm_step,
            translation_dict=filtered_species_names,
            author=author,
            author_email=author_email,
        )
