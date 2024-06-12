import os
import shutil
from pathlib import Path
from platform import python_version
import sys
from time import localtime, perf_counter, asctime, sleep
from loguru import logger
from xspect import file_io
from xspect.definitions import get_xspect_tmp_path
from xspect.file_io import concatenate_meta
from xspect.models.probabilistic_filter_svm_model import ProbabilisticFilterSVMModel
from xspect.models.probabilistic_single_filter_model import (
    ProbabilisticSingleFilterModel,
)
from xspect.train_filter.ncbi_api import (
    ncbi_assembly_metadata,
    ncbi_taxon_metadata,
    ncbi_children_tree,
    download_assemblies,
)
from xspect.train_filter import (
    create_svm,
    html_scrap,
    extract_and_concatenate,
    get_paths,
    interface_XspecT,
)


def check_user_input(user_input: str):
    """The given input of the user will be checked. The input has to be a genus in NCBI.

    :return: The genus name.
    """
    taxon_metadata = ncbi_taxon_metadata.NCBITaxonMetadata([user_input])
    all_metadata = taxon_metadata.get_metadata()
    for metadata in all_metadata.values():
        sci_name = metadata["sci_name"]
        tax_id = metadata["tax_id"]
        rank = metadata["rank"]
        lineage = metadata["lineage"]
        bacteria_id = 2
        if not sci_name == user_input and not tax_id == user_input:
            print(
                f"{get_current_time()}| The given genus: {user_input} was found as genus: {sci_name} "
                f"ID: {tax_id}"
            )
            print(f"{get_current_time()}| Using {sci_name} as genus name.")
        if rank == "GENUS":
            if bacteria_id not in lineage:
                print(f"{get_current_time()}| The given genus is not a bacteria.")
                print(f"{get_current_time()}| Do you want to continue: [y/n]")
                choice = input("-> ").lower()
                if choice == "y":
                    return str(sci_name)
                else:
                    print(f"{get_current_time()}| Exiting...")
                    exit()
            else:
                return str(sci_name)
        else:
            print(f"{get_current_time()}| {user_input} is rank {rank} and not genus.")
            exit()


def copy_custom_data(bf_path: str, svm_path: str, dir_name: str):
    """

    :param bf_path:
    :param svm_path:
    :param dir_name:
    :return:
    """
    path = Path(os.getcwd()) / "genus_metadata" / dir_name
    new_bf_path = path / "concatenate"
    new_svm_path = path / "training_data"

    # Make the new directories.
    path.mkdir(exist_ok=True)
    new_bf_path.mkdir(exist_ok=True)
    new_svm_path.mkdir(exist_ok=True)

    # Move bloomfilter files.
    bf_files = os.listdir(bf_path)
    for file in bf_files:
        file_path = Path(bf_path) / file
        new_file_path = new_bf_path / file
        shutil.copy2(file_path, new_file_path)

    # Move svm files.
    svm_files = os.listdir(svm_path)
    for file in svm_files:
        file_path = Path(svm_path) / file
        new_file_path = new_svm_path / file
        shutil.copy2(file_path, new_file_path)


def set_logger(dir_name: str):
    """Sets the logger parameters.

    :param dir_name: Name of the folder where the log should be saved.
    """
    genus = dir_name.split("_")[0]

    # Starting logger.
    logger.remove()
    logger.add(sys.stderr, format="{time:HH:mm:ss} | {level} | {message}", level="INFO")
    log_path = get_xspect_tmp_path() / dir_name / (genus + ".log")
    logger.add(log_path, format="{time:HH:mm:ss} | {level} | {message}", level="DEBUG")


def create_translation_dict(dir_name: str) -> dict[str, str]:
    """Create a translation dictionary to translate the taxon ID to its scientific name from the file names.

    :param dir_name: Directory name for current genus.
    :return: The created translation dictionary.
    """
    path = get_xspect_tmp_path() / dir_name / "concatenate"
    files = os.listdir(path)
    translation_dict = dict()
    for file in files:
        file_split = file.split(".")[0].split("_")
        tax_id = file_split[0]
        final_file_name = tax_id + ".fasta"
        name = file_split[1]
        translation_dict[final_file_name] = name

    return translation_dict


def change_bf_assembly_file_names(dir_name: str):
    """Change all concatenated assembly names to only the taxon ID.

    :param dir_name: Directory name for current genus.
    """
    path = get_xspect_tmp_path() / dir_name / "concatenate"
    files = os.listdir(path)
    for file in files:
        file_split = file.split(".")[0].split("_")
        tax_id = file_split[0]
        new_file_name = f"{tax_id}.fasta"
        os.rename((path / file), (path / new_file_name))


def get_current_time():
    """Returns the current time in the form hh:mm:ss."""
    return asctime(localtime()).split()[3]


def train_ncbi(genus: str, svm_step: int = 1):
    """Train genus and species models with NCBI assemblies from the given genus."""

    if not isinstance(genus, str):
        raise TypeError("genus must be a string")

    # Check user input.
    genus = check_user_input(user_input=genus)

    # The directory name is defined in the following format: 'genus'_DD_MM_YYYY_hh-mm-ss
    curr_time = localtime()
    dir_name = f"{genus}_{curr_time[2]}_{curr_time[1]}_{curr_time[0]}_{curr_time[3]}-{curr_time[4]}-{curr_time[5]}"

    # Set the logger.
    set_logger(dir_name)

    # Time for the whole program.
    start_all = perf_counter()

    # Search for every defined species of the genus.
    start_tax = perf_counter()
    logger.info("Getting all species of the genus")
    children_ids = ncbi_children_tree.NCBIChildrenTree(genus).children_ids()
    species_dict = ncbi_taxon_metadata.NCBITaxonMetadata(children_ids).get_metadata()
    end_tax = perf_counter()

    # Get all gcf accessions that have Taxonomy check result OK.
    logger.info("Checking ANI data for updates")
    ani = html_scrap.TaxonomyCheck()
    ani_gcf = ani.ani_gcf()

    # Look for up to 8 assembly accessions per species.
    start_meta = perf_counter()
    logger.info("Getting assembly metadata")
    all_metadata = ncbi_assembly_metadata.NCBIAssemblyMetadata(
        all_metadata=species_dict, ani_gcf=ani_gcf, count=8, contig_n50=10000
    )
    all_metadata = all_metadata.get_all_metadata()
    logger.info("Finished metadata collecting\n")
    end_meta = perf_counter()

    # Ensure that the genus has at least one species with accessions.
    if not all_metadata:
        raise ValueError("No species with accessions found")

    # Download the chosen assemblies.
    # One file for each species with it's downloaded assemblies in zip format.
    start_download = perf_counter()

    # Iterate through all species.
    logger.info("Downloading assemblies for bloomfilter training")
    for metadata in all_metadata.values():
        # Only try to download when the species has accessions.
        if len(metadata["accessions"]) >= 1:
            sleep(5)
            species_name = metadata["sci_name"]
            tax_id = metadata["tax_id"]
            logger.info("Downloading {id}_{name}", id=tax_id, name=species_name)
            file_name = f"{tax_id}_{species_name}.zip"

            # Selecting the first 4 assemblies for training the filters.
            accessions = metadata["accessions"][:4]

            download_assemblies.download_assemblies(
                accessions=accessions,
                dir_name=dir_name,
                target_folder="zip_files",
                zip_file_name=file_name,
            )
    logger.info("Downloads finished\n")
    end_download = perf_counter()

    # Concatenate all assemblies of each species.
    start_concatenate = perf_counter()
    extract_bf = extract_and_concatenate.ExtractConcatenate(
        dir_name=dir_name, delete=True
    )
    extract_bf.bf()
    concatenate_meta(get_xspect_tmp_path() / dir_name, genus)
    logger.info("Finished extracting and concatenating\n")
    end_concatenate = perf_counter()

    # Download assemblies for svm creation.
    start_svm_dl = perf_counter()
    logger.info("Downloading assemblies for support-vector-machine training")
    accessions = {}
    for metadata in all_metadata.values():
        # Only add taxon with accessions.
        if len(metadata["accessions"]) >= 1:
            accessions[metadata["tax_id"]] = metadata["accessions"]

    # Downloading assemblies.
    create_svm.get_svm_assemblies(all_accessions=accessions, dir_name=dir_name)
    logger.info("Finished downloading\n")

    # Extracting assemblies.
    extract_svm = extract_and_concatenate.ExtractConcatenate(
        dir_name=dir_name, delete=True
    )
    extract_svm.svm(species_accessions=accessions)

    end_svm_dl = perf_counter()

    # Make dictionary for translating taxon ID to scientific name.
    translation_dict = create_translation_dict(dir_name)
    change_bf_assembly_file_names(dir_name)

    # Train new Bloomfilters with concatenated files of each species.
    start_bf = perf_counter()

    # Train Bloomfilters of species.
    logger.info("Training bloomfilters")
    species_files_path, species_result_path = interface_XspecT.make_paths(
        dir_name, genus
    )

    # Train Bloomfilter for complete genus.
    logger.info("Training metagenome bloomfilter")
    mg_files_path = Path(get_paths.get_current_dir_file_path(dir_name))

    genus_model = ProbabilisticSingleFilterModel(
        k=21,
        model_display_name=genus,
        author="Test",
        author_email="test@example.com",
        model_type="Genus",
        base_path=Path(species_result_path).parent,
    )
    genus_model.fit(mg_files_path / f"{genus}.fasta", genus)
    genus_model.save()

    # Delete metagenome file.
    os.remove(mg_files_path / f"{genus}.fasta")

    end_bf = perf_counter()

    # Create support vector machine.
    start_svm = perf_counter()
    logger.info("Training species filters with support-vector-machine")

    species_model = ProbabilisticFilterSVMModel(
        k=21,
        model_display_name=genus,
        author="Test",
        author_email="test@example.com",
        model_type="Species",
        base_path=Path(species_result_path).parent,
        kernel="rbf",
        c=1.0,
    )
    svm_dir = get_xspect_tmp_path() / dir_name / "training_data"
    species_model.fit(
        Path(species_files_path),
        svm_dir,
        display_names=translation_dict,
        svm_step=svm_step,
    )
    species_model.save()

    # Delete concatenated assemblies.
    # Delete species files.
    shutil.rmtree(species_files_path)

    # Delete used assemblies.
    assemblies_path = get_paths.get_current_dir_file_path(dir_name) / "training_data"
    shutil.rmtree(assemblies_path)

    end_svm = perf_counter()
    end_all = perf_counter()

    logger.info(
        "Program runtime: {time} m", time=(round((end_all - start_all) / 60, 2))
    )

    # Print and save collected statistics.
    logger.info("Saving collected runtime statistics")
    time_print = (
        f"Python version: {python_version()} \n"
        f"All time: {(end_all-start_all)/60:.2f} m\n"
        f"Tax time: {(end_tax-start_tax):.2f} s\n"
        f"Meta time: {(end_meta-start_meta)/60:.2f} m\n"
        f"Download time: {(end_download-start_download)/60:.2f} m\n"
        f"Concatenate time: {(end_concatenate-start_concatenate):.2f} s\n"
        f"Training time: {(end_bf-start_bf)/60:.2f} m\n"
        f"Support vector machine time: {((end_svm+end_svm_dl)-(start_svm+start_svm_dl))/60:.2f} m\n"
    )

    # Save time measurements.
    interface_XspecT.save_time_stats(time_print, dir_name)

    # Save translation dict
    # save_translation_dict(dir_name, translation_dict)

    logger.info("XspecT-trainer is finished.")


def train_from_directory(display_name: str, dir_path: Path, meta: bool = False):
    """Train the gene family and gene filter.

    :param display_name: Name of the model.
    :param dir: Input directory.
    """

    if not isinstance(display_name, str):
        raise TypeError("display_name must be a string")

    if not isinstance(dir_path, Path) and dir_path.exists() and dir_path.is_dir():
        raise ValueError("dir must be Path object to a valid directory")

    # check if the directory contains the necessary files
    # copy to temp path
    # check if svm training data exists
    # train model, with svm data if it exists
    # add display names
    # train metagenome model
    # clean up temp path
    pass
