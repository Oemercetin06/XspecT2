"""Module for utility functions used in other modules regarding MLST."""

__author__ = "Cetin, Oemer"

import requests
import json
from io import StringIO
from pathlib import Path
from Bio import SeqIO
from xspect.definitions import get_xspect_model_path


def create_fasta_files(locus_path: Path, fasta_batch: str) -> None:
    """
    Create Fasta-Files for every allele of a locus.

    This function creates a fasta file for each record in the batch-string of a locus.
    The batch originates from an API-GET-request to PubMLST.
    The files are named after the record ID.
    If a fasta file already exists, it will be skipped.

    Args:
        locus_path (Path): The directory where the fasta-files will be saved.
        fasta_batch (str): A string containing every record of a locus from PubMLST.
    """
    # fasta_batch = full string of a fasta file containing every allele sequence of a locus
    for record in SeqIO.parse(StringIO(fasta_batch), "fasta"):
        number = record.id.split("_")[-1]  # example id = Oxf_cpn60_263
        output_fasta_file = locus_path / f"Allele_ID_{number}.fasta"
        if output_fasta_file.exists():
            continue  # Ignore existing ones
        with open(output_fasta_file, "w") as allele:
            SeqIO.write(record, allele, "fasta")


def pick_species_number_from_db(available_species: dict) -> str:
    """
    Get the chosen species from all available ones in the database.

    This function lists all available species of PubMLST.
    The user is then asked to pick a species by its associated number.

    Args:
        available_species (dict): A dictionary storing all available species.

    Returns:
        str: The name of the chosen species.

    Raises:
        ValueError: If the user input is not valid.
    """
    # The "database" string can look like this: pubmlst_abaumannii_seqdef
    for counter, database in available_species.items():
        print(str(counter) + ":" + database.split("_")[1])
    print("\nPick one of the above databases")
    while True:
        try:
            choice = input("Choose a species by selecting the corresponding number:")
            if int(choice) in available_species.keys():
                chosen_species = available_species.get(int(choice))
                return chosen_species
            else:
                print(
                    "Wrong input! Try again with a number that is available in the list above."
                )
        except ValueError:
            print(
                "Wrong input! Try again with a number that is available in the list above."
            )


def pick_scheme_number_from_db(available_schemes: dict) -> str:
    """
    Get the chosen scheme from all available ones of a species.

    This function lists all available schemes of a species.
    The user is then asked to pick a scheme by its associated number.

    Args:
        available_schemes (dict): A dictionary storing all available schemes.

    Returns:
        str: The name of the chosen scheme.

    Raises:
        ValueError: If the user input is not valid.
    """
    # List all available schemes of a species database
    for counter, scheme in available_schemes.items():
        print(str(counter) + ":" + scheme[0])
    print("\nPick any available scheme that is listed for download")
    while True:
        try:
            choice = input("Choose a scheme by selecting the corresponding number:")
            if int(choice) in available_schemes.keys():
                chosen_scheme = available_schemes.get(int(choice))[1]
                return chosen_scheme
            else:
                print(
                    "Wrong input! Try again with a number that is available in the above list."
                )
        except ValueError:
            print(
                "Wrong input! Try again with a number that is available in the above list."
            )


def scheme_list_to_dict(scheme_list: list[str]):
    """
    Converts the scheme list into a dictionary.

    Args:
        scheme_list (list[str]): A list storing all chosen schemes.

    Returns:
        dict: The converted dictionary.
    """
    return dict(zip(range(1, len(scheme_list) + 1), scheme_list))


def pick_scheme_from_models_dir() -> Path:
    """
    Get the chosen scheme from models that have been fitted prior.

    This function creates a dictionary containing all trained models.
    The dictionary is used as an argument for the "pick_scheme" function.

    Returns:
        Path: The path to the chosen model (trained).
    """
    schemes = {}
    counter = 1
    for entry in sorted((get_xspect_model_path() / "MLST").iterdir()):
        schemes[counter] = entry
        counter += 1
    return pick_scheme(schemes)


def pick_scheme(available_schemes: dict) -> Path:
    """
    Get the chosen scheme from the scheme dictionary.

    This function lists all available schemes of a species that have been downloaded.
    The user is then asked to pick a scheme by its associated number.

    Args:
        available_schemes (dict): A dictionary storing all available schemes.

    Returns:
        Path: The path to the chosen model (trained).

    Raises:
        ValueError: If the user input is not valid or if no scheme was downloaded prior.
    """
    if not available_schemes:
        raise ValueError("No scheme has been chosen for download yet!")

    if len(available_schemes.items()) == 1:
        return next(iter(available_schemes.values()))

    # List available schemes
    for counter, scheme in available_schemes.items():
        # For Strain Typing with an API-POST Request to the db
        if str(scheme).startswith("http"):
            scheme_json = requests.get(scheme).json()
            print(str(counter) + ":" + scheme_json["description"])

        # To pick a scheme after download for fitting
        else:
            print(str(counter) + ":" + str(scheme).split("/")[-1])

    print("\nPick a scheme for strain type prediction")
    while True:
        try:
            choice = input("Choose a scheme by selecting the corresponding number:")
            if int(choice) in available_schemes.keys():
                chosen_scheme = available_schemes.get(int(choice))
                return chosen_scheme
            else:
                print(
                    "Wrong input! Try again with a number that is available in the above list."
                )
        except ValueError:
            print(
                "Wrong input! Try again with a number that is available in the above list."
            )


class MlstResult:
    """Class for storing MLST results."""

    def __init__(
        self,
        scheme_model: str,
        steps: int,
        hits: dict[str, list[dict]],
        input_source: str = None,
    ):
        """Initialise an MlstResult object."""
        self.scheme_model = scheme_model
        self.steps = steps
        self.hits = hits
        self.input_source = input_source

    def get_results(self) -> dict:
        """
        Stores the result of a prediction in a dictionary.

        Returns:
            dict: The result dictionary with s sequence ID as key and the Strain type as value.
        """
        results = {seq_id: result for seq_id, result in self.hits.items()}
        return results

    def to_dict(self) -> dict:
        """
        Converts all attributes into one dictionary.

        Returns:
            dict: The dictionary containing all metadata of a run.
        """
        result = {
            "Scheme": self.scheme_model,
            "Steps": self.steps,
            "Results": self.get_results(),
            "Input_source": self.input_source,
        }
        return result

    def save(self, output_path: Path | str) -> None:
        """
        Saves the result as a JSON file.

        Args:
            output_path (Path,str): The path where the results are saved.
        """

        if isinstance(output_path, str):
            output_path = Path(output_path)

        output_path.parent.mkdir(exist_ok=True, parents=True)
        json_object = json.dumps(self.to_dict(), indent=4)

        with open(output_path, "w", encoding="utf-8") as file:
            file.write(json_object)
