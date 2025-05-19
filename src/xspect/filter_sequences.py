from pathlib import Path
from xspect.model_management import get_genus_model, get_species_model
from xspect.file_io import filter_sequences


def filter_species(
    model_genus: str,
    model_species: str,
    input_path: Path,
    output_path: Path,
    threshold: float,
):
    """Filter sequences by species.
    This function filters sequences from the input file based on the species model.
    It uses the genus model to identify the genus of the sequences and then applies
    the species model to filter the sequences.

    Args:
        model_genus (str): The genus model slug.
        model_species (str): The species model slug.
        input_path (Path): The path to the input file containing sequences.
        output_path (Path): The path to the output file where filtered sequences will be saved.
        threshold (float): The threshold for filtering sequences. Only sequences with a score
            above this threshold will be included in the output file.
    """
    species_model = get_species_model(model_genus)
    result = species_model.predict(input_path)
    included_ids = result.get_filtered_subsequence_labels(model_species, threshold)
    if not included_ids:
        print("No sequences found for the given species.")
        return
    filter_sequences(
        input_path,
        output_path,
        included_ids,
    )


def filter_genus(
    model_genus: str,
    input_path: Path,
    output_path: Path,
    threshold: float,
):
    genus_model = get_genus_model(model_genus)
    result = genus_model.predict(Path(input_path))
    included_ids = result.get_filtered_subsequence_labels(model_genus, threshold)
    if not included_ids:
        print("No sequences found for the given genus.")
        return

    filter_sequences(
        input_path,
        output_path,
        included_ids,
    )
