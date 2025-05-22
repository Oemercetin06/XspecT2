from pathlib import Path
from xspect.model_management import get_genus_model, get_species_model
from xspect.file_io import filter_sequences
from xspect.definitions import fasta_endings, fastq_endings


def filter_species(
    model_genus: str,
    model_species: str,
    input_path: Path,
    output_path: Path,
    threshold: float,
    classification_output_path: Path | None = None,
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
            above this threshold will be included in the output file. A threshold of -1 will
            include only sequences if the species score is the highest among the
            available species scores.
        classification_output_path (Path): Optional path to save the classification results.
        threshold (float): The threshold for filtering sequences. Only sequences with a score
            above this threshold will be included in the output file. A threshold of -1 will
            include only sequences if the species score is the highest among the
            available species scores.
    """
    species_model = get_species_model(model_genus)

    input_paths = []
    input_is_dir = input_path.is_dir()
    ending_wildcards = [f"*.{ending}" for ending in fasta_endings + fastq_endings]

    if input_is_dir:
        input_paths = [p for e in ending_wildcards for p in input_path.glob(e)]
    elif input_path.is_file():
        input_paths = [input_path]

    for idx, current_path in enumerate(input_paths):
        result = species_model.predict(current_path)
        result.input_source = current_path.name

        if classification_output_path:
            classification_output_name = (
                f"{classification_output_path.stem}_{idx+1}{classification_output_path.suffix}"
                if input_is_dir
                else classification_output_path.name
            )
            result.save(classification_output_path.parent / classification_output_name)
            print(
                f"Saved classification results from {current_path.name} as {classification_output_name}"
            )

        included_ids = result.get_filtered_subsequence_labels(model_species, threshold)
        if not included_ids:
            print(f"No sequences found for the given species in {current_path.name}.")
            continue
        output_name = (
            f"{output_path.stem}_{idx+1}{output_path.suffix}"
            if input_is_dir
            else output_path.name
        )
        filter_sequences(
            current_path,
            output_path.parent / output_name,
            included_ids,
        )
        print(f"Saved filtered sequences from {current_path.name} as {output_name}")


def filter_genus(
    model_genus: str,
    input_path: Path,
    output_path: Path,
    threshold: float,
    classification_output_path: Path | None = None,
):
    """Filter sequences by genus.
    This function filters sequences from the input file based on the genus model.
    It uses the genus model to identify the genus of the sequences and then applies
    the filtering based on the provided threshold.

    Args:
        model_genus (str): The genus model slug.
        input_path (Path): The path to the input file containing sequences.
        output_path (Path): The path to the output file where filtered sequences will be saved.
        threshold (float): The threshold for filtering sequences. Only sequences with a score
            above this threshold will be included in the output file.
        classification_output_path (Path): Optional path to save the classification results.

    """
    genus_model = get_genus_model(model_genus)

    input_paths = []
    input_is_dir = input_path.is_dir()
    ending_wildcards = [f"*.{ending}" for ending in fasta_endings + fastq_endings]

    if input_is_dir:
        input_paths = [p for e in ending_wildcards for p in input_path.glob(e)]
    elif input_path.is_file():
        input_paths = [input_path]

    for idx, current_path in enumerate(input_paths):
        result = genus_model.predict(current_path)
        result.input_source = current_path.name

        if classification_output_path:
            classification_output_name = (
                f"{classification_output_path.stem}_{idx+1}{classification_output_path.suffix}"
                if input_is_dir
                else classification_output_path.name
            )
            result.save(classification_output_path.parent / classification_output_name)
            print(
                f"Saved classification results from {current_path.name} as {classification_output_name}"
            )

        included_ids = result.get_filtered_subsequence_labels(model_genus, threshold)
        if not included_ids:
            print(f"No sequences found for the given genus in {current_path.name}.")
            continue
        output_name = (
            f"{output_path.stem}_{idx+1}{output_path.suffix}"
            if input_is_dir
            else output_path.name
        )
        filter_sequences(
            current_path,
            output_path.parent / output_name,
            included_ids,
        )
        print(f"Saved filtered sequences from {current_path.name} as {output_name}")
