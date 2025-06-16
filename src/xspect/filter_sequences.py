from pathlib import Path
from xspect.model_management import get_genus_model, get_species_model
from xspect.file_io import filter_sequences, prepare_input_output_paths


def filter_species(
    model_genus: str,
    model_species: str,
    input_path: Path,
    output_path: Path,
    threshold: float,
    classification_output_path: Path | None = None,
    sparse_sampling_step: int = 1,
):
    """
    Filter sequences by species.

    This function filters sequences from the input file based on the species model.
    It uses the species model to identify the species of individual sequences and then applies
    a threshold filter the sequences.

    Args:
        model_genus (str): The genus of the species model.
        model_species (str): The species to filter by.
        input_path (Path): The path to the input file containing sequences.
        output_path (Path): The path to the output file where filtered sequences will be saved.
        classification_output_path (Path): Optional path to save the classification results.
        threshold (float): The threshold for filtering sequences. Only sequences with a score
            above this threshold will be included in the output file. A threshold of -1 will
            include only sequences if the species score is the highest among the
            available species scores.
        sparse_sampling_step (int): The step size for sparse sampling. Defaults to 1.
    """
    species_model = get_species_model(model_genus)
    input_paths, get_output_path = prepare_input_output_paths(input_path)

    for idx, current_path in enumerate(input_paths):
        result = species_model.predict(current_path, step=sparse_sampling_step)
        result.input_source = current_path.name

        if classification_output_path:
            cls_out = get_output_path(idx, classification_output_path)
            result.save(cls_out)
            print(
                f"Saved classification results from {current_path.name} as {cls_out.name}"
            )

        included_ids = result.get_filtered_subsequence_labels(model_species, threshold)
        if not included_ids:
            print(f"No sequences found for the given species in {current_path.name}.")
            continue

        filter_output_path = get_output_path(idx, output_path)
        filter_sequences(current_path, filter_output_path, included_ids)
        print(
            f"Saved filtered sequences from {current_path.name} as {filter_output_path.name}"
        )


def filter_genus(
    model_genus: str,
    input_path: Path,
    output_path: Path,
    threshold: float,
    classification_output_path: Path | None = None,
    sparse_sampling_step: int = 1,
):
    """
    Filter sequences by genus.

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
        sparse_sampling_step (int): The step size for sparse sampling. Defaults to 1.

    """
    model = get_genus_model(model_genus)
    input_paths, get_output_path = prepare_input_output_paths(input_path)

    for idx, current_path in enumerate(input_paths):
        result = model.predict(current_path, step=sparse_sampling_step)
        result.input_source = current_path.name

        if classification_output_path:
            cls_out = get_output_path(idx, classification_output_path)
            result.save(cls_out)
            print(
                f"Saved classification results from {current_path.name} as {cls_out.name}"
            )

        included_ids = result.get_filtered_subsequence_labels(model_genus, threshold)
        if not included_ids:
            print(f"No sequences found for the given genus in {current_path.name}.")
            continue

        filter_output_path = get_output_path(idx, output_path)
        filter_sequences(current_path, filter_output_path, included_ids)
        print(
            f"Saved filtered sequences from {current_path.name} as {filter_output_path.name}"
        )
