from pathlib import Path
from xspect.mlst_feature.mlst_helper import pick_scheme_from_models_dir
import xspect.model_management as mm
from xspect.models.probabilistic_filter_mlst_model import (
    ProbabilisticFilterMlstSchemeModel,
)
from xspect.file_io import prepare_input_output_paths


def classify_genus(
    model_genus: str, input_path: Path, output_path: Path, step: int = 1
):
    """
    Classify the genus of sequences.

    This function classifies input files using the genus model.
    The input path can be a file or directory

    Args:
        model_genus (str): The genus model slug.
        input_path (Path): The path to the input file/directory containing sequences.
        output_path (Path): The path to the output file where results will be saved.
        step (int): The amount of kmers to be skipped.
    """
    model = mm.get_genus_model(model_genus)
    input_paths, get_output_path = prepare_input_output_paths(input_path)

    for idx, current_path in enumerate(input_paths):
        result = model.predict(current_path, step=step)
        result.input_source = current_path.name
        cls_path = get_output_path(idx, output_path)
        result.save(cls_path)
        print(f"Saved result as {cls_path.name}")


def classify_species(
    model_genus: str, input_path: Path, output_path: Path, step: int = 1
):
    """
    Classify the species of sequences.

    This function classifies input files using the species model.
    The input path can be a file or directory

    Args:
        model_genus (str): The genus model slug.
        input_path (Path): The path to the input file/directory containing sequences.
        output_path (Path): The path to the output file where results will be saved.
        step (int): The amount of kmers to be skipped.
    """
    model = mm.get_species_model(model_genus)
    input_paths, get_output_path = prepare_input_output_paths(input_path)

    for idx, current_path in enumerate(input_paths):
        result = model.predict(current_path, step=step)
        result.input_source = current_path.name
        cls_path = get_output_path(idx, output_path)
        result.save(cls_path)
        print(f"Saved result as {cls_path.name}")


def classify_mlst(input_path: Path, output_path: Path, limit: bool):
    """
    Classify the strain type using the specific MLST model.

    Args:
        input_path (Path): The path to the input file/directory containing sequences.
        output_path (Path): The path to the output file where results will be saved.
        limit (bool): A limit for the highest allele_id results that are shown.
    """

    scheme_path = pick_scheme_from_models_dir()
    model = ProbabilisticFilterMlstSchemeModel.load(scheme_path)
    input_paths, get_output_path = prepare_input_output_paths(input_path)
    for idx, current_path in enumerate(input_paths):
        result = model.predict(scheme_path, current_path, step=1, limit=limit)
        result.input_source = current_path.name
        cls_path = get_output_path(idx, output_path)
        result.save(cls_path)
        print(f"Saved result as {cls_path.name}")
