from pathlib import Path
from xspect.mlst_feature.mlst_helper import pick_scheme_from_models_dir
import xspect.model_management as mm
from xspect.models.probabilistic_filter_mlst_model import (
    ProbabilisticFilterMlstSchemeModel,
)
from xspect.definitions import fasta_endings, fastq_endings


def classify_genus(
    model_genus: str, input_path: Path, output_path: Path, step: int = 1
):
    """Classify the input file using the genus model."""
    model = mm.get_genus_model(model_genus)

    input_paths = []
    input_is_dir = input_path.is_dir()
    ending_wildcards = [f"*.{ending}" for ending in fasta_endings + fastq_endings]

    if input_is_dir:
        input_paths = [p for e in ending_wildcards for p in input_path.glob(e)]
    elif input_path.is_file():
        input_paths = [input_path]

    for idx, current_path in enumerate(input_paths):
        result = model.predict(current_path, step=step)
        result.input_source = current_path.name
        output_name = (
            f"{output_path.stem}_{idx+1}{output_path.suffix}"
            if input_is_dir
            else output_path.name
        )
        result.save(output_path.parent / output_name)
        print(f"Saved result as {output_name}")


def classify_species(model_genus, input_path, output_path, step=1):
    """Classify the input file using the species model."""
    model = mm.get_species_model(model_genus)

    input_paths = []
    input_is_dir = input_path.is_dir()
    ending_wildcards = [f"*.{ending}" for ending in fasta_endings + fastq_endings]

    if input_is_dir:
        input_paths = [p for e in ending_wildcards for p in input_path.glob(e)]
    elif input_path.is_file():
        input_paths = [input_path]

    for idx, current_path in enumerate(input_paths):
        result = model.predict(current_path, step=step)
        result.input_source = current_path.name
        output_name = (
            f"{output_path.stem}_{idx+1}{output_path.suffix}"
            if input_is_dir
            else output_path.name
        )
        result.save(output_path.parent / output_name)
        print(f"Saved result as {output_name}")


def classify_mlst(input_path, output_path):
    """Classify the input file using the MLST model."""
    scheme_path = pick_scheme_from_models_dir()
    model = ProbabilisticFilterMlstSchemeModel.load(scheme_path)
    result = model.predict(scheme_path, input_path)
    result.save(output_path)
