from pathlib import Path
from xspect.mlst_feature.mlst_helper import pick_scheme_from_models_dir
import xspect.model_management as mm
from xspect.models.probabilistic_filter_mlst_model import (
    ProbabilisticFilterMlstSchemeModel,
)


def classify_genus(
    model_genus: str, input_path: Path, output_path: Path, step: int = 1
):
    """Classify the input file using the genus model."""
    model = mm.get_genus_model(model_genus)
    result = model.predict(input_path, step=step)
    result.input_source = input_path.name
    result.save(output_path)


def classify_species(model_genus, input_path, output_path, step=1):
    """Classify the input file using the species model."""
    model = mm.get_species_model(model_genus)
    result = model.predict(input_path, step=step)
    result.input_source = input_path.name
    result.save(output_path)


def classify_mlst(input_path, output_path):
    """Classify the input file using the MLST model."""
    scheme_path = pick_scheme_from_models_dir()
    model = ProbabilisticFilterMlstSchemeModel.load(scheme_path)
    result = model.predict(scheme_path, input_path)
    result.save(output_path)
