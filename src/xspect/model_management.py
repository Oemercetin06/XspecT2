""" This module contains functions to manage models. """

from json import loads
from pathlib import Path
from xspect.models.probabilistic_single_filter_model import (
    ProbabilisticSingleFilterModel,
)
from xspect.models.probabilistic_filter_svm_model import ProbabilisticFilterSVMModel
from xspect.definitions import get_xspect_model_path


def get_genus_model(genus):
    """Get a metagenomic model for the specified genus."""
    genus_model_path = get_xspect_model_path() / (genus.lower() + "-genus.json")
    genus_filter_model = ProbabilisticSingleFilterModel.load(genus_model_path)
    return genus_filter_model


def get_species_model(genus):
    """Get a species classification model for the specified genus."""
    species_model_path = get_xspect_model_path() / (genus.lower() + "-species.json")
    species_filter_model = ProbabilisticFilterSVMModel.load(species_model_path)
    return species_filter_model


def get_model_metadata(path: Path):
    """Get the metadata of a model."""
    if not isinstance(path, Path):
        raise TypeError("Path must be a pathlib.Path object.")
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f"Model file {path} does not exist.")

    with open(path, "r", encoding="utf-8") as file:
        model_json = loads(file.read())
        return model_json


def get_models():
    """Get a list of all available models in a dictionary by type."""
    model_dict = {}
    for model_file in get_xspect_model_path().glob("*.json"):
        model_metadata = get_model_metadata(model_file)
        model_type = model_metadata["model_type"]
        model_dict.setdefault(model_type, []).append(
            model_metadata["model_display_name"]
        )
    return model_dict
