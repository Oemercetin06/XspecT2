"""This module contains functions to manage models."""

from json import loads, dumps
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


def get_model_metadata(model: str | Path):
    """Get the metadata of a model."""
    if isinstance(model, str):
        model_path = get_xspect_model_path() / (model.lower() + ".json")
    elif isinstance(model, Path):
        model_path = model
    else:
        raise ValueError("Model must be a string (slug) or a Path object.")

    if not model_path.exists() or not model_path.is_file():
        raise ValueError(f"Model at {model_path} does not exist.")

    with open(model_path, "r", encoding="utf-8") as file:
        model_json = loads(file.read())
        return model_json


def update_model_metadata(model_slug: str, author: str, author_email: str):
    """Update the metadata of a model."""
    model_metadata = get_model_metadata(model_slug)
    model_metadata["author"] = author
    model_metadata["author_email"] = author_email

    model_path = get_xspect_model_path() / (model_slug + ".json")
    with open(model_path, "w", encoding="utf-8") as file:
        file.write(dumps(model_metadata, indent=4))


def update_model_display_name(model_slug: str, filter_id: str, display_name: str):
    """Update the display name of a filter in a model."""
    model_metadata = get_model_metadata(model_slug)
    model_metadata["display_names"][filter_id] = display_name

    model_path = get_xspect_model_path() / (model_slug + ".json")
    with open(model_path, "w", encoding="utf-8") as file:
        file.write(dumps(model_metadata, indent=4))


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


def get_model_display_names(model_slug: str):
    """Get the display names included in a model."""
    model_metadata = get_model_metadata(model_slug)
    return list(model_metadata["display_names"].values())
