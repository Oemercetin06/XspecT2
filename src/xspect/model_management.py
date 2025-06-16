"""This module contains functions to manage models."""

from json import loads, dumps
from pathlib import Path
from xspect.models.probabilistic_single_filter_model import (
    ProbabilisticSingleFilterModel,
)
from xspect.models.probabilistic_filter_svm_model import ProbabilisticFilterSVMModel
from xspect.definitions import get_xspect_model_path


def get_genus_model(genus) -> ProbabilisticSingleFilterModel:
    """
    Get a genus model for the specified genus.

    This function retrieves a pre-trained genus classification model based on the provided genus name.

    Args:
        genus (str): The genus name for which the model is to be retrieved.

    Returns:
        ProbabilisticSingleFilterModel: An instance of the genus classification model.
    """
    genus_model_path = get_xspect_model_path() / (genus.lower() + "-genus.json")
    genus_filter_model = ProbabilisticSingleFilterModel.load(genus_model_path)
    return genus_filter_model


def get_species_model(genus) -> ProbabilisticFilterSVMModel:
    """
    Get a species classification model for the specified genus.

    This function retrieves a pre-trained species classification model based on the provided genus name.

    Args:
        genus (str): The genus name for which the species model is to be retrieved.

    Returns:
        ProbabilisticFilterSVMModel: An instance of the species classification model.
    """
    species_model_path = get_xspect_model_path() / (genus.lower() + "-species.json")
    species_filter_model = ProbabilisticFilterSVMModel.load(species_model_path)
    return species_filter_model


def get_model_metadata(model: str | Path) -> dict:
    """
    Get metadata of a specified model.

    This function retrieves the metadata of a model from its JSON file.

    Args:
        model (str | Path): The slug of the model (as a string) or the path to the model JSON file.

    Returns:
        dict: A dictionary containing the model metadata.

    Raises:
        ValueError: If the model does not exist or is not a valid file.
    """
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


def update_model_metadata(model_slug: str, author: str, author_email: str) -> None:
    """
    Update the metadata of a model.

    This function updates the author and author email in the model's metadata JSON file.

    Args:
        model_slug (str): The slug of the model to update.
        author (str): The name of the author to set in the metadata.
        author_email (str): The email of the author to set in the metadata.
    """
    model_metadata = get_model_metadata(model_slug)
    model_metadata["author"] = author
    model_metadata["author_email"] = author_email

    model_path = get_xspect_model_path() / (model_slug + ".json")
    with open(model_path, "w", encoding="utf-8") as file:
        file.write(dumps(model_metadata, indent=4))


def update_model_display_name(
    model_slug: str, filter_id: str, display_name: str
) -> None:
    """
    Update the display name of a filter in a model.

    This function updates the display name of a specific filter in the model's metadata JSON file.

    Args:
        model_slug (str): The slug of the model to update.
        filter_id (str): The ID of the filter whose display name is to be updated.
        display_name (str): The new display name for the filter.
    """
    model_metadata = get_model_metadata(model_slug)
    model_metadata["display_names"][filter_id] = display_name

    model_path = get_xspect_model_path() / (model_slug + ".json")
    with open(model_path, "w", encoding="utf-8") as file:
        file.write(dumps(model_metadata, indent=4))


def get_models() -> dict[str, list[dict]]:
    """
    Get a list of all available models in a dictionary by type.

    This function scans the model directory for JSON files and organizes them by their model type.

    Returns:
        dict[str, list[dict]]: A dictionary where keys are model types and values are lists of model display names.
    """
    model_dict = {}
    for model_file in get_xspect_model_path().glob("*.json"):
        model_metadata = get_model_metadata(model_file)
        model_type = model_metadata["model_type"]
        model_dict.setdefault(model_type, []).append(
            model_metadata["model_display_name"]
        )
    return model_dict


def get_model_display_names(model_slug: str) -> list[str]:
    """
    Get the display names included in a model.

    This function retrieves the display names of individual filters from the model's metadata.

    Args:
        model_slug (str): The slug of the model for which to retrieve display names.

    Returns:
        list[str]: A list of display names for the individual filters in the model.
    """
    model_metadata = get_model_metadata(model_slug)
    return list(model_metadata["display_names"].values())
