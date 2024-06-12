"""Tests for the model_management module."""

import xspect.model_management as mm
from xspect.definitions import get_xspect_model_path


def test_get_genus_model():
    """Test the get_genus_model function."""
    aci = mm.get_genus_model("Acinetobacter")
    sal = mm.get_genus_model("Salmonella")

    assert aci.model_display_name == "Acinetobacter"
    assert aci.model_type == "Genus"

    assert sal.model_display_name == "Salmonella"
    assert sal.model_type == "Genus"


def test_get_species_model():
    """Test the get_species_model function."""
    aci = mm.get_species_model("Acinetobacter")
    sal = mm.get_species_model("Salmonella")

    assert aci.model_display_name == "Acinetobacter"
    assert aci.model_type == "Species"

    assert sal.model_display_name == "Salmonella"
    assert sal.model_type == "Species"


def test_get_model_metadata():
    """Test the get_model_metadata function."""
    aci_path = get_xspect_model_path() / "acinetobacter-genus.json"
    aci_meta = mm.get_model_metadata(aci_path)

    assert aci_meta["model_display_name"] == "Acinetobacter"


def test_get_models():
    """Test the get_models function."""
    model_dict = mm.get_models()
    assert "Genus" in model_dict
    assert "Species" in model_dict
    assert "Acinetobacter" in model_dict["Genus"]
    assert "Salmonella" in model_dict["Genus"]
    assert "Acinetobacter" in model_dict["Species"]
    assert "Salmonella" in model_dict["Species"]
