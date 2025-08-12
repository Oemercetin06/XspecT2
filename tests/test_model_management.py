"""Tests for the model_management module."""

import xspect.model_management as mm


def test_get_genus_model_path():
    """Test the get_genus_model function."""
    aci = mm.get_genus_model_path("Acinetobacter")
    sal = mm.get_genus_model_path("Salmonella")

    assert "acinetobacter-genus.json" in str(aci)
    assert "salmonella-genus.json" in str(sal)


def test_get_species_model_path():
    """Test the get_species_model path function."""
    aci = mm.get_species_model_path("Acinetobacter")
    sal = mm.get_species_model_path("Salmonella")

    assert "acinetobacter-species.json" in str(aci)
    assert "salmonella-species.json" in str(sal)


def test_get_model_metadata():
    """Test the get_model_metadata function."""
    aci_meta = mm.get_model_metadata("acinetobacter-genus")

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
