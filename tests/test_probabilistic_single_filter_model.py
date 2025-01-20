"""Tests for the ProbabilisticFilterModel class."""

# pylint: disable=redefined-outer-name, line-too-long, protected-access

from pathlib import Path
import pytest
from Bio.Seq import Seq
from xspect.models.probabilistic_single_filter_model import (
    ProbabilisticSingleFilterModel,
)


@pytest.fixture
def filter_model(tmpdir):
    """Fixture for the ProbabilisticFilterModel class."""
    base_path = tmpdir.mkdir("xspect_data")
    return ProbabilisticSingleFilterModel(
        k=21,
        model_display_name="Test Filter",
        author="John Doe",
        author_email="john.doe@example.com",
        model_type="Species",
        base_path=Path(base_path),
    )


@pytest.fixture
def trained_filter_model(filter_model, concatenated_assembly_file_path):
    """Fixture for the ProbabilisticFilterModel class."""
    filter_model.fit(Path(concatenated_assembly_file_path), "Acinetobacter baumannii")
    return filter_model


def test_fit(trained_filter_model):
    """Test the fit method."""
    assert trained_filter_model.bf is not None
    assert trained_filter_model.display_names == {
        "concatenated_assembly": "Acinetobacter baumannii"
    }


def test_calculate_hits(trained_filter_model):
    """Test the calculate_hits method."""
    hits = trained_filter_model.calculate_hits(Seq("TAAATAAATTTATATAGCTAAA"))
    assert hits == {"concatenated_assembly": 2}


def test_save_and_load(trained_filter_model):
    """Test the save and load methods."""
    trained_filter_model.save()
    loaded_filter_model = ProbabilisticSingleFilterModel.load(
        trained_filter_model.base_path / (trained_filter_model.slug() + ".json")
    )
    assert loaded_filter_model.to_dict() == trained_filter_model.to_dict()
    hits = trained_filter_model.calculate_hits(Seq("TAAATAAATTTATATAGCTAAA"))
    assert hits == {"concatenated_assembly": 2}
