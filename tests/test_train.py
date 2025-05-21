"""
Training script tests
"""

import os
import pytest
from xspect import train
from xspect.definitions import get_xspect_model_path


def test_invalid_taxonomy_check():
    """
    Test if a ValueError is thrown when attempting to train a genus
    where species do not fulfill taxonomy check requirements.
    """
    ncbi_api_key = os.environ.get("NCBI_API_KEY")

    with pytest.raises(ValueError):
        train.train_from_ncbi("Amnimonas", ncbi_api_key=ncbi_api_key)


def test_train_from_ncbi():
    """
    Test the train_from_ncbi function with a valid genus name.
    """
    ncbi_api_key = os.environ.get("NCBI_API_KEY")

    train.train_from_ncbi("Salmonella", ncbi_api_key=ncbi_api_key)
    model_path = get_xspect_model_path()
    salmonella_paths = [
        model_path / "salmonella-species.json",
        model_path / "salmonella-genus.json",
        model_path / "salmonella-species" / "scores.csv",
        model_path / "salmonella-species" / "index.cobs_classic",
        model_path / "salmonella-genus" / "filter.bloom",
    ]
    for salmonella_path in salmonella_paths:
        assert salmonella_path.exists()
