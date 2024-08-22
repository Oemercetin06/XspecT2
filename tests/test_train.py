"""
Training script tests
"""

import pytest
from src.xspect import train


def test_invalid_taxonomy_check():
    """
    Test if a ValueError is thrown when attempting to train a genus
    where species do not fulfill taxonomy check requirements.
    """
    with pytest.raises(ValueError):
        train.train_ncbi("Amnimonas")
