"""
Training script tests
"""

import pytest
import src.xspect.train as train


def test_invalid_taxonomy_check():
    """
    Test if a ValueError is thrown when attempting to train a genus
    where species do not fulfill taxonomy check requirements.
    """
    with pytest.raises(ValueError):
        train.train("Amnimonas", "1", False, "", "", "")
