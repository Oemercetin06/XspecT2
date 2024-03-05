"""
Training script tests
"""

import pytest
import src.xspect.train as train


def test_invalid_taxonomy_check():
    with pytest.raises(ValueError):
        train.train("Amnimonas", "1", False, "", "", "")
