"""Scaffold smoke test for the pytest harness (Module 0).

Module 4 replaces the body with real BAP1-classifier assertions
(non-circular labels, CV + held-out AUROC).
"""
import numpy as np


def test_numpy_rng_is_reproducible():
    # Arrange / Act
    first = np.random.default_rng(0).integers(0, 10, size=3)
    second = np.random.default_rng(0).integers(0, 10, size=3)
    # Assert
    assert np.array_equal(first, second)
    assert first.shape == (3,)


def test_fixture_dir_exists(fixture_dir):
    # Assert the subsampled-fixture convention directory is present
    assert fixture_dir.is_dir()
