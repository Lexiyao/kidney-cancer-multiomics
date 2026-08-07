"""Shared pytest fixtures for the Python suite.

Heavy scanpy / mofapy2 paths are marked `heavy` and skipped in CI
(run only when HEAVY_PULL=true).
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

FIXTURE_DIR = Path(__file__).resolve().parents[1] / "fixtures"


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "heavy: triggers a heavy pull (mofapy2/scanpy); skipped in CI",
    )


@pytest.fixture
def fixture_dir() -> Path:
    return FIXTURE_DIR


@pytest.fixture
def heavy_pull_enabled() -> bool:
    return os.environ.get("HEAVY_PULL", "false").lower() == "true"


@pytest.fixture
def bap1_fixture():
    """Subsampled expression with a planted BAP1-linked signal (non-circular:
    labels are the mutation status, not derived from the features)."""
    rng = np.random.default_rng(20160128)
    n, g = 120, 50
    labels = np.array([1] * 40 + [0] * 80)
    base = rng.normal(size=(n, g))
    # genes 0..4 carry a mean shift only in BAP1-mutant samples
    base[labels == 1, :5] += 1.5
    expr = pd.DataFrame(base, columns=[f"gene{i}" for i in range(g)])
    return expr, labels


@pytest.fixture
def bap1_null_fixture():
    """Same shape as `bap1_fixture` but with NO planted signal: the labels are
    unrelated to the features. Used to prove the reported AUROC comes from the
    data rather than from the label leaking into the scoring rows."""
    rng = np.random.default_rng(1867)
    n, g = 120, 50
    labels = np.array([1] * 40 + [0] * 80)
    expr = pd.DataFrame(rng.normal(size=(n, g)),
                        columns=[f"gene{i}" for i in range(g)])
    return expr, labels
