"""Shared pytest fixtures for the Python suite.

Heavy scanpy / mofapy2 paths are marked `heavy` and skipped in CI
(run only when HEAVY_PULL=true).
"""
import os
from pathlib import Path

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
