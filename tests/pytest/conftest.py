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


# Two synthetic populations x recognisable ccRCC markers + MT genes.
# Marker POSITIONS are load-bearing: `synthetic_10x_h5` boosts columns 3:7
# (tumour epithelial) and 7:11 (T cell), so never reorder or prepend here.
_FIXTURE_MARKER_GENES = [
    "MT-CO1", "MT-ND1", "MT-CYB",          # mitochondrial (drive pct_counts_mt)
    "CA9", "NDUFA4L2", "VEGFA", "EPCAM",   # tumour epithelial
    "CD3D", "CD3E", "CD8A", "IL7R",        # T cell
    "CD14", "LYZ", "CD68",                 # myeloid
    "PECAM1", "VWF", "CLDN5",              # endothelial
    "ACTA2", "COL1A1",                     # fibroblast
    "CD79A", "MS4A1",                      # B cell
]

# `process_singlecell` runs QC at PRODUCTION defaults (SC_MIN_GENES_PER_CELL =
# 200), so a fixture carrying only the 22 marker genes would lose every cell to
# `sc.pp.filter_cells` and collapse the matrix to 0x0. Pad with inert synthetic
# genes to 250 so the production QC thresholds are genuinely exercisable. The
# filler carries no biology and is deliberately NOT named after real symbols.
_FIXTURE_N_GENES = 250
_FIXTURE_GENES = _FIXTURE_MARKER_GENES + [
    f"SYNTH{i:04d}" for i in range(_FIXTURE_N_GENES - len(_FIXTURE_MARKER_GENES))
]


def _write_10x_v3_h5(path, dense_cells_by_genes, gene_names, barcodes):
    # h5py and scipy are imported HERE, not at module scope: conftest.py is
    # imported for the whole Python suite, and CI's py-checks job installs only
    # numpy/pandas/scikit-learn. A top-level `import h5py` would turn a missing
    # single-cell dependency into a collection error for every pytest file,
    # including the Module 4 BAP1 tests.
    import h5py
    import scipy.sparse as sp

    # sc.read_10x_h5 expects a CSC matrix of shape (n_genes, n_cells).
    mat = sp.csc_matrix(dense_cells_by_genes.T.astype(np.int32))
    n_genes, n_cells = mat.shape
    with h5py.File(path, "w") as f:
        g = f.create_group("matrix")
        g.create_dataset("data", data=mat.data)
        g.create_dataset("indices", data=mat.indices)
        g.create_dataset("indptr", data=mat.indptr)
        g.create_dataset("shape", data=np.array([n_genes, n_cells], dtype=np.int64))
        g.create_dataset("barcodes",
                         data=np.array([b.encode() for b in barcodes]))
        feat = g.create_group("features")
        ids = np.array([f"ENSG{i:011d}".encode() for i in range(n_genes)])
        names = np.array([n.encode() for n in gene_names])
        ftype = np.array([b"Gene Expression"] * n_genes)
        genome = np.array([b"GRCh38"] * n_genes)
        feat.create_dataset("id", data=ids)
        feat.create_dataset("name", data=names)
        feat.create_dataset("feature_type", data=ftype)
        feat.create_dataset("genome", data=genome)


@pytest.fixture(scope="session")
def synthetic_10x_h5(tmp_path_factory):
    rng = np.random.default_rng(0)
    n_cells, n_genes = 120, len(_FIXTURE_GENES)
    # baseline counts + a marker-driven two-population structure
    counts = rng.poisson(3, size=(n_cells, n_genes)).astype(np.int32)
    counts[:60, 3:7] += rng.poisson(30, size=(60, 4))   # pop A: tumour markers
    counts[60:, 7:11] += rng.poisson(30, size=(60, 4))  # pop B: T-cell markers
    counts[:, 0:3] += rng.poisson(2, size=(n_cells, 3))  # low MT => cells survive QC
    barcodes = [f"CELL{i:04d}-1" for i in range(n_cells)]
    path = str(tmp_path_factory.mktemp("sc") / "gse159115_synthetic.h5")
    _write_10x_v3_h5(path, counts, _FIXTURE_GENES, barcodes)
    return path
