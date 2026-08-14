"""Module 6 QC tests — run on a synthetic 10x H5, never on GSE159115.

`importorskip` guards the whole module so a host without the scanpy stack
skips instead of failing collection; ci.yml's py-checks installs
scanpy/anndata/h5py, so the skip LIFTS there and these tests run in CI.
No heavy pull is involved: the fixture is generated in-process.
"""
import pytest

pytest.importorskip("scanpy")
pytest.importorskip("h5py")

# Import deliberately below the importorskip guards, not at the top of file.
from python.singlecell_qc import (
    load_h5,
    load_h5ad,
    normalise,
    process_singlecell,
    run_qc,
)


def test_load_h5_reads_10x_matrix(synthetic_10x_h5):
    # Act
    adata = load_h5(synthetic_10x_h5)

    # Assert
    assert adata.n_obs == 120
    assert "CA9" in adata.var_names


def test_run_qc_flags_mito_and_filters(synthetic_10x_h5):
    # Arrange
    adata = load_h5(synthetic_10x_h5)

    # Act
    qc = run_qc(adata, min_genes=3, min_cells=1, max_pct_mt=90.0)

    # Assert: QC metrics attached, no cell exceeds the mito ceiling
    assert "pct_counts_mt" in qc.obs.columns
    assert (qc.obs["pct_counts_mt"] < 90.0).all()
    assert qc.n_obs <= adata.n_obs


def test_normalise_keeps_counts_and_flags_hvgs(synthetic_10x_h5):
    # Arrange
    qc = run_qc(load_h5(synthetic_10x_h5),
                min_genes=3, min_cells=1, max_pct_mt=90.0)

    # Act
    norm = normalise(qc)

    # Assert: raw counts preserved in a layer, log1p applied, HVGs flagged
    assert "counts" in norm.layers
    assert "log1p" in norm.uns
    assert "highly_variable" in norm.var.columns
    assert norm.var["highly_variable"].sum() > 0


def test_process_singlecell_returns_qcd_normalised_anndata(synthetic_10x_h5):
    # Act
    adata = process_singlecell(synthetic_10x_h5)

    # Assert: QC ran (mito metrics attached) and normalisation ran (counts
    # layer + log1p) on a non-empty matrix. There is deliberately NO cluster
    # column: the k-means step was removed because its output had no consumer
    # and its k no justification (see the module docstring).
    assert adata.n_obs > 0
    assert "pct_counts_mt" in adata.obs.columns
    assert "counts" in adata.layers
    assert "cluster" not in adata.obs.columns


def test_load_h5ad_reads_written_processed_object(synthetic_10x_h5, tmp_path):
    # Arrange: process then persist to .h5ad (mirrors the sc_object target)
    adata = process_singlecell(synthetic_10x_h5)
    out = str(tmp_path / "processed.h5ad")
    adata.write_h5ad(out)

    # Act: the 10x reader cannot parse this; load_h5ad must
    reloaded = load_h5ad(out)

    # Assert
    assert reloaded.n_obs == adata.n_obs
    assert "pct_counts_mt" in reloaded.obs.columns
    assert "counts" in reloaded.layers
