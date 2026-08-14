"""Module 6 annotation + gated-mapping tests — synthetic 10x H5 only.

Same `importorskip` gate as test_singlecell_qc.py: CI's py-checks job has no
scanpy/h5py, so this module must SKIP there rather than fail collection.
"""
import numpy as np
import pytest

pytest.importorskip("scanpy")
pytest.importorskip("h5py")

# Imports deliberately below the importorskip guards, not at the top of file.
from python.singlecell_annotate import annotate_celltypes, map_bulk_signature
from python.singlecell_qc import load_h5, normalise, run_qc


def _prep(path):
    return normalise(run_qc(load_h5(path),
                            min_genes=3, min_cells=1, max_pct_mt=90.0))


def test_annotate_assigns_a_celltype_to_every_cell(synthetic_10x_h5):
    # Arrange
    adata = _prep(synthetic_10x_h5)

    # Act
    annotated = annotate_celltypes(adata)

    # Assert: every cell labelled, and BOTH planted populations appear. The
    # old assertion ("Tumor_epithelial in set") was satisfiable by a constant
    # labeller that named every one of the 120 cells Tumor_epithelial.
    assert "cell_type" in annotated.obs.columns
    assert annotated.obs["cell_type"].notna().all()
    assert {"Tumor_epithelial", "T_cell"} <= set(annotated.obs["cell_type"])


def test_annotate_recovers_both_planted_populations(synthetic_10x_h5):
    """The fixture plants tumour markers in cells 0:60 and T-cell markers in
    cells 60:120 (conftest.synthetic_10x_h5). The annotator must recover that
    structure for a substantial majority of each population -- a regression
    that collapses the labels into one class, or a scoring change that lets
    panel size rather than signal decide the argmax, fails here."""
    # Arrange
    adata = _prep(synthetic_10x_h5)

    # Act
    annotated = annotate_celltypes(adata)

    # Assert: barcodes are CELL<i:04d>-1, so the planted population is
    # recoverable from the obs name even after QC filtering.
    idx = np.array([int(name[4:8]) for name in annotated.obs_names])
    labels = np.asarray(annotated.obs["cell_type"])
    pop_a = labels[idx < 60]
    pop_b = labels[idx >= 60]
    assert len(pop_a) > 0 and len(pop_b) > 0
    assert (pop_a == "Tumor_epithelial").mean() >= 0.9
    assert (pop_b == "T_cell").mean() >= 0.9


def test_annotate_refuses_a_marker_set_with_no_measurable_gene(synthetic_10x_h5):
    """A panel with zero genes in var_names must RAISE, not score 0.0: a
    zero-filled column enters a fabricated number into the argmax."""
    # Arrange
    adata = _prep(synthetic_10x_h5)
    sets = {"Ghost": ["NOT_A_GENE_A", "NOT_A_GENE_B"],
            "T_cell": ["CD3D", "CD3E"]}

    # Act / Assert
    with pytest.raises(ValueError, match="Ghost"):
        annotate_celltypes(adata, marker_sets=sets)


def test_annotate_records_marker_coverage(synthetic_10x_h5):
    # Arrange
    adata = _prep(synthetic_10x_h5)

    # Act
    annotated = annotate_celltypes(adata)

    # Assert: genes-used / genes-requested travels with the object
    cov = annotated.uns["marker_coverage"]
    assert cov["Tumor_epithelial"] == {"n_genes_requested": 4,
                                       "n_genes_used": 4}
    assert set(cov) == {"Tumor_epithelial", "T_cell", "Myeloid",
                        "Endothelial", "Fibroblast", "B_cell"}


def test_map_bulk_signature_reframes_when_purity_confounded(synthetic_10x_h5):
    # Arrange
    adata = annotate_celltypes(_prep(synthetic_10x_h5))
    sigs = {"ccA": ["CA9", "EPCAM"], "ccB": ["VWF", "PECAM1"]}

    # Act
    res = map_bulk_signature(adata, sigs, purity_confounded=True)

    # Assert
    assert res["purity_confounded"] is True
    assert "proxy" in res["interpretation"].lower()
    assert set(res["scores_by_celltype"].keys()) == {"ccA", "ccB"}


@pytest.mark.parametrize("bad_verdict", [None, 1, 0, "TRUE", [True]])
def test_map_bulk_signature_refuses_a_non_boolean_gate_verdict(synthetic_10x_h5,
                                                               bad_verdict):
    """A missing/NA gate verdict must NOT silently read as 'gate passed'.

    Task 6.8 hands `purity_check$is_purity_proxy` across the reticulate
    boundary; an absent list element arrives as None and `bool(None)` is False,
    i.e. exactly the unsafe direction (mapping presented as cell-intrinsic).
    """
    # Arrange
    adata = annotate_celltypes(_prep(synthetic_10x_h5))
    sigs = {"ccA": ["CA9", "EPCAM"]}

    # Act / Assert
    with pytest.raises(TypeError, match="purity_confounded"):
        map_bulk_signature(adata, sigs, purity_confounded=bad_verdict)


def test_map_bulk_signature_accepts_a_numpy_bool_gate_verdict(synthetic_10x_h5):
    """reticulate/numpy hand a np.bool_, not a builtin bool — must be allowed."""
    # Arrange
    adata = annotate_celltypes(_prep(synthetic_10x_h5))
    sigs = {"ccA": ["CA9", "EPCAM"]}

    # Act
    res = map_bulk_signature(adata, sigs, purity_confounded=np.bool_(True))

    # Assert
    assert res["purity_confounded"] is True
    assert "proxy" in res["interpretation"].lower()


def test_map_bulk_signature_states_a_non_firing_gate_narrowly(synthetic_10x_h5):
    """A gate that did not fire is NOT a positive licence. The old pass-branch
    text ('interpretable as cell-intrinsic') overclaimed on a test that only
    establishes purity/immune infiltration did not reach the stated effect
    size; the reworded branch must say what was tested and keep the scores
    exploratory."""
    # Arrange
    adata = annotate_celltypes(_prep(synthetic_10x_h5))
    sigs = {"ccA": ["CA9", "EPCAM"], "ccB": ["VWF", "PECAM1"]}

    # Act
    res = map_bulk_signature(adata, sigs, purity_confounded=False)

    # Assert
    assert res["purity_confounded"] is False
    assert "did not fire" in res["interpretation"].lower()
    assert "exploratory" in res["interpretation"].lower()
    assert "not evidence" in res["interpretation"].lower()
    assert "interpretable as cell-intrinsic" not in res["interpretation"].lower()
    # The exploratory framing travels WITH the object, not only in prose.
    assert res["status"] == "exploratory"


def test_map_bulk_signature_reports_nan_and_coverage_for_missing_genes(synthetic_10x_h5):
    """An unmeasurable signature must be NaN, never an exact 0.0 -- kable
    prints 0.000 identically for 'measured no enrichment' and 'no genes
    matched'. Coverage counts must accompany every score so a partial overlap
    cannot print with full-panel authority."""
    # Arrange
    adata = annotate_celltypes(_prep(synthetic_10x_h5))
    sigs = {"ccX": ["NOT_A_GENE_A", "NOT_A_GENE_B"],
            "ccA": ["CA9", "EPCAM"]}

    # Act
    res = map_bulk_signature(adata, sigs, purity_confounded=False)

    # Assert: the unmeasurable row is all-NaN, the measurable one is not
    ccx = list(res["scores_by_celltype"]["ccX"].values())
    cca = list(res["scores_by_celltype"]["ccA"].values())
    assert all(np.isnan(v) for v in ccx)
    assert not any(np.isnan(v) for v in cca)
    assert res["coverage"]["ccX"] == {"n_genes_requested": 2,
                                      "n_genes_used": 0,
                                      "low_coverage": True}
    assert res["coverage"]["ccA"] == {"n_genes_requested": 2,
                                      "n_genes_used": 2,
                                      "low_coverage": False}
