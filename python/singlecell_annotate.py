"""Module 6 (v1.1): cell-type annotation and GATED bulk->single-cell mapping.

map_bulk_signature must be called with purity_confounded taken from the
fn_subtype_purity_test gate (purity_check$is_purity_proxy). When the gate
failed, the mapping is re-framed: single-cell scores reflect cell-type
composition, not a tumour-cell-intrinsic subtype program (spec section 10).

Marker provenance (all canonical, none invented):
  Tumor_epithelial  CA9 / NDUFA4L2 are the two classic VHL-HIF-driven ccRCC
                    tumour-cell markers; VEGFA is a HIF target; EPCAM is a
                    generic epithelial marker. CAVEAT: EPCAM is frequently LOW
                    in ccRCC, and this panel has no proximal-tubule set, so
                    adjacent-normal tubule cells in GSE159115 can be pulled
                    into Tumor_epithelial. CAVEAT: NDUFA4L2 is also a
                    documented PERICYTE marker (Mesa-Ciller et al., J Cereb
                    Blood Flow Metab 2023, PMID 35929074); ccRCC is among the
                    most vascularised solid tumours, so pericyte signal is
                    split between this set and the ACTA2-containing
                    Fibroblast set -- one of them the tumour label. CAVEAT:
                    VEGFA is a hypoxia-response gene expressed broadly by
                    myeloid, stromal and tubular cells, not tumour-restricted,
                    so angiogenic macrophages are pulled toward this label.
                    See the limitation note below.
  T_cell            CD3D / CD3E (TCR-CD3 complex), CD8A (cytotoxic), IL7R
                    (memory/CD4; also marks ILCs).
  Myeloid           CD14 (monocyte/macrophage), LYZ (lysozyme), CD68.
  Endothelial       PECAM1 (CD31), VWF, CLDN5.
  Fibroblast        COL1A1 (collagen I); ACTA2 (alpha-SMA) also marks pericytes
                    and vascular smooth muscle, so this set is not pure.
  B_cell            CD79A, MS4A1 (CD20).

LIMITATION (must be stated wherever these labels are shown): assignment is
argmax over per-cell marker scores with no confidence threshold and no
proximal-tubule / NK / plasma-cell set, so every cell gets one of six labels
whether or not it belongs to any of them. The scores entering the argmax are
standardised against a per-set PERMUTATION NULL first: sc.tl.score_genes
returns mean(marker set) - mean(control genes), whose sampling variance
scales as 1/len(set), so RAW scores from panels of different size are not
comparable -- measured on pure noise, the 2-gene panels won the raw argmax
1.76x as often as the 4-gene ones, monotone in panel size. Each observed
score is therefore z-scored against SC_N_NULL_SETS random same-size gene
sets (per cell), which removes the size bias without penalising real signal
the way a plain across-cell z-score does. It does not add a confidence
threshold.
"""
from __future__ import annotations

import numpy as np
import scanpy as sc
from anndata import AnnData

CCRCC_MARKER_SETS = {
    "Tumor_epithelial": ["CA9", "NDUFA4L2", "VEGFA", "EPCAM"],
    "T_cell": ["CD3D", "CD3E", "CD8A", "IL7R"],
    "Myeloid": ["CD14", "LYZ", "CD68"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5"],
    "Fibroblast": ["ACTA2", "COL1A1"],
    "B_cell": ["CD79A", "MS4A1"],
}

# Policy floor, not a measured threshold: a signature scored from fewer than
# half its requested genes is flagged low_coverage in the returned coverage
# table. The exact genes-used / genes-requested counts are always reported
# beside every score, so the flag never substitutes for the number.
SC_MIN_SIGNATURE_COVERAGE = 0.5

# Permutation-null standardisation of marker scores (see the module
# docstring). 50 same-size random gene sets per panel: enough to estimate a
# per-cell null mean/sd, cheap enough that six panels stay interactive.
# Seeded so annotate_celltypes is deterministic.
SC_N_NULL_SETS = 50
SC_NULL_RANDOM_STATE = 0


def _null_standardised_score(adata: AnnData, genes: list, score_name: str,
                             rng: np.random.Generator,
                             n_null: int = SC_N_NULL_SETS) -> np.ndarray:
    """Score `genes`, then z-score per cell against same-size random gene sets.

    A plain across-cell z-score would also remove the size bias under the
    null, but it penalises REAL signal: a bimodal signal column has a large
    empirical sd, capping its z-range, while pure-noise columns reach higher
    z by chance (verified on the synthetic two-population fixture, where it
    halves recovery). The permutation null estimates the spread each panel
    SIZE produces by chance, per cell, and leaves genuine enrichment intact.
    """
    sc.tl.score_genes(adata, genes, score_name=score_name)
    observed = adata.obs[score_name].to_numpy(copy=True)
    pool = np.asarray(adata.var_names)
    null = np.empty((adata.n_obs, n_null))
    for j in range(n_null):
        draw = rng.choice(pool, size=len(genes), replace=False)
        sc.tl.score_genes(adata, list(draw), score_name="_score_null")
        null[:, j] = adata.obs["_score_null"].to_numpy()
    del adata.obs["_score_null"]
    spread = null.std(axis=1)
    spread[spread == 0.0] = 1.0
    z = (observed - null.mean(axis=1)) / spread
    adata.obs[score_name] = z
    return z


def annotate_celltypes(adata: AnnData, marker_sets: dict = CCRCC_MARKER_SETS) -> AnnData:
    adata = adata.copy()
    cell_types = list(marker_sets.keys())
    coverage = {}
    zscores = []
    rng = np.random.default_rng(SC_NULL_RANDOM_STATE)
    for cell_type in cell_types:
        col = f"score_{cell_type}"
        genes = marker_sets[cell_type]
        present = [g for g in genes if g in adata.var_names]
        if not present:
            # Refuse, never zero-fill: a panel with no measurable gene scored
            # as 0.0 enters a fabricated number into the argmax and silently
            # competes for (or concedes) every cell's label.
            raise ValueError(
                f"no marker gene of {cell_type} ({', '.join(genes)}) is "
                "present in adata.var_names: an unmeasurable panel cannot be "
                "scored. Drop the panel explicitly or fix the gene space; "
                "refusing to zero-fill."
            )
        zscores.append(_null_standardised_score(adata, present, col, rng))
        coverage[cell_type] = {
            "n_genes_requested": len(genes),
            "n_genes_used": len(present),
        }
    best = np.argmax(np.column_stack(zscores), axis=1)
    adata.obs["cell_type"] = [cell_types[i] for i in best]
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
    adata.uns["marker_coverage"] = coverage
    return adata


def _validate_gate_verdict(purity_confounded) -> bool:
    """Reject anything that is not an actual boolean gate verdict.

    `bool(None)` is False, so a missing or NA `purity_check$is_purity_proxy`
    coming across the reticulate boundary would otherwise be read as "the
    purity gate passed" and the mapping would be presented as cell-intrinsic.
    That is the unsafe direction, so fail loudly instead. numpy booleans are
    accepted because reticulate/numpy do not hand back builtin bools.
    """
    if not isinstance(purity_confounded, (bool, np.bool_)):
        raise TypeError(
            "purity_confounded must be a boolean taken from the purity gate "
            "(purity_check$is_purity_proxy); got "
            f"{type(purity_confounded).__name__}. Refusing to treat a missing "
            "or NA gate verdict as 'gate passed'."
        )
    return bool(purity_confounded)


def map_bulk_signature(adata: AnnData, signature_sets: dict,
                       purity_confounded: bool) -> dict:
    """Score bulk subtype signatures on annotated cells, gated by purity_check.

    Returns per-signature COVERAGE (n_genes_requested / n_genes_used /
    low_coverage) beside every score: the signatures are derived on bulk TCGA
    expression and intersected with the post-QC scRNA gene space, so partial
    dropout is the EXPECTED case, and a score computed from 2 of 20 genes must
    not print with the authority of one computed from all 20. A signature with
    NO gene present scores NaN, never 0.0 -- an exact zero is
    indistinguishable from a measured absence of enrichment.
    """
    purity_confounded = _validate_gate_verdict(purity_confounded)
    adata = adata.copy()
    scores_by_celltype = {}
    coverage = {}
    for subtype, genes in signature_sets.items():
        col = f"bulk_{subtype}_score"
        present = [g for g in genes if g in adata.var_names]
        n_requested = len(genes)
        used_fraction = len(present) / n_requested if n_requested else 0.0
        coverage[subtype] = {
            "n_genes_requested": n_requested,
            "n_genes_used": len(present),
            "low_coverage": used_fraction < SC_MIN_SIGNATURE_COVERAGE,
        }
        if present:
            sc.tl.score_genes(adata, present, score_name=col)
        else:
            adata.obs[col] = float("nan")
        scores_by_celltype[subtype] = (
            adata.obs.groupby("cell_type", observed=True)[col].mean().to_dict()
        )
    if purity_confounded:
        interpretation = (
            "Bulk subtypes are a tumour-purity / immune-infiltration proxy "
            "(purity_check failed): these single-cell scores reflect cell-type "
            "composition, NOT a tumour-cell-intrinsic subtype program."
        )
    else:
        # Deliberately narrow. Clearing one two-arm confound test does not
        # establish that a bulk signature scored on single cells is
        # tumour-cell-intrinsic, and the effect-size arm can stay silent on a
        # real but moderate gradient when the smallest subtype is small -- a
        # gate that did not fire is not a positive licence.
        interpretation = (
            "Purity confound gate did not fire: neither tumour purity nor "
            "immune infiltration separates the bulk subtypes at the tested "
            "effect size. That is NOT evidence the subtypes are "
            "tumour-cell-intrinsic; these single-cell scores remain "
            "EXPLORATORY."
        )
    return {
        "scores_by_celltype": scores_by_celltype,
        "coverage": coverage,
        "purity_confounded": purity_confounded,
        # The exploratory framing travels WITH the object, not only in the
        # dashboard prose around it: anything reading sc_mapping straight out
        # of the _targets store gets the caveat too.
        "status": "exploratory",
        "interpretation": interpretation,
    }
