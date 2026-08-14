"""Module 6 (v1.1): scanpy QC and normalisation for GSE159115.

Dataset: GSE159115 (Zhang et al., PNAS 2021, PMID 34099557,
DOI 10.1073/pnas.2103240118) -- per-sample 10x H5 files; this module
processes ONE specimen's H5 (see config/params.yml).

NO CLUSTERING STEP, deliberately. An earlier revision ran sklearn KMeans at a
hard-coded k = 10 on the PCA embedding, but the resulting obs['cluster']
column had NO consumer -- sc_mapping groups by cell_type, and the dashboard
page has no cluster section -- and the k had no selection criterion, no
stability check and no documented rationale (on the synthetic two-population
fixture it produced 10 clusters including singletons). A partition nobody
reads, computed with an arbitrary parameter, is exactly the kind of
unexamined number this project refuses elsewhere, so it was removed rather
than justified after the fact. If clustering is ever needed, it must arrive
with a consumer and a stated k-selection rule.

All functions return a new AnnData; inputs are never mutated in place.
"""
from __future__ import annotations

import scanpy as sc
from anndata import AnnData

SC_MIN_GENES_PER_CELL = 200
SC_MIN_CELLS_PER_GENE = 3
SC_MAX_PCT_MT = 20.0
SC_N_TOP_HVG = 2000
SC_TARGET_SUM = 1e4


def load_h5(path: str) -> AnnData:
    adata = sc.read_10x_h5(path)
    adata.var_names_make_unique()
    return adata


def load_h5ad(path: str) -> AnnData:
    """Read a written AnnData .h5ad (e.g. the processed sc_object).

    sc.read_10x_h5 only parses 10x CellRanger HDF5, not an AnnData .h5ad,
    so re-reading the processed object requires sc.read_h5ad.
    """
    return sc.read_h5ad(path)


def run_qc(adata: AnnData,
           min_genes: int = SC_MIN_GENES_PER_CELL,
           min_cells: int = SC_MIN_CELLS_PER_GENE,
           max_pct_mt: float = SC_MAX_PCT_MT) -> AnnData:
    adata = adata.copy()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                               log1p=False, inplace=True)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    return adata[adata.obs["pct_counts_mt"] < max_pct_mt, :].copy()


def normalise(adata: AnnData,
              n_top_hvg: int = SC_N_TOP_HVG,
              target_sum: float = SC_TARGET_SUM) -> AnnData:
    adata = adata.copy()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=min(n_top_hvg, adata.n_vars))
    return adata


def process_singlecell(path: str) -> AnnData:
    return normalise(run_qc(load_h5(path)))
