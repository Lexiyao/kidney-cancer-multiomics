# kidney-cancer-multiomics

Reproducible TCGA-KIRC somatic multi-omics pipeline (ingest → MOFA2/SNF
integration → literature positive controls → low-dimensional survival + BAP1
classifier → live-updating Quarto/Plotly dashboard).

Rendered site: <https://lexiyao.github.io/kidney-cancer-multiomics>

## Status

Phases 0-1 complete and verified on the real data; Modules 2-6 not yet built.

Verified by real runs against the frozen `curatedTCGAData` KIRC snapshot inside
`bioconductor/bioconductor_docker:RELEASE_3_23`:

- **Container chain** (run `30570220145`): MOFA2 trains via
  `run_mofa(use_basilisk = FALSE)` against the system Python — no conda pulled.
- **Cohort census** (run `30642823359`, re-verified with zero drift by
  `30708943504`): main cohort **524** cases (RNA + methylation(HM27 ∪ HM450) +
  CNV); **413** with mutation; mutation MAF covers **417**. These are snapshot
  numbers and differ from the live GDC API — see the design spec §2.
- **Survival census** (run `30708943504`): **173** OS events (33.1%, median
  follow-up 1188 d; **148** within 5 years), so the EPV-10 predictor budget is
  17 (14 restricted). The pipeline computes this at fit time; Phase 4 wires 5.
- **Module 1 materialised** (run `30708943504`): `rna_mat` 5000 × 524,
  `methyl_mat` 5000 × 524, `cnv_mat` 24776 × 524, `mut_annot` 417 × 7.
- **`renv.lock`**: a complete 218-package lock snapshotted from a machine with
  every Import installed — not hand-authored.

**Not yet done:** MOFA2 integration (Module 2), the literature positive-control
suite (Module 3), the survival model and BAP1 classifier (Module 4), the
dashboard (Module 5) and single-cell (Module 6). No factor, subtype, survival
or discrimination result is claimed anywhere.

## Reproducibility scope (read before trusting the CI badge)

- The research core is a **frozen 2016 snapshot** (`curatedTCGAData` 1.34.0,
  snapshot `20160128`, hg19). It never updates.
- **CI does NOT run the full pipeline.** CI lints, runs unit tests on
  subsampled fixtures, and renders the dashboard from cached `_targets` /
  release-asset results. A green CI badge is **not** full reproduction.
- The full pipeline runs locally once (see `scripts/run_full_pipeline.R`);
  expected runtime and hardware are documented in `docs/runtime.md`.

## Skills mapping

_Populated in Module 5._

## Licence

MIT — see `LICENSE`.
