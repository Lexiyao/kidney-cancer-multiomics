# kidney-cancer-multiomics

Reproducible TCGA-KIRC somatic multi-omics pipeline (ingest → MOFA2/SNF
integration → literature positive controls → low-dimensional survival + BAP1
classifier → live-updating Quarto/Plotly dashboard).

Rendered site: <https://lexiyao.github.io/kidney-cancer-multiomics>

## Status

Phases 0-3 have run on the real data. Module 4 (survival models + BAP1
classifier) is **implemented and unit-tested on fixtures but has never been
executed against the frozen snapshot**; Modules 5-6 are not built. **No
survival model has been fitted, so no C-index, calibration, hazard ratio or
discrimination result is claimed anywhere.**

**Two of the Module 3 literature anchors came back RED, and they stay red.**
They are on the front page rather than only in the design docs, because a
status section reporting only the green ones would misrepresent what this
pipeline found.

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

- **Module 2 integration** (run `30718392588`): MOFA2 trained on the three
  views; subtypes imbalanced (S1=20, S2=306, S3=76, S4=122); MOFA-vs-SNF
  concordance ARI **0.351** — moderate, not high. Raw output at
  `docs/results/module2-run-30718392588.txt`.
- **Module 3 credibility anchors** (run `30840373033`), three green, one red:

  | check | verdict | measured |
  |---|---|---|
  | `mutation_freq` | **PASS** | VHL 44.8%, PBRM1 30.5%, SETD2 10.1%, BAP1 8.6% — all inside published ranges, n=417 |
  | `ccab_signature` | **PASS** | ccA/ccB anti-correlation rho **−0.354**, p 6.5e-17, full 6+6 panels |
  | `bap1_survival` | **DIRECTIONALLY RIGHT, UNDERPOWERED** | HR **1.584**, 95% CI 0.967–2.595, p 0.068, n=417. This is a positive control, not a result: the direction matches the literature, the significance is out of this cohort's reach (Schoenfeld needs ~470 events) |
  | `methyl_strata` | **RED — a real negative result** | silhouette 0.1197, Kruskal p 1.3e-82, but cluster-vs-platform ARI **0.583** against a 0.25 ceiling |

- **Platform confound** (run `30911448546`): the cohort is **214 HM27 / 310
  HM450**, merged with no batch correction. The merged methylation partition
  tracks the **assay**, not the biology — and the within-platform 4-means
  silhouettes (HM27 0.0858, HM450 0.0489) are BOTH *below* the merged 0.1197,
  so the merge is what manufactures the apparent structure. Decision taken:
  keep all 524 cases, do **not** restrict to one platform, and apply **no**
  batch correction (too few cases are assayed on both platforms for ComBat to
  be validated, and the probe sets differ); adjust for platform as a covariate
  and take predictors only from the platform-clean factors instead.
  `SANITY_MAX_PLATFORM_ARI` stays at 0.25 and the m1–m4 anchor stays failing —
  its job now is to fail *informatively*. Raw output at
  `docs/results/platform-diagnosis-run-30911448546.txt`.
- **What the confound does NOT touch**: the MOFA subtypes are platform-clean
  (ARI **0.0058**), and the mutation-frequency and ccA/ccB anchors read no
  methylation matrix at all.

**Built but never run on the snapshot:** Module 4 — the Cox / penalised-Cox /
RSF survival arms (`R/functions_survival.R`), the from-scratch C-index and
calibration (`R/functions_model_eval.R`), the BAP1-from-expression classifier
(`python/bap1_classifier.py`) and the six DAG targets that wire them. The code
exists and its unit tests pass on synthetic fixtures; the targets have not been
built in any workflow, so **no C-index, calibration, hazard ratio,
discrimination or AUROC result is claimed anywhere.**

**Not yet done:** the dashboard (Module 5) and single-cell (Module 6).

**A limit on what the held-out figures will mean when they are produced.** The
MOFA factors are fitted on all 524 cases, and the 5000-gene variance filter is
computed over all samples, before either model draws its train/test split. Both
steps are outcome-blind — no label leak, and the BAP1 task stays non-circular —
but the test rows did help define the latent axes and the feature set, so the
held-out C-index and AUROC are *unsupervised-transductive*, not fully
out-of-sample. The reported optimism bounds the supervised component only.

## Reproducibility scope (read before trusting the CI badge)

- The research core is a **frozen 2016 snapshot**: `curatedTCGAData` **data
  version 2.0.1** (the `version=` argument, `CURATED_VERSION` in
  `R/constants.R`), snapshot `20160128`, hg19, served by the
  `curatedTCGAData` **package** 1.34.0 pinned in `renv.lock`. The two numbers
  are different things and both appear in the design docs; the one that
  determines the data is 2.0.1. It never updates.
- **CI does NOT run the full pipeline.** CI lints, runs unit tests on
  subsampled fixtures, and renders the dashboard from cached `_targets` /
  release-asset results. A green CI badge is **not** full reproduction.
- The full pipeline runs locally once (see `scripts/run_full_pipeline.R`);
  expected runtime and hardware are documented in `docs/runtime.md`.

## Skills mapping

_Populated in Module 5._

## Licence

MIT — see `LICENSE`.
