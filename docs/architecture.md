# Architecture — module contracts

`targets` orchestrates seven stage-scoped modules. Pure functions live in
`R/functions_*.R` / `python/*.py`; all state is cached in `_targets/`.

| Module | Key targets | Owner files |
|---|---|---|
| 0 Scaffold | `scaffold_env_check` | `_targets.R`, `R/constants.R` |
| 1 Ingest/preprocess | `mae_raw`, `mae_qc`, `rna_mat`, `methyl_mat`, `cnv_mat`, `mut_annot` | `R/functions_ingest.R`, `R/functions_preprocess.R`, `R/functions_utils.R` |
| 2 Integrate | `mofa_model`, `mofa_factors`, `mofa_varexp`, `subtypes_mofa`, `snf_clusters`, `concordance` | `R/functions_integrate.R` |
| 3 Sanity | `sanity_results` | `R/functions_sanity.R` |
| 4 Model | `cox_fit`, `survival_metrics`, `bap1_auroc` | `R/functions_survival.R`, `R/functions_model_eval.R`, `python/bap1_classifier.py` |
| 5 Dashboard | `gdc_live_panel`, `dashboard_site` | `R/functions_gdc_live.R`, `dashboard/*.qmd` |
| 6 Single-cell (v1.1) | `purity_check`, `sc_object` | `R/functions_purity.R`, `python/singlecell_*.py` |

## Top-level globals (defined in `_targets.R`)

- `config` — parsed `config/params.yml`; Module 6 reads `config$singlecell`.
- `HEAVY_PULL` — bool; `config$heavy_pull` OR `HEAVY_PULL=true`. Module 1's
  `mae_raw` gates its `tar_cue` on it; every heavy pull checks it.

## Hard invariants (see design spec globalConstraints)

- Methylation = HM27+HM450 merged on common CpGs; never HM450 alone.
- RNA = `log2(x+1)` on RSEM-normalised values; never `DESeq2::vst`.
- Mutation is annotation on the n=417 subset, never a MOFA view.
- Survival model is low-dimensional (EPV cap = 10 events per variable; the
  predictor cap is derived at fit time as `floor(observed_OS_events / EPV_CAP)`,
  never hard-coded, and feature selection is never genome-wide); BAP1 classifier
  is non-circular (expression → BAP1 status).
- MOFA2 runs with `run_mofa(use_basilisk = FALSE)` against `RETICULATE_PYTHON`.
- CI never runs the full pipeline; heavy pulls are HEAVY_PULL-gated.

## Measured status (not aspiration)

Measured by GitHub Actions runs 30708943504 (2026-08-01), 30718392588 and
31375702141 (2026-08-10), all in `bioconductor/bioconductor_docker:RELEASE_3_23`
against the frozen curatedTCGAData 2.0.1 KIRC snapshot 20160128. Each bullet
names the run it comes from; nothing here is claimed without one.

- **Module 1 is fully materialised on the real snapshot** (`tar_make` with
  `HEAVY_PULL=true`, all targets built): `cohort_n` = 524; `rna_mat` 5000 × 524;
  `methyl_mat` 5000 × 524; `cnv_mat` 24776 × 524; `mut_annot` 417 × 7
  (`sample_id` character, plus one integer 0/1 column per driver gene). The
  canonical `mut_annot` contract therefore holds on real data, not only on
  synthetic fixtures.
- **Survival census measured** from `colData(mae)`: main cohort 173 OS events
  (33.1%, median follow-up 1188 d) → EPV-10 predictor cap 17; 5-year-restricted
  (1825 d) 148 events → cap 14. Phase 4 wires 5 predictors and does not spend
  the cap.
- **`renv.lock` is a complete, measured 218-package lock** — resolved from a
  real install, not partial or hand-authored.
- **Module 2 is built on real data** (run 30718392588): MOFA2 trained on the
  three views, 15 factors on n = 524; subtypes S1 20 / S2 306 / S3 76 / S4 122;
  MOFA-vs-SNF concordance ARI **0.351** (moderate, not high). Raw output at
  `docs/results/module2-run-30718392588.txt`.
- **Module 3 anchors have run** (run 31375702141): four of the five literature
  checks PASS; `methyl_strata` (m1–m4) is **RED** and stays red — it carries the
  suite's only 2 red expectations, in 2 anchor tests, out of 203. Pinned as an
  expected failure: a new red fails the job, and an m1–m4 anchor that starts
  passing fails it too.
- **Module 4 is fitted on real data** (run 31375702141): survival frame 519
  rows, 171 OS events (32.9 %), train 364 / test 155 with 124 training events →
  measured EPV-10 cap 12, 5 predictors used. Held-out C-index Cox **0.7486**
  (penalised 0.7492, RSF 0.7524), Cox optimism 0.0125; BAP1-from-expression
  AUROC **0.960** held-out on n = 413 with 36 mutants. In that fit only
  `stage_num` and `age_years` are significant, so the C-index is a
  stage-and-age model — see README.md before quoting it. Raw output at
  `docs/results/module4-run-31375702141.txt`.
- **Module 5 is built and rendering, and is NOT published.** Five Quarto pages
  render from the cached store; `.github/workflows/pages.yml` is
  `workflow_dispatch`-only on purpose, so no site is deployed.
- **Not yet done:** Module 6 (single-cell GSE159115 + the ESTIMATE purity gate)
  is **NOT BUILT** — no code for it exists in this repository, and no
  single-cell or purity claim is made anywhere.
