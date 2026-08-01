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

Measured by GitHub Actions run 30708943504, 2026-08-01, in
`bioconductor/bioconductor_docker:RELEASE_3_23` against the frozen
curatedTCGAData 2.0.1 KIRC snapshot 20160128.

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
- **Not yet done:** the survival model has NOT been fitted, and Module 2
  (MOFA2 integration) has NOT been run. Modules 2–6 remain unmaterialised.
