# kidney-cancer-multiomics

Reproducible somatic **multi-omics integration on TCGA-KIRC** (kidney renal
clear cell carcinoma), built as a reproducible analytical pipeline: ingest →
MOFA2 / SNF integration → literature positive controls → low-dimensional
survival + a non-circular BAP1-from-expression classifier → an interactive
Quarto/Plotly dashboard, orchestrated by `targets` and pinned by `renv` +
Docker.

**Dashboard status: built, not published, and the store it needs does not exist
yet.** The five-page site renders from a cached `_targets` store, and
`.github/workflows/pages.yml` is `workflow_dispatch`-only on purpose, so nothing
is deployed at `lexiyao.github.io/kidney-cancer-multiomics`. Publishing is
**not** a two-line change: `pages.yml` begins by running `gh release download
targets-store`, and this repository currently has **zero releases**, so a
dispatch would fail at its first step. The prerequisite chain is: run the full
pipeline (Docker path below) → `Rscript scripts/freeze_release_assets.R` →
upload the `targets-store` asset → then dispatch `pages.yml`. Until that is
done, render locally with `quarto render dashboard`, which produces the site
with a stated gap in place of every store-backed figure.

## Why this repo exists — skills → module map

Each row names the code **and the run that executed it**. Rows with nothing
behind them say so rather than being quietly dropped.

| Skill | Where it lives | Status — and the evidence |
|---|---|---|
| Somatic multi-omics integration | `R/functions_integrate.R` (MOFA2 + SNF) | **Run on real data** — 15 factors on n=524, four subtypes, MOFA-vs-SNF ARI 0.351 (run `30718392588`) |
| Transcriptomic / epigenomic / CNV curation | `R/functions_ingest.R`, `R/functions_preprocess.R` | **Run on real data** — `rna_mat` 5000×524, `methyl_mat` 5000×524, `cnv_mat` 24776×524, `mut_annot` 417×7 (run `30708943504`, transcript not committed) |
| Reproducible pipeline + version control + env manager | `_targets.R`, `renv.lock` (218 pkgs), `Dockerfile`, GitHub Actions | **Run** — every result below was produced inside `bioconductor/bioconductor_docker:RELEASE_3_23`, with the research packages installed by `BiocManager`; see "What the lockfile does and does not pin" below |
| R + Python in one pipeline | MOFA2 via `reticulate`/`basilisk`; `python/bap1_classifier.py` | **Run on real data** — MOFA2 trains against system Python, no conda pulled (run `30570220145`, transcript not committed) |
| ML + survival analysis | `R/functions_survival.R`, `R/functions_model_eval.R` (C-index + calibration from scratch) | **Fitted on real data** (run `31375702141`) — read the caveat under "Module 4" before quoting the C-index |
| Literature positive controls as tests | `R/functions_sanity.R`, `tests/testthat/test-sanity.R` | **Run on real data** — 24 anchor tests (201 expectations green, 2 red), 4 of 5 checks pass, the 4-cluster methylation check pinned as an expected red |
| Data viz (Plotly) + communicating to a non-technical reader | `dashboard/*.qmd` | **Built and rendering** — 5 pages, each degrading to a stated gap without the store |
| Public data → regularly updated tool | `R/functions_gdc_live.R`, `.github/workflows/cron.yml` | **BUILT, NOT YET RUN.** The weekly cron is wired to refresh the live GDC panel only, but no execution of it exists yet (the two scheduled runs on record are the earlier placeholder job), and it needs the `targets-store` release first |
| Single-cell | — | **NOT BUILT.** Module 6 (GSE159115 + ESTIMATE purity gate) is v1.1 and non-blocking; no code exists yet |
| Publication-quality static figures | — | **NOT BUILT.** Figures here are interactive HTML; there is no manuscript-figure export path |

## One-command run

The full pipeline is **intended** to be run once, with the `_targets/` store
then frozen as a GitHub release asset that Pages and the cron restore. **That
asset has not been produced yet** — `gh release list` on this repository is
empty — so the two commands below are the path to creating it, not a record of
it having been created:

```bash
docker build -t kidney-cancer-multiomics .
docker run --rm -v "$PWD":/work -w /work kidney-cancer-multiomics \
  Rscript scripts/run_full_pipeline.R
Rscript scripts/freeze_release_assets.R
```

Measured runtime: **596 s of target compute (≈ 10 min), 82 % of it MOFA2
training**, in a 15 min 28 s CI job. Per-stage table and reference hardware
(8-core / 16 GB / Docker) in [docs/runtime.md](docs/runtime.md), together with
what each workflow actually does.

## Cohort reality (read before quoting sample sizes)

The multi-omics n is **not ~530**. Main analysis cohort =
RNA + Methylation(HM27+HM450 merged on common CpGs) + CNV = **524**; the fully
intersected RNA+Methyl+CNV+Mutation set is **413**, and RNA+HM450+CNV+Mutation
is only **241**. Mutation is an **annotation** (BAP1/PBRM1 factor labels on the
n=417 subset), never a MOFA view. RNA is RSEM upper-quartile normalised Level-3
data — `log2(x+1)` + variable-gene filtering, **not** `vst`.

These counts are **measured** on the frozen 2016 `curatedTCGAData` snapshot
(data version 2.0.1, snapshot `20160128`, hg19) — primary tumours only — not
queried from the live GDC API. The two genuinely differ in both directions (the
2016 legacy MAF covers 417 cases, today's masked-somatic-mutation MAF fewer),
so every quoted n must say which source it came from. Only the dashboard's
[live GDC panel](dashboard/live-gdc.qmd) reports current-portal numbers, and
they are never mixed with snapshot numbers.

## What has actually run

**On the run ids below.** Four runs have their raw output committed under
`docs/results/`: `30718392588`, `30840373033`, `30911448546`, `31375702141`.
The Module 1 / container / census runs — `30570220145`, `30642823359`,
`30708943504` — do not; they are marked *(transcript not committed)* where they
appear, because a run id a reader cannot open is weaker evidence than one they
can.

Verified by real runs against the frozen snapshot inside
`bioconductor/bioconductor_docker:RELEASE_3_23`. **One of the five Module 3
literature checks came back RED — the merged-methylation 4-cluster check — and
it stays red** (it carries the suite's only 2 red expectations, out of 203, in 2 anchor
tests; both are that one finding). A status section reporting only the green
ones would misrepresent what this pipeline found.

- **Container chain** (run `30570220145`, transcript not committed): MOFA2
  trains via `run_mofa(use_basilisk = FALSE)` against the system Python — no
  conda pulled.
- **Cohort census** (run `30642823359`, re-verified with zero drift by
  `30708943504`; neither transcript committed): main cohort **524** cases;
  **413** with mutation; mutation MAF covers **417**. The **241** RNA + HM450 +
  CNV + Mutation figure quoted above comes from the same census.
- **Survival census** (run `30708943504`, transcript not committed): **173** OS
  events (33.1%, median follow-up 1188 d; **148** within 5 years), so the EPV-10
  predictor budget is 17 (14 restricted). This is the `colData` census over the
  522-case cohort, **not** the fitted frame — Module 4 below reports 171 events
  over the 519 rows that reach the fit, and the two are different denominators,
  not a contradiction. The pipeline recomputes this at fit time; Phase 4 wires 5.
- **Module 1 materialised** (run `30708943504`, transcript not committed):
  `rna_mat` 5000 × 524, `methyl_mat` 5000 × 524, `cnv_mat` 24776 × 524,
  `mut_annot` 417 × 7.
- **`renv.lock`**: a complete 218-package lock snapshotted from a machine with
  every Import installed — not hand-authored.
- **Module 2 integration** (run `30718392588`): MOFA2 trained on the three
  views; subtypes imbalanced (S1=20, S2=306, S3=76, S4=122); MOFA-vs-SNF
  concordance ARI **0.351** — moderate, not high. Raw output at
  `docs/results/module2-run-30718392588.txt`.
- **Module 3 credibility anchors** (run `31375702141` — the run that printed
  all five verdicts), four green, one red. The earlier run `30840373033`
  recorded only **four** checks (`mutation_freq`, `bap1_survival`,
  `methyl_strata`, `ccab_signature`) at **11 anchors / 57 passed / 4 failed**,
  and ended `conclusion: failure`; `subtype_platform` did not exist yet
  (`fn_check_subtype_platform` was added in response to run `30911448546`, which
  is later). Every measured value the two runs share is byte-identical; what
  changed between them is the **test definition**, not the data — the BAP1
  anchor's `p < 0.05` and `ci_low > 1` assertions were retired as
  mis-specified (`R/constants.R`), so its verdict moved from red to green while
  HR, CI and p stayed exactly the same. Run `31375702141` records 24 anchors /
  201 passed / 2 failed.

  | check | verdict | measured |
  |---|---|---|
  | `mutation_freq` | **PASS** | VHL 44.8%, PBRM1 30.5%, SETD2 10.1%, BAP1 8.6% — all inside published ranges, n=417 |
  | `ccab_signature` | **PASS** | ccA/ccB anti-correlation rho **−0.354**, p 6.5e-17, full 6+6 panels |
  | `subtype_platform` | **PASS** | MOFA subtypes vs assay platform ARI **0.0058** (first measured in run `30911448546`), p 0.53 |
  | `bap1_survival` | **DIRECTIONALLY RIGHT, UNDERPOWERED** | HR **1.584**, 95% CI 0.967–2.595, p 0.068, n=417. A positive control, not a result: the direction matches the literature, the significance is out of this cohort's reach — Schoenfeld needs ~470 events and the subset **records 147, 18 of them in the mutant arm** |
  | `methyl_strata` (4 clusters in merged HM27+HM450) | **RED — a real negative result** | silhouette 0.1197, Kruskal p 1.3e-82, but cluster-vs-platform ARI **0.583** against a 0.25 ceiling |

- **Platform confound** (run `30911448546`): the cohort is **214 HM27 / 310
  HM450**, merged with no batch correction. The merged methylation partition
  tracks the **assay**, not the biology — and the within-platform 4-means
  silhouettes (HM27 0.0858, HM450 0.0489) are BOTH *below* the merged 0.1197,
  so the merge is what manufactures the apparent structure. Decision taken:
  keep all 524 cases, do **not** restrict to one platform, and apply **no**
  batch correction (only **3** cases are assayed on both platforms and the probe
  sets differ, so ComBat could be neither validated nor trusted); adjust for
  platform as a covariate and take predictors only from the platform-clean
  factors instead. `SANITY_MAX_PLATFORM_ARI` stays at 0.25 and the methylation
  cluster anchor stays failing — its job now is to fail *informatively*. Raw output at
  `docs/results/platform-diagnosis-run-30911448546.txt`.
- **What the confound does NOT touch**: the MOFA subtypes are platform-clean
  (ARI **0.0058**), and the mutation-frequency and ccA/ccB anchors read no
  methylation matrix at all.
- **Module 4, fitted on real data** (run `31375702141`, raw output at
  `docs/results/module4-run-31375702141.txt`): survival frame 519 rows, **171 OS
  events (32.9%)** — the FITTED frame, against the 173-event `colData` census
  above — predictors `Factor1, Factor4, age_years, stage_num,
  platform`; Cox train 364 / test 155 with **124 training events**, so the
  measured EPV-10 cap is 12 predictors and 5 are used. Held-out C-index **Cox
  0.7486 / penalised Cox 0.7492 / RSF 0.7524**, Cox optimism (apparent −
  held-out) **0.0125**. BAP1-from-expression AUROC **0.960 held-out / 0.958
  cross-validated** on n=413 with 36 mutants.

  **What that C-index is not.** In the fitted Cox model only the clinical terms
  reach significance — `stage_num` p = 9.6e-17, `age_years` p = 9.7e-06 — while
  **neither MOFA factor does**: `Factor1` HR 0.982 (p = 0.60), `Factor4` HR
  0.966 (p = 0.24). A C-index of 0.75 is a **stage-and-age** model. On this
  cohort the multi-omics factors add nothing detectable to it, and that is the
  result, not a presentation problem. The three arms landing within 0.004 of
  each other says the same thing from another direction.

  **The BAP1 AUROC rests on 36 mutants** out of 413, so the held-out split holds
  only a handful of positives and the interval around 0.96 is wide. The task is
  genuinely non-circular — the label is an external mutation call — but the
  point estimate is not a validated classifier.

**A limit on every held-out figure above.** The MOFA factors are fitted on all
524 cases and the 5000-gene variance filter is computed over all samples, before
either model draws its train/test split. Both steps are outcome-blind — no label
leak, and the BAP1 task stays non-circular — but the test rows helped define the
latent axes and the feature set, so the held-out C-index and AUROC are
*unsupervised-transductive*, not fully out-of-sample. The reported optimism
bounds the supervised component only.

**Not built:** single cell (Module 6, GSE159115) and its ESTIMATE purity gate.
No single-cell or purity claim is made anywhere in this repository.

## What CI does (and does not) do

- **`ci.yml`** lints and runs `testthat` / `pytest` on subsampled fixtures inside
  the Bioconductor image, asserting installed versions match `renv.lock`, and a
  third job renders all five dashboard pages twice — once with no store
  (asserting every page degrades to stated gaps) and once against a **synthetic
  fixture store** (asserting every page renders with no gaps left). It does
  **not** restore the `targets-store` release asset, so the figures it renders
  are invented stand-ins and a green badge is not **end-to-end reproduction**.
- **`pages.yml`** *would* restore the `targets-store` release asset, render
  `dashboard/`, and deploy to Pages. **Manual trigger only, and never yet
  dispatched**; with no release present it would currently fail on its first
  step.
- **`cron.yml`** is **wired but has never run in its current form.** As
  committed it refreshes **only** the live GDC statistics panel, asserts that no
  frozen target rebuilt, re-freezes the store, and in a separate job rebuilds
  the container as the dependency-drift check. It does **not** deploy Pages and
  the 2016/hg19 research core never updates. The two scheduled executions that
  exist — runs `30795010494` (2026-08-03) and `31359626479` (2026-08-10) — ran
  the earlier `drift-check` scaffold, whose only step was a placeholder; the
  live-panel job landed on `main` afterwards. The next scheduled run will also
  fail until the `targets-store` release exists, because it restores the store
  before doing anything else.
- **Three workflows build research targets on real data** (`HEAVY_PULL=true`),
  all manual, and between them they produced every measured number here:
  **`verify-module2.yml`** (Modules 1–4 and the anchor suite — runs
  `30718392588`, `30840373033`, `31375702141`), **`heavy-pull.yml`** (the
  cohort and survival census, the Module 1 matrix dimensions and the
  218-package lock — runs `30642823359`, `30708943504`), and
  **`diagnose-platform.yml`** (the platform confound: the 214/310 split, every
  factor AUC, the within-platform silhouettes, `methyl_platform_overlap` — run
  `30911448546`). Push/PR CI never builds them.

See [docs/runtime.md](docs/runtime.md) for the full workflow table.

## Reproducibility scope (read before trusting the CI badge)

- The research core is a **frozen 2016 snapshot**: `curatedTCGAData` **data
  version 2.0.1** (the `version=` argument, `CURATED_VERSION` in
  `R/constants.R`), snapshot `20160128`, hg19, served by the
  `curatedTCGAData` **package** 1.34.0 pinned in `renv.lock`. The two numbers
  are different things and both appear in the design docs; the one that
  determines the data is 2.0.1. It never updates.
- The full pipeline runs once (see `scripts/run_full_pipeline.R`); measured
  runtime and reference hardware are in [docs/runtime.md](docs/runtime.md).
- **What the lockfile does and does not pin.** `renv.lock` is a complete
  218-package lock, but it is not what produced the numbers above, and saying
  otherwise would overclaim. The three real-data workflows install their
  research packages with `BiocManager::install()` from an unpinned name list
  inside the pinned image — same container, not the same resolved versions.
  `renv.lock` is *generated* by `heavy-pull.yml`, *restored* by the
  `Dockerfile`'s `renv::restore()`, and *asserted against* by `ci.yml` (which
  produces no research numbers). So: right image, `BiocManager`-resolved
  packages, with the lockfile pinning the environment the container builds
  rather than the environment the figures came out of. Closing that gap means
  adding a `renv::restore()` and a version assertion to `verify-module2.yml`
  and re-running.

## Known limitations (stated openly)

- **n is 241–524**, not ~530, set by modality intersection. The survival model
  runs on the main cohort (519 rows reach the fit) and is kept low-dimensional
  under an EPV-10 budget measured on the *training* events (124 → cap 12, 5
  used). No genome-wide feature selection, ever.
- **The methylation platform confound is the largest single limitation, and it
  is measured, not suspected.** `methyl_mat` is `cbind(HM27, HM450)` on common
  CpGs with no batch correction. Factor2 (AUC 0.888) and Factor5 (0.818) are
  substantially assay effects; only Factor1 (0.500) and Factor4 (0.535) among
  the wired predictors carry no detectable platform signal, which is why they —
  and not the higher-variance factors — are the ones in the model. The choice
  was made **outcome-blind**, on platform cleanliness and variance explained.
- **The merged-methylation 4-cluster positive control FAILS** and is kept
  failing, pinned as an expected red in the container run (`verify-module2.yml`,
  manual dispatch): a new red fails the job, and this anchor starting to *pass*
  also fails the job, because that would mean the recorded negative no longer
  holds. **What it is not:** this is not a failure to reproduce published
  methylation biology. TCGA KIRC (Nature 2013) defines m1–m4 as mRNA
  **expression** subtypes and reports no four DNA-methylation strata, so `k = 4`
  here is a design choice carried over from those expression subtypes. The red
  is a finding about *this* merged HM27/HM450 matrix — it partitions by assay —
  not a contradiction of the paper.
- **The survival C-index is carried by stage and age**, not by the multi-omics
  factors (see Module 4 above). Reporting 0.75 without that sentence would be
  the central overclaim this project is designed not to make.
- **BAP1 survival is underpowered by construction**: ~470 events needed for 80 %
  power, 147 recorded (18 in the mutant arm). Directionally consistent with the
  literature; neither a confirmation nor a negative finding.
- **Held-out is not out-of-sample.** Unsupervised-transductive, as above, and
  there is no external RCC validation cohort.
- **Frozen 2016 / hg19 snapshot.** Genuine live update is limited to the GDC
  statistics panel, and the published page advances only when Pages is
  dispatched.
- **Bulk → single-cell mapping is not attempted**, so its purity/immune confound
  is currently an unaddressed design risk rather than a checked one. The
  ESTIMATE gate is specified (design spec §10) and not implemented.
- **CI does not run the full pipeline**; it lints and tests on fixtures, and the
  site renders from cached results.

## Licence

MIT — see `LICENSE`.
