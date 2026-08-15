# Runtime and hardware

The full pipeline is **designed** to be run once (`scripts/run_full_pipeline.R`)
with the resulting `_targets/` store published as a GitHub **release asset**
(`scripts/freeze_release_assets.R`), so that Pages and the cron can restore that
asset and render from it. **The asset does not exist yet** - this repository has
no releases - so `pages.yml` and `cron.yml`, both of which begin with `gh
release download targets-store` under `set -euo pipefail`, would currently fail
at their first step. A green CI badge does **not** mean CI reproduced the
analysis end to end.

## Measured runtime - full pipeline, GitHub Actions run `31375702141`

These are **measured**, not estimated: they are the `seconds` field that
`targets` recorded for each target during run `31375702141` (2026-08-10,
`verify-module2.yml`, `HEAVY_PULL=true`, inside
`bioconductor/bioconductor_docker:RELEASE_3_23` on a GitHub-hosted
`ubuntu-latest` runner). The store's `meta` file was uploaded as a run artifact,
which is where the table comes from - **that artifact is not committed to this
repository**, so unlike the Module 2/3/4 figures under `docs/results/`, this
table cites a run id a reader cannot open. It is recorded here rather than
dropped, and labelled rather than left to look like the others.

| Stage | Targets | Measured wall time |
|---|---|---|
| Ingest (`curatedTCGAData`, snapshot 20160128) | `mae_raw`, `mae_qc`, `clinical`, `rna_full`, `cohort_n` | 49.7 s |
| Preprocess (RNA / HM27+HM450 merge / CNV / mutation) | 9 targets | 31.1 s |
| **MOFA2 training** | `mofa_model` alone is **486.1 s**; with factors, varexp, subtypes, mutation annotation and the platform diagnostic | **490.8 s** |
| SNF sensitivity + concordance | `snf_clusters`, `concordance` | 6.4 s |
| Literature positive controls | `sanity_results` | 9.1 s |
| Module 4 (Cox / penalised Cox / RSF / C-index / BAP1 classifier) | 7 targets | 8.8 s |
| **Total of the rows above** | 30 targets | **595.9 s ≈ 9.9 min** |

**A discrepancy this table used to hide.** The total row previously read
"36 timed targets | 596 s". The six stage rows enumerate **30** targets, and
their measured times sum to 595.9 s - so for both figures to hold, six further
timed targets would have to account for ~0.1 s between them. That is not
credible: run `31375702141` ran a bare `targets::tar_make()`
(`.github/workflows/verify-module2.yml`), which also builds `gdc_live_panel` -
a live HTTPS query to `api.gdc.cancer.gov` - plus `subtypes_df`,
`km_subtype_df`, `sanity_table` and `dashboard_site`. `_targets.R` declares 41
targets in total. Either the count or the partition was wrong, and the `meta`
artifact that would settle it is not committed, so the **target count is
withdrawn** and only the measured times are kept. The stage times themselves
are unaffected - they are the per-target `seconds` fields as recorded.

The **job** took 15 min 28 s wall clock (09:40:25 → 09:55:53 UTC). The
difference - about 5.5 minutes - is container pull, package installation and the
test suite, none of which is pipeline compute.

Two things this table corrects about earlier planning estimates:

- **MOFA2 training is 82 % of the pipeline.** Everything else together is under
  two minutes.
- **The ingest download was not the bottleneck.** `mae_raw` took 44.9 s on that
  runner. Download time depends on network and on ExperimentHub cache state, so
  a cold first pull on a slow link will be slower than this - but the
  "10–25 minutes dominated by the HM450 HDF5 download" figure that appears in
  the plan was never measured and is not what happened.

Nothing about the survival or classifier stage is expensive: the model is
low-dimensional by design (5 predictors), so Module 4 costs under 9 seconds.

## Site render

Measured locally (macOS, R 4.6.0, Quarto CLI 1.8.27), rendering all five pages
of `dashboard/` **without** the `_targets` store restored: **10.4 s**, every page
produced, each showing a stated `PENDING` block where a target was unavailable.
With the store restored the render is slower, because the Plotly figures are
then actually built; that figure has not been measured here and is not guessed
at.

## Reference hardware

8-core CPU, 16 GB RAM, Docker. No GPU required. **This is the specification the
container targets, not the record of a local run** - the full pipeline has not
been executed on a local workstation in this project's history; every measured
figure above comes from CI containers.

## One-command run

```bash
docker build -t kidney-cancer-multiomics .
docker run --rm -v "$PWD":/work -w /work kidney-cancer-multiomics \
  Rscript scripts/run_full_pipeline.R
```

Then freeze the store as a release asset for CI/Pages:

```bash
Rscript scripts/freeze_release_assets.R
```

## CI scope (honest)

This section describes the workflow files **as committed**, which is not what
the plan text for Phase 5 describes. Where they differ, the files are what runs.

| Workflow | Trigger | What it actually does |
|---|---|---|
| `ci.yml` | push to `main`, PR | Three jobs. `r-checks` and `py-checks` lint (R + Python) and run `testthat` / `pytest` on subsampled fixtures, inside the Bioconductor image, asserting the installed versions match `renv.lock`. `render-smoke` renders all five `.qmd` pages twice - once with **no** store, asserting every page degrades to stated gaps, and once against a **synthetic fixture store** (`scripts/fixture_store_pipeline.R`), asserting every page renders its figures with no gaps left. **It does not restore the `targets-store` release asset**, so it is not reproduction and the numbers it renders are invented stand-ins; the deploying render is `pages.yml`. |
| `pages.yml` | **`workflow_dispatch` only** | Restore the `targets-store` release asset, render `dashboard/`, deploy `_site/` to GitHub Pages. Deliberately manual: publishing makes the site public, and that decision is the repo owner's, not a side effect of a merge. The site is **not currently deployed**. |
| `cron.yml` | Weekly (Mon 04:17 UTC) + manual | Rebuild **only** `gdc_live_panel` (and assert that nothing else rebuilt), re-freeze the store, and separately rebuild the container as the dependency-drift check. **It does not deploy Pages** and it does not touch the frozen research core. |
| `verify-module2.yml` | manual | Builds the research targets on real data (`HEAVY_PULL=true`): Modules 1-4 and the credibility-anchor suite. Runs `30718392588`, `30840373033`, `31375702141`. |
| `heavy-pull.yml` | manual | Also builds research targets on real data (`HEAVY_PULL=true`, bare `tar_make()`). Produced the 524/413/417/241 cohort census, the 173-event survival census, the Module 1 matrix dimensions, and the generated 218-package `renv.lock`. Runs `30642823359`, `30708943504`. |
| `diagnose-platform.yml` | manual | Also builds research targets on real data (`HEAVY_PULL=true`, bare `tar_make()`). Produced the whole platform-confound block: the 214/310 split, every factor-vs-platform AUC, the within-platform silhouettes, and `methyl_platform_overlap = 3`. Run `30911448546`. |
| `validate-container.yml` | manual | Container/toolchain validation only; builds no research targets. |

- `render-smoke` is a separate job on purpose: `tests/testthat.R` refuses to
  skip the store-backed literature anchors against a populated store, so a
  fixture store in the same workspace would make the real-data anchors error
  against synthetic values.
- CI **never** runs the full pipeline on push or PR. The HM450 pull, MOFA2
  training and the survival fits are guarded behind `HEAVY_PULL` and run only in
  the manual container workflows above.
- Consequence for the "regularly updated tool" claim, stated plainly: the cron
  keeps the **cached** live-GDC panel current, but the **published** page
  advances only when `pages.yml` is dispatched. See `dashboard/live-gdc.qmd`,
  which says the same thing to the reader.
