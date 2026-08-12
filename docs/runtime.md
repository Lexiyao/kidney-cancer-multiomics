# Runtime and hardware

The full pipeline is run **once** (`scripts/run_full_pipeline.R`) and the
resulting `_targets/` store is published as a GitHub **release asset**
(`scripts/freeze_release_assets.R`). CI and the dashboard restore that asset and
render from it. A green CI badge does **not** mean CI reproduced the analysis end
to end.

## Measured runtime — full pipeline, GitHub Actions run `31375702141`

These are **measured**, not estimated: they are the `seconds` field that
`targets` recorded for each target during run `31375702141` (2026-08-10,
`verify-module2.yml`, `HEAVY_PULL=true`, inside
`bioconductor/bioconductor_docker:RELEASE_3_23` on a GitHub-hosted
`ubuntu-latest` runner). The store's `meta` file was uploaded as a run artifact,
which is where the table comes from — **that artifact is not committed to this
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
| **Total target compute** | 36 timed targets | **596 s ≈ 9.9 min** |

The **job** took 15 min 28 s wall clock (09:40:25 → 09:55:53 UTC). The
difference — about 5.5 minutes — is container pull, package installation and the
test suite, none of which is pipeline compute.

Two things this table corrects about earlier planning estimates:

- **MOFA2 training is 82 % of the pipeline.** Everything else together is under
  two minutes.
- **The ingest download was not the bottleneck.** `mae_raw` took 44.9 s on that
  runner. Download time depends on network and on ExperimentHub cache state, so
  a cold first pull on a slow link will be slower than this — but the
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
container targets, not the record of a local run** — the full pipeline has not
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
| `ci.yml` | push to `main`, PR | Three jobs. `r-checks` and `py-checks` lint (R + Python) and run `testthat` / `pytest` on subsampled fixtures, inside the Bioconductor image, asserting the installed versions match `renv.lock`. `render-smoke` renders all five `.qmd` pages twice — once with **no** store, asserting every page degrades to stated gaps, and once against a **synthetic fixture store** (`scripts/fixture_store_pipeline.R`), asserting every page renders its figures with no gaps left. **It does not restore the `targets-store` release asset**, so it is not reproduction and the numbers it renders are invented stand-ins; the deploying render is `pages.yml`. |
| `pages.yml` | **`workflow_dispatch` only** | Restore the `targets-store` release asset, render `dashboard/`, deploy `_site/` to GitHub Pages. Deliberately manual: publishing makes the site public, and that decision is the repo owner's, not a side effect of a merge. The site is **not currently deployed**. |
| `cron.yml` | Weekly (Mon 04:17 UTC) + manual | Rebuild **only** `gdc_live_panel` (and assert that nothing else rebuilt), re-freeze the store, and separately rebuild the container as the dependency-drift check. **It does not deploy Pages** and it does not touch the frozen research core. |
| `verify-module2.yml` | manual | The only workflow that builds the research targets on real data; `HEAVY_PULL=true`. This is where every measured number in this repository comes from. |
| `heavy-pull.yml`, `diagnose-platform.yml`, `validate-container.yml` | manual | One-off verification runs. |

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
