# Design Spec — Reproducible Multi-Omics Pipeline for Kidney Cancer (TCGA-KIRC)

- **Date:** 2026-07-30
- **Status:** Design approved. Corrected against an independent feasibility review (2026-07-30).
- **Author:** Zixi (Lexi) Yao — github.com/Lexiyao
- **Repo:** `kidney-cancer-multiomics` → rendered site at `lexiyao.github.io/kidney-cancer-multiomics`

## Purpose

A portfolio project demonstrating **reproducible somatic multi-omics integration** on public
kidney-cancer data, built as a proper reproducible analytical pipeline (RAP) with an
interactive dashboard. It exists to fill specific, verifiable gaps between the author's
current CV (germline / statistical-genetics genomics: genetic association, 100k Genomes,
PGS, reproducible R pipelines) and the somatic / functional multi-omics skills required by
target roles:

- **Cambridge Early Cancer Institute — Mitchell Group RA (RE50394):** kidney-cancer somatic
  multi-omics integration (WGS / transcriptomic / epigenomic / single-cell), reproducible
  pipelines, clear visualisation.
- **Imperial IGHI — data-analyst-type roles:** public data → regularly updated tool, R + Python,
  Plotly/Tableau-class visualisation, RAP + version control + environment managers,
  communicating to non-technical audiences.
- **Kar Lab (cancer genetic epidemiology):** stays on the author's cancer-genomics main line.

## Design principle: honesty over impressiveness

This project is a credibility artifact. Every claim it makes must survive an interviewer's
follow-up question. The feasibility review found several places where the original design
would not survive that test; those are corrected here, and the residual limitations are
stated openly (see §12) because stating them **is** the maturity signal.

---

## 1. Skills → module traceability

| Required skill | Where it lives |
|---|---|
| Somatic multi-omics integration (Mitchell core) | Modules 1–3 |
| Transcriptomic / epigenomic / CNV curation | Modules 1–2 |
| Single-cell (Mitchell) | Module 7 (v1.1, non-blocking) |
| Reproducible pipeline + version control + env manager (both roles, explicit) | Scaffold |
| R + Python (both roles) | Module 3 (MOFA2 is genuinely R+Python), Module 5 (sklearn), Module 7 (scanpy) |
| ML + survival / epi (both + author strength) | Module 5 |
| Data viz (Plotly) + communicate to non-technical (Imperial) | Module 6 |
| Public data → regularly updated tool (Imperial) | Module 6 live-GDC panel + CI cron (§9) |
| Publication-quality visualisation (Mitchell) | Module 6 |

---

## 2. Data — corrected sample sizes and access

Cohort: **TCGA-KIRC** (kidney renal clear cell carcinoma, ccRCC). Access via
`curatedTCGAData` 1.34.0 (Bioconductor 3.23) — versioned `MultiAssayExperiment`, snapshot
`20160128`, hg19 legacy. This is the reproducible choice and is retained.

**The cohort is NOT ~530 for multi-omics.** Complete-case intersection (case/patient level,
GDC open-access):

| Modality combination | cases |
|---|---|
| RNA-seq | 533 |
| Methylation (any platform) | 535 — HM450: 319, HM27: 219, overlap: 3 |
| CNV (gene-level) | 532 |
| Masked somatic mutation | 374 |
| RNA + Methylation(any) + CNV | ~528 |
| RNA + Methylation(any) + CNV + Mutation | 370 |
| RNA + HM450 + CNV + Mutation | 267 |

**Two bottlenecks:** (1) open-access masked MAF covers only 374 cases; (2) KIRC methylation is
split across two nearly non-overlapping platforms (HM450 319 / HM27 219). Hard-coding "450k"
would silently halve the cohort.

**Decision — main analysis cohort:** `RNA + Methylation(HM27+HM450 merged on common CpGs) + CNV
≈ 528`. Mutation is **not** a modeling view; it is used as annotation on its `n=374` subset,
reported with stratification.

**Model-complexity implication (load-bearing):** at `n≈270` the ccRCC 5-year OS event count is
~90–110. At events-per-variable ≈ 10 that caps the effective survival model at ~10 predictors.
Module 5 is therefore a **low-dimensional model on factors/subtypes + a few clinical variables**,
never a genome-wide feature selection.

### 2a. RNA-seq object — do not call it `vst`

`curatedTCGAData`'s `KIRC_RNASeq2GeneNorm` is TCGA Level-3 **RSEM upper-quartile normalised
values**, not raw integer counts. `DESeq2::vst` assumes raw counts and is methodologically
wrong here. **Decision:** keep the normalised object, apply `log2(x + 1)` + variable-gene
filtering, and label it exactly that ("log-transformed normalised expression"). MOFA2 does not
need a count scale, so no data-source change is required. (Alternative, not chosen: pull GDC
harmonized STAR counts via `TCGAbiolinks` 2.40.0 for true `vst`, at the cost of losing the
versioned ExperimentHub snapshot.)

---

## 3. Architecture — modules (orchestrated by `targets`)

Each module has one purpose and explicit inputs/outputs; each is independently testable.

1. **ingest** (R): `curatedTCGAData` → `MultiAssayExperiment`; harmonise sample IDs; QC;
   cache versioned RDS. → `MultiAssayExperiment`, MultiAssayExperiment 1.38.0.
2. **preprocess** (R): per-omics transforms aligned on the ~528 common samples:
   - RNA: `log2(x+1)` on RSEM-normalised, top-variable genes (Gaussian).
   - Methylation: β → M-values; drop SNP-adjacent and sex-chromosome probes; merge HM27/HM450
     on common CpGs; top-variable CpGs (Gaussian).
   - CNV: GISTIC gene-level thresholded, continuous/ordinal (Gaussian).
   - Mutation: **excluded as a view** — see §6a.
3. **integrate** (R, **genuinely R+Python**): `MOFA2` 1.22.0 latent factors as the main method;
   `SNFtool` 2.3.1 as a cheap sensitivity analysis. Deliverable: factor matrix,
   variance-explained per omics per factor, subtype assignments, and a **two-method concordance
   check** (MOFA-derived clusters vs SNF clusters). **`iClusterPlus` is cut** (hours-to-days on
   4 omics × 300+ samples, redundant with MOFA2).
4. **sanity** (R, credibility anchor — built early, see §7): literature positive controls as
   real `testthat` assertions.
5. **model** (R + Python):
   - Survival (R): low-dimensional Cox / penalised Cox (`glmnet` 5.0) / random survival forest
     (`randomForestSRC` 3.6.2) on factors/subtypes + clinical; **held-out split or nested CV
     with optimism reported** (§6c); C-index + calibration **reusing the author's
     `model-evaluation-from-scratch` code**.
   - Classifier (Python `scikit-learn`): a **non-circular** supervised task — predict **BAP1
     mutation status** from expression (literature-anchored ccRCC signal; labels from the `n=374`
     mutation subset, features from expression). CV AUROC + held-out. This is the R+Python proof
     and avoids the tautology in §6b.
6. **report + dashboard** (Quarto + Plotly): Quarto site (publication-quality figures) plus an
   interactive dashboard (subtypes, survival curves, factor loadings, gene views) + a **live GDC
   panel** (§9). Deployed to GitHub Pages.
7. **singlecell** (Python `scanpy`, **v1.1, non-blocking**): GSE159115 (10x H5); QC, clustering,
   annotation; map bulk subtype signatures onto single-cell programs — **after** a purity/immune
   confound check (§10).

**Reproducibility scaffold:** `renv` (R lock) + `targets` (DAG) + `Dockerfile` +
GitHub Actions CI. See §5, §8.

---

## 4. Data flow

```
curatedTCGAData(20160128)
  → ingest (MultiAssayExperiment)
  → preprocess (aligned RNA / Methyl / CNV matrices, n≈528)
  → integrate (MOFA2 factors + subtypes; SNF sensitivity)
  → { sanity (positive-control assertions) ,
      model (survival held-out + BAP1 classifier) }
  → report/dashboard (+ live GDC panel)
  → GitHub Pages
singlecell (GSE159115) — separate branch, joins dashboard at v1.1
```

`targets` caches each stage; only changed stages re-run.

---

## 5. Reproducibility scaffold — the real build blocker is MOFA2

MOFA2 is **not pure R**: `SystemRequirements: Python (>=3), numpy, pandas, h5py, scipy, sklearn,
mofapy2`, called through `reticulate` + `basilisk`. `basilisk` builds its own conda env, which
collides with `renv` (R-only) + Docker layering + a reticulate-pointed Python. **This is the most
likely place to lose a day, so it is solved in the scaffold phase, not deferred to Module 3.**

Resolution: the Dockerfile installs Python + `mofapy2` explicitly, sets `RETICULATE_PYTHON`, and
forces basilisk to an external/system env (`basilisk.useSystemDir` / `BASILISK_EXTERNAL_DIR`) so
it never downloads conda inside the container.

Reframing bonus: because MOFA2 is itself R+Python, the "bilingual environment" is **genuinely
necessary**, not a contrived demo — a stronger story than an artificial R/Python split.

---

## 6. Methodology guards

**(a) Mutation is annotation, not a MOFA view.** At `n≈370`, only ~a dozen KIRC genes exceed 5%
mutation frequency (VHL, PBRM1, SETD2, BAP1, MTOR, KDM5C…). A sparse 10–15-column Bernoulli view
contributes ~nothing and MOFA2's Bernoulli likelihood converges worse than Gaussian. Instead,
mutation status is an **external label for factor interpretation** ("which factor tracks
BAP1/PBRM1 status?") — dodges the numerical problem and yields tellable biology. CNV as
GISTIC-thresholded gene-level (continuous/ordinal) remains a valid view.

**(b) No tautological classifier.** Deriving subtype labels from omics and then predicting them
from the same omics is circular; high accuracy is guaranteed and proves nothing. Module 5's
classifier is a **non-circular** target — BAP1 status from expression (§3.5). (Optional future
variant: minimal transferable signature distillation validated on an **external** RCC cohort.)

**(c) Subtype → survival needs held-out.** Computing log-rank p on the same samples that defined
the subtypes is a textbook optimism trap. Minimum bar: split subtype discovery from survival
testing, or nested CV with optimism reported. Better (future): external RCC cohort validation.

---

## 7. Testing — positive controls as the credibility anchor

Built **early** (Module 4), not at the end. Real `testthat` / `pytest` assertions, not just
figures:

- VHL / PBRM1 / SETD2 / BAP1 mutation frequencies fall within published ccRCC ranges.
- BAP1-mutant tumours show worse OS.
- Recovery of TCGA KIRC methylation strata (m1–m4).
- ccA / ccB expression-signature separation.

Unit tests cover deterministic helpers (ID harmonisation, matrix alignment, transforms) on
subsampled fixtures. Target ~80% coverage on the helper/utility layer; the analysis itself is
validated by these literature positive controls. **This suite is the single most persuasive part
of the repo.**

---

## 8. CI scope — stated honestly

A GitHub Actions standard runner (2 core / 7 GB / 6 h) **cannot** run the full pipeline: HM450 is
a large HDF5 download (ExperimentHub cache unreliable across runs), and MOFA2 + survival + scanpy
need real compute.

- **CI job =** lint + unit tests on subsampled fixtures + **render the dashboard from cached
  `_targets` / release-asset results**. Heavy re-pull is guarded behind a flag.
- **Full pipeline =** run locally once; outputs frozen as release assets / committed `_targets`
  cache. README states expected **runtime and hardware**.
- The README must not imply that a green CI badge equals full reproduction.

---

## 9. The "regularly updated tool" claim — honest resolution

The multi-omics core is a frozen 2016 snapshot; a weekly cron re-running it produces
byte-identical output. Claiming "public data → regularly updated tool" on frozen data would break
under one interview question. Resolution:

- Keep the research core on the frozen, versioned snapshot (reproducibility win).
- Add **one dashboard panel driven by a live GDC API query** (current TCGA-KIRC sample counts /
  clinical distribution), refreshed by the weekly cron. This genuinely updates and honestly
  demonstrates "pull public data on a schedule → update a tool."
- The cron's stated job is: (i) refresh the live GDC panel, (ii) dependency/environment drift
  detection + container rebuild. It does **not** claim the frozen research data updates.

---

## 10. Single-cell confound check (Module 7)

Feasible locally (16 GB RAM, ~50k cells in scanpy). Dataset: **GSE159115** (10x H5, ccRCC +
microenvironment; direct H5 read).

Scientific guard: if bulk subtypes are driven mainly by **tumour purity / immune infiltration**
(very common in ccRCC), mapping them onto single cells merely rediscovers "immune-cell
proportion," not a cellular basis for the subtypes. So **before** mapping, run an
ESTIMATE/purity check to test whether the subtypes are a purity proxy; if they are, re-frame the
module's conclusion. That check is itself tellable content.

---

## 11. Scope & MVP phasing

Even after cutting iClusterPlus and the mutation view, this is a heavy project. Delivery is
ordered so the first six modules form a complete, releasable, tellable project; single-cell is an
independent increment that cannot block release.

1. **Scaffold** — `renv` + `targets` + Docker, **including the MOFA2 basilisk/reticulate/Docker
   resolution (§5)**, CI hello-world, Quarto → Pages deploy.
2. **Ingest + preprocess** — RNA + methylation + CNV, `n≈528`, correct terminology (no `vst` on
   RSEM).
3. **Integrate** — MOFA2 main + SNF sensitivity; mutation as annotation.
4. **Sanity-check suite** — known ccRCC positive controls as `testthat` assertions (credibility
   anchor, deliberately early).
5. **Model** — survival (low-dim, held-out) + BAP1-status classifier (R+Python); reuse
   `model-evaluation-from-scratch`.
6. **Dashboard + README** — skills-mapping table, live GDC panel, honest CI/runtime docs.
7. **Single-cell** (GSE159115) — **v1.1, non-blocking.**

---

## 12. Known limitations (state openly in README + interview)

- **n is 267–370**, not ~530, set by modality intersection; the model is kept low-dimensional
  accordingly.
- **Frozen 2016 / hg19 snapshot**; genuine live-update is limited to the GDC statistics panel.
- **Subtype → survival optimism** is controlled by held-out / nested CV, not (yet) by an external
  cohort.
- **Bulk → single-cell mapping** is confounded by purity / immune infiltration; this is checked,
  not assumed.
- **CI does not run the full pipeline**; it tests and renders from cached results.

---

## 13. Deliverable shape

- New public GitHub repo `kidney-cancer-multiomics` → GitHub Pages site.
- README: motivation, **explicit skills-mapping table (§1)**, one-command run (Docker/`targets`)
  with runtime + hardware, results summary, figures, and the §12 limitations.
- A short interview talking-points doc: end-to-end narrative plus the honest caveats (which
  demonstrate methodological maturity).

## Decisions locked

1. Dashboard = Quarto Dashboard + Plotly (GitHub Pages, no server). **Confirmed.**
2. Repo = `kidney-cancer-multiomics` at `~/Desktop/CV/kidney-cancer-multiomics`.
3. Single-cell = **included, last, independently releasable** (delivery-risk isolation).
