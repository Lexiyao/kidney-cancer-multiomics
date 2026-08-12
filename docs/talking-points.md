# Interview talking points

Every number below is tied to a named GitHub Actions run, and the run id is
quoted beside each figure. `docs/results/*.txt` carries the **raw output for the
Module 2 / 3 / 4 figures** — runs `30718392588`, `30840373033`, `30911448546`,
`31375702141`. The Module 1 cohort and census figures, and the runtime table in
`docs/runtime.md`, cite runs whose transcripts are **not committed**; those are
marked "(run id only; transcript not committed)" where they appear, because a
run id a reader cannot open is weaker evidence than one they can, and saying so
is cheaper than being caught not saying it. If a question goes somewhere no run
has been, the answer is "not measured" — that is the whole posture of this
project.

## The 60-second story

TCGA-KIRC somatic multi-omics, built as a reproducible analytical pipeline.
`curatedTCGAData` gives a versioned `MultiAssayExperiment` (frozen 2016
snapshot, data version 2.0.1, hg19). I align RNA (log-normalised, **not** vst —
it is RSEM Level-3 data), methylation (HM27 + HM450 merged on common CpGs —
hard-coding "450k" would silently halve the cohort from 524 to 241 — run
`30708943504`, run id only; transcript not committed), and CNV on the 524
common cases. MOFA2 — genuinely R + Python, via reticulate against the
container's system Python so no conda is pulled at runtime — gives 15 shared
latent factors; SNF is a cheap sensitivity check and I report the two-method
concordance (ARI 0.351, moderate, not high). Mutation is annotation, not a view:
I ask which factor tracks BAP1/PBRM1 status. Survival is low-dimensional under
an EPV cap measured on the training events, on a held-out split with optimism
reported, and the Python classifier is non-circular — it predicts BAP1 status
from expression. The whole thing is pinned by `renv` + Docker and shipped as a
Quarto/Plotly dashboard.

## The parts I lead with

- **The positive-control test suite.** VHL/PBRM1/SETD2/BAP1 mutation frequencies
  inside published ccRCC ranges (44.8 / 30.5 / 10.1 / 8.6 %), BAP1-mutant worse
  OS, ccA/ccB separation (rho −0.354), MOFA subtypes independent of assay
  platform (ARI 0.0058), four stable clusters in the merged methylation matrix
  — as **24 real `testthat`
  anchor tests** (201 expectations green, 2 red) that execute against the frozen
  store in the **`verify-module2.yml` container run `31375702141`**, not as
  figures I eyeballed. Push/PR CI (`ci.yml`) lints and runs the fixture tests
  only; it does not restore the store, so these anchors SKIP there — measured on
  a clean checkout with no store, `Rscript tests/testthat.R` gives **0 fail, 26
  skip, every store-backed anchor among the skips** (705 passing expectations at
  the time of writing; the pass count moves whenever a test is added, so the
  falsifiable claim is the 0-fail / 26-skip structure, not the total). This is a
  local measurement on a clean tree, not a figure from a committed run. This is the most
  persuasive part of the repo, and it is persuasive because of the container
  run, not because of the badge.

- **The failing one is pinned, not silenced.** The 4-cluster methylation check
  comes back RED: the merged methylation partition tracks the HM27/HM450 assay
  split (ARI 0.583 against a 0.25 veto). I am careful about what that red means
  — TCGA's m1–m4 are mRNA *expression* subtypes and no four DNA-methylation
  strata are published for KIRC, so `k = 4` is my design choice and the red is a
  finding about this assay merge, not a failed replication. That container run holds an expected-failure ledger — recorded as
  `expected red : 2 / observed red : 2`, "the only red anchors are the recorded
  m1-m4 negative result" — so a *new* red fails the job, and this anchor
  starting to *pass* also fails the job, because that would mean the
  recorded negative no longer holds. Neither widening a published
  range nor lowering a threshold was on the table.

- **Diagnosing the confound rather than reporting the failure.** I scored every
  MOFA factor against the assay split (run `30911448546`): Factor2 AUC 0.888,
  Factor5 0.818 — the two most methylation-loaded factors are substantially
  batch. So the survival model takes only the platform-clean factors, Factor1
  (0.500) and Factor4 (0.535), chosen **outcome-blind** on cleanliness and
  variance explained, with platform as an adjustment covariate. No batch
  correction: only **3** cases are assayed on both platforms and the probe sets
  differ, so a ComBat-style correction could be neither validated nor trusted.

- **The MOFA2 build blocker solved in the scaffold** — basilisk vs renv vs
  Docker vs reticulate, resolved with an external system env. The bilingual
  environment is genuinely necessary, not decorative.

## The caveats I raise before I'm asked (this is the maturity signal)

- **The held-out C-index is a stage-and-age model.** 0.7486 for Cox (0.7492
  penalised, 0.7524 RSF; run `31375702141`) — but in that fit only `stage_num`
  (p = 9.6e-17) and `age_years` (p = 9.7e-06) are significant, while `Factor1`
  (HR 0.982, p = 0.60) and `Factor4` (HR 0.966, p = 0.24) are not. On this
  cohort the multi-omics factors add nothing detectable to clinical staging. I
  say that before anyone reads 0.75 as an integration result.

- **Held-out is not out-of-sample.** MOFA is fitted on all 524 cases and the
  5000-gene variance filter is computed over all samples before the split. Both
  are outcome-blind — no label leak, the BAP1 task stays non-circular — but the
  test rows helped define the axes and the features, so the figures are
  *unsupervised-transductive*. The recorded optimism (0.0125) bounds the
  supervised component only, and there is no external RCC cohort.

- **The BAP1 AUROC of 0.960 rests on 36 mutants** out of 413. The label is
  external so the task is genuinely non-circular, but a handful of held-out
  positives cannot support a claim stronger than "the expression signature is
  plausibly there".

- **BAP1 survival is underpowered by construction.** HR 1.584 (95 % CI
  0.967–2.595, p = 0.068). Schoenfeld needs ≈ 470 events for 80 % power; the
  mutation subset **records 147, 18 of them in the mutant arm**. Directionally
  consistent with the literature, and neither a confirmation nor a negative
  finding. I also quote snapshot numbers, never live-GDC numbers — the frozen
  2016 MAF and today's masked MAF cover different case sets, and conflating them
  is an easy silent error this project has already made once and corrected.

- **Effective n is 241–524**, not ~530, set by modality intersection (241 and
  the 524/413/417 census come from run `30708943504`; run id only, transcript
  not committed). The survival fit reaches 519 rows with 171 OS events — that is
  the FITTED frame, recorded in `docs/results/module4-run-31375702141.txt`, and
  it is not the same census as the 173 OS events counted over the 522-case
  `colData` cohort. The EPV budget is measured on the 124 *training* events
  (cap 12), and 5 predictors are used.

- **"Regularly updated tool" is scoped to one panel.** The research core is a
  frozen 2016/hg19 snapshot; a weekly re-run reproduces it byte-for-byte. What
  genuinely updates is the live GDC API panel, which a weekly cron is **wired**
  to refresh while asserting that no frozen target rebuilt. I say "wired", not
  "does": that cron has never executed in its current form — the two scheduled
  runs on record (`30795010494`, `31359626479`) ran an earlier placeholder job —
  and it cannot succeed until the `targets-store` release asset it restores
  actually exists, which it does not. And the *published* page only advances
  when the Pages workflow is dispatched, which is manual. So the honest claim is
  "a scheduled refresh of a cached live panel, built and not yet run", not "a
  live site".

- **Single cell is not built.** Module 6 (GSE159115) and its ESTIMATE purity
  gate are specified and non-blocking; no code for either exists yet, so the
  bulk → single-cell purity confound is an unaddressed design risk, not a
  checked one. I would rather say that than imply a check that has never run.

- **CI does not reproduce the analysis.** It lints and runs unit tests on
  subsampled fixtures. The deploying render is *designed* to read a cached
  `_targets` store restored from a release asset — that asset has not been
  produced yet, so the site has never been rendered against real numbers, and a
  local render shows a stated gap wherever a store-backed figure would go. A
  green badge means the tests passed, nothing more.

## If asked "what would you do next?"

In order: an external RCC validation cohort (the single biggest weakness);
implement the ESTIMATE purity gate before anything single-cell touches the
subtypes; and re-run the methylation arm restricted to HM450 alone at n = 241 to
see whether four clusters survive when the assay confound is removed by design
rather
than adjusted for — accepting the cohort loss as the price of the answer.
