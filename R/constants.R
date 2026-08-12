# All magic values for the pipeline live here (coding-style: no magic numbers).

# --- Frozen data snapshot (curatedTCGAData 1.34.0) ---
SNAPSHOT_DATE <- "20160128"
GENOME_BUILD <- "hg19"

# --- ccRCC somatic driver genes (mutation used as annotation only, spec 6a) ---
DRIVER_GENES <- c("VHL", "PBRM1", "SETD2", "BAP1", "MTOR", "KDM5C")

# --- Published ccRCC mutation-frequency ranges: see PUBLISHED_MUT_FREQ_RANGES
# in the Module 3 block below. The earlier `MUTATION_FREQ_RANGES` here was a
# second, unnamed-bounds copy of the same published anchor with different
# bounds (PBRM1 0.30-0.45, SETD2/BAP1 0.08-0.15) and ZERO consumers. Two names
# for one literature anchor is the TOP_VARIABLE_GENES/N_TOP_GENES drift hazard
# again, so it is DELETED rather than aliased. Do not reintroduce it. ---

# --- Methylation platforms (merged on common CpGs; never HM450 alone) ---
METHYL_PLATFORMS <- c("HM27", "HM450")

# --- Model-complexity cap (spec section 2). EPV_CAP is EVENTS PER VARIABLE, not
# a predictor count: the pipeline derives the predictor cap AT FIT TIME as
# floor(observed_OS_events / EPV_CAP) from the events actually present in the
# fitted cohort. It is deliberately NOT hard-coded here.
#
# The survival model runs on the MAIN cohort (RNA ∩ Methyl(any) ∩ CNV), not the
# mutation subset; an earlier comment here reasoned from n~270 and an assumed
# 90-110 OS events, which was both the wrong cohort and a number nobody had
# measured. The event count is no longer pending: it is MEASURED.
#
# MEASURED, NOT ASSUMED: GitHub Actions run 30708943504, 2026-08-01, inside
# bioconductor/bioconductor_docker:RELEASE_3_23, from colData(mae) on the FROZEN
# curatedTCGAData 2.0.1 KIRC snapshot 20160128 via this repo's own fn_load_mae /
# fn_qc_mae / fn_harmonise_ids. Event = vital_status in {dead, deceased, 1};
# time = days_to_death if event else days_to_last_followup.
#
#   MAIN cohort n = 524 (522 usable): 173 OS events, 33.1% event rate,
#     median follow-up 1188 d  ->  EPV-10 predictor cap 17
#   MAIN, 5-year restricted (1825 d): 148 events            ->  cap 14
#   (all cases n = 536: 177 events, cap 17; + Mutation n = 413: cap 14)
#
# The measurement LICENSES a cap of 14-17; it does not oblige the design to
# spend it. Phase 4 wires 5 predictors, comfortably under both, and the
# "never genome-wide feature selection" rule is unchanged. ---
EPV_CAP <- 10L

# Minimum OS events before the Module 3 BAP1 positive control will fit a Cox
# model at all. This is a REFUSAL FLOOR, not an anchor: it never touches the
# published direction (HR > 1) or SANITY_MAX_P, so it cannot turn a failing
# result green — it only stops the check reporting a verdict it has no
# information for.
#
# Why it is needed: survival::coxph does NOT error on a degenerate design.
# MEASURED on 300 constructed cases, zero OS events returns hr = NA, p = NA and
# `pass = NA` with no error and no warning, and a SINGLE event returns
# hr = 3.7e+09, ci = [0, Inf], p = 0.999 and `pass = TRUE` — a fit carrying no
# information reporting GREEN. That is the silent-green failure mode the
# VITAL_STATUS_DEAD_VALUES comment below describes, and the guard belongs in the
# check as well as in the decode.
#
# 20 is far below the MEASURED cohort (173 OS events among 522 usable
# main-cohort cases, run 30708943504), so a legitimate run can never trip it.
MIN_OS_EVENTS <- 20L

# --- Feature-selection sizes: see N_TOP_GENES / N_TOP_CPGS in the Module 1
# block below. They are the single source of truth used by _targets.R; do not
# reintroduce a second pair here.

# --- Cohort sizes at case level (primary tumours only, barcode 14-15 == "01").
# MEASURED, NOT ASSUMED: verified 2026-07-31 by a real HEAVY_PULL run (GitHub
# Actions run 30642823359, bioconductor/bioconductor_docker:RELEASE_3_23) using
# this repo's own fn_load_mae / fn_qc_mae / fn_harmonise_ids against the FROZEN
# curatedTCGAData 2.0.1 KIRC snapshot 20160128 that the pipeline consumes.
#
# These are snapshot numbers, NOT live-GDC numbers, and the two disagree. The
# earlier values here (528 / 370 / 374) came from the current GDC API, where
# somatic calls are the harmonised "Masked Somatic Mutation" set; the 2016
# legacy MAF in the snapshot covers MORE cases (417 vs 374), while snapshot
# GISTIC coverage is slightly smaller (528 vs 532 cases). Any figure quoted from
# the live API therefore does not describe this pipeline.
#
# Honest effective n is 241-524: 524 with mutation dropped, 413 with mutation
# required, and only 241 if methylation is (wrongly) restricted to HM450 alone.
#
# Documented expected values; the build-time tolerance band that `cohort_n`
# actually asserts against is COHORT_MIN / COHORT_MAX in the Module 1 block.
# That band [520, 535] brackets the verified 524, so the contract held.
#
# RE-VERIFIED with ZERO DRIFT by run 30708943504 (2026-08-01): RNA ∩ Methyl(any)
# ∩ CNV = 524, + Mutation = 413, RNA ∩ HM450 ∩ CNV ∩ Mutation = 241. ---
COHORT_SIZES <- list(
  rna_methyl_cnv = 524L,          # RNA ∩ Methyl(HM27 ∪ HM450) ∩ CNV: main cohort
  rna_methyl_cnv_mutation = 413L, # ... ∩ Mutation: mutation-annotated subset
  mutation_subset = 417L          # Mutation alone (2016 legacy MAF, not masked)
)

# --- Module 1: ingest + preprocess ---------------------------------------
CURATED_CANCER      <- "KIRC"                    # TCGA kidney renal clear cell
CURATED_VERSION     <- "2.0.1"                   # curatedTCGAData 1.34.0 data version
CURATED_ASSAYS      <- c(                        # exact ExperimentHub assay stubs
  "RNASeq2GeneNorm",                             # RSEM upper-quartile normalised (NOT counts)
  "Methylation_methyl27",
  "Methylation_methyl450",
  "GISTIC_ThresholdedByGene",                    # CNV: gene-level thresholded (-2..2)
  "Mutation"
)
PRIMARY_TUMOUR_CODE <- "01"                      # TCGA barcode chars 14-15 = primary solid tumour
PATIENT_BARCODE_LEN <- 12L                       # TCGA-XX-XXXX
SAMPLE_BARCODE_LEN  <- 15L                       # TCGA-XX-XXXX-01
SEX_CHROMOSOMES     <- c("chrX", "chrY")         # dropped from methylation
METHYL_ANNO_PKG     <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"  # minfi manifest
MVALUE_CLAMP_EPS    <- 1e-3                       # bound beta into [eps, 1-eps] before logit
N_TOP_GENES         <- 5000L                     # top-variable RNA genes (Gaussian view)
N_TOP_CPGS          <- 5000L                     # top-variable merged CpGs (Gaussian view)
COHORT_MIN          <- 520L                      # RNA ∩ Methyl(any) ∩ CNV lower bound (524)
COHORT_MAX          <- 535L                      # ... upper bound (verified n = 524)
SILENT_CLASSES      <- c(                        # variant classes treated as non-driver events
  "Silent", "Intron", "IGR", "RNA",
  "3'UTR", "5'UTR", "3'Flank", "5'Flank"
)

# --- Module 2: integration (MOFA2 / SNFtool) ---
MOFA_N_FACTORS <- 15L        # upper bound; MOFA prunes inactive factors
MOFA_MAXITER   <- 1000L      # convergence_mode = "fast" stops earlier
MOFA_SEED      <- 42L
K_SUBTYPES     <- 4L         # KIRC has 4 documented methylation/expression strata
SUBTYPE_SEED   <- 42L
SNF_K          <- 20L        # SNF K-nearest-neighbours
SNF_ALPHA      <- 0.5        # SNF affinity hyperparameter (sigma)
SNF_T          <- 20L        # SNF diffusion iterations
MIN_MUT_ANNOT_SAMPLES <- 50L # guard: mutation subset must overlap factors
MIN_SNF_COMPLETE_FRAC <- 0.5 # guard: fraction of SNF features that must be finite

# --- Module 3: sanity-check positive controls -------------------------------

# ccRCC somatic driver panel (TCGA KIRC, Nature 2013; PBRM1, Varela 2011).
# COLLISION RESOLUTION: `DRIVER_GENES` (defined at the top of this file) is the
# CANONICAL panel — fn_extract_mutation_status, fn_annotate_mutation, _targets.R
# and the fixture generator all consume it. The plan's Phase-3 name is kept as
# an ALIAS, not a second literal, so the two can never drift apart.
DRIVER_GENE_PANEL <- DRIVER_GENES

# Published somatic mutation-frequency ranges for ccRCC (fraction of tumours).
# Sources: TCGA KIRC (Nature 2013), COSMIC, Ricketts et al. (Cell Rep 2018).
# INDEPENDENT anchor: these bounds are set from the literature and must never
# be widened to make an observed frequency fall inside them.
PUBLISHED_MUT_FREQ_RANGES <- list(
  VHL   = c(low = 0.40, high = 0.60),
  PBRM1 = c(low = 0.28, high = 0.45),
  SETD2 = c(low = 0.08, high = 0.18),
  BAP1  = c(low = 0.06, high = 0.18)
)

# ccA/ccB proxy signature (Brannon et al. 2010). Simplified marker panels, NOT
# the full ClearCode34 classifier: ccB = proliferative / worse OS, ccA =
# angiogenic. Limitation stated openly in dashboard + talking-points.
CCB_PROLIFERATION_MARKERS <- c("MKI67", "CCNB1", "CCNB2", "BUB1", "TOP2A", "FOXM1")
CCA_ANGIOGENESIS_MARKERS  <- c("VEGFA", "CA9", "EPAS1", "ANGPT2", "KDR", "FLT1")

# Number of TCGA KIRC DNA-methylation strata (m1-m4).
METHYL_N_STRATA <- 4L

# Maximum adjusted Rand index between the methylation k-means partition and the
# ASSAY PLATFORM. methyl_mat is cbind(HM27, HM450) with NO batch correction
# (fn_merge_methyl_platforms only column-binds), and platform is the strongest
# single axis in merged 27k/450k M-values — so without this term the green light
# "TCGA KIRC methylation resolves into m1-m4 strata" could be produced entirely
# by the assay rather than by the published biology.
#
# MEASURED on constructed data (500 CpGs, 200 HM27-like + 324 HM450-like
# samples): with iid noise, NO biological strata and only a per-platform mean
# offset, the check returned pass = TRUE from offset 1.5 SD upward
# (silhouette 0.122 / 0.164 / 0.222 at 1.5 / 2.0 / 3.0 SD, kruskal p = 2.7e-80),
# with ARI(cluster, platform) = 0.504.
#
# CALIBRATION (four true strata spread evenly over both platforms): when k-means
# recovered the true strata (ARI vs truth 1.000) the platform ARI was -0.003;
# when it locked onto the assay instead (ARI vs truth 0.299) the platform ARI
# was 0.532. The two regimes are cleanly bimodal, and 0.25 sits in the gap.
#
# This term makes the check STRICTER: it can only turn a green verdict red.
SANITY_MAX_PLATFORM_ARI <- 0.25

# Minimum published markers that must survive per ccA/ccB panel. A REFUSAL
# FLOOR, not an anchor: it never touches SANITY_MIN_SILHOUETTE, SANITY_MAX_P or
# the required direction of the ccA/ccB opposition, so it cannot turn a failing
# check green. The previous floor of 2 let the "Brannon 2010 proxy" degrade to a
# 2-gene-vs-2-gene comparison while returning an object indistinguishable from
# the full 6-vs-6 run — and MEASURED on 20 structureless matrices, the false-
# green rate rises as the panel shrinks (6+6: 0/20 pass; 2+2: 1/20).
SANITY_MIN_MARKERS_PER_PANEL <- 4L

# Sanity-check decision thresholds.
#
# SANITY_MIN_SILHOUETTE applies to the HIGH-DIMENSIONAL methylation clustering
# only (a ~4568-dimensional 4-means, where iid noise gives ~0.005 — see the
# fn_check_methyl_strata negative control). It must NOT be reused for the 2-D
# ccA/ccB clustering: the two have completely different null distributions, and
# at 0.10 the ccA/ccB silhouette assertion could not fail on any input the
# function accepts, while still being carried in `pass` and shown on the
# dashboard as if it were evidence.
SANITY_MIN_SILHOUETTE <- 0.10

# Separately CALIBRATED threshold for the 2-D (ccB-score, ccA-score) 2-means.
# In two dimensions k-means always splits a blob, so this statistic is large
# under the null and needs its own floor, set ABOVE the measured null ceiling.
#
# MEASURED, NOT ASSUMED (R 4.6.0, this repo's own fn_check_ccab_signature,
# 2026-08-02): 200 structureless matrices at each of n = 40 / 100 / 524 — all 12
# published markers present, iid noise, NO ccA/ccB structure whatsoever:
#
#   n =  40  silhouette min 0.284  median 0.350  q99 0.434  MAX 0.466
#   n = 100  silhouette min 0.295  median 0.334  q99 0.379  MAX 0.387
#   n = 524  silhouette min 0.298  median 0.318  q99 0.341  MAX 0.349
#
# Genuine ccA/ccB opposition, weakest effect tested (a 1-SD-unit shift between
# the two halves, n = 200, 40 replicates): silhouette 0.490-0.562.
#
# 0.50 therefore sits above the null ceiling (0.466) and below the weakest real
# structure (0.490). It makes the check STRICTER, never more permissive.
SANITY_MIN_SILHOUETTE_2D <- 0.50

SANITY_MAX_P          <- 0.05
SANITY_SEED           <- 42L

# --- BAP1 survival: what this cohort can and cannot be asked to show --------
#
# Published effect band for BAP1-mutant vs wild-type ccRCC overall survival.
# Sources: Kapur et al. (Lancet Oncol 2013), Hakimi et al. (Eur Urol 2013),
# TCGA KIRC (Nature 2013). Reported adjusted/unadjusted HRs across those
# cohorts sit roughly in 1.2-3.0; the anchor requires the point estimate to be
# BOTH in the published DIRECTION (HR > 1) and of a plausible MAGNITUDE.
#
# This band is a TIGHTENING, not a loosening. The previous anchor asserted only
# `hr > 1`, which an HR of 1.0001 or of 40 (a merge blow-up, or a fit on a
# handful of cases) would both satisfy. An INDEPENDENT literature anchor: never
# widen it to make an observation pass.
PUBLISHED_BAP1_HR_RANGE <- c(low = 1.2, high = 3.0)

# Target power for the Schoenfeld event-requirement arithmetic below.
SURVIVAL_TARGET_POWER <- 0.80

# Widest 95% CI the BAP1 anchor will accept, as the ratio ci_high / ci_low.
#
# THIS IS NOT A SIGNIFICANCE DEMAND. `ci_high / ci_low` = exp(2 * 1.96 * SE) is
# a pure PRECISION statistic: it says nothing about where the interval sits
# relative to 1, so it cannot be satisfied or violated by moving the point
# estimate, and it can never turn a red check green.
#
# WHY IT EXISTS. Removing `expect_gt(ci_low, 1)` (see the block above) also
# removed the only assertion that the fit carried INFORMATION, and nothing
# replaced it: `hr > 1` plus the published band plus an UPPER bound on the
# mutant event count are all satisfied by a fit on a collapsed exposed arm.
# VERIFIED by replaying the anchor bodies against a stub: hr 1.6, n_mutant 5,
# n_events_mutant 2, CI 0.35-7.3 (ratio 21x), p 0.51 passed every BAP1 anchor.
#
# WHERE 5 COMES FROM. The interval ratio MEASURED on the frozen snapshot (run
# 30840373033) is 2.595 / 0.967 = 2.68, i.e. SE(log HR) = 0.252. A two-group
# Cox with d1 events in the exposed arm and d0 in the reference has
# SE(log HR) ~= sqrt(1/d1 + 1/d0), so a ceiling of 5 (SE 0.41) is slack by a
# wide margin for any fit this cohort can produce, while the 21x degenerate
# case above is rejected by a factor of four. Tighten it if the fit ever gets
# more precise; never widen it to admit a fit that got worse.
BAP1_MAX_CI_RATIO <- 5
#
# WHY THE ANCHOR NO LONGER DEMANDS p < 0.05 FOR BAP1 — the arithmetic, in full,
# so that this reads as a MIS-SPECIFIED REQUIREMENT being corrected and not as a
# threshold being quietly relaxed to turn a red check green.
#
# MEASURED on the frozen snapshot (GitHub Actions run 30840373033), i.e. printed
# in docs/results/phase3-anchors-run-30840373033.txt and checkable there:
#   HR 1.584, 95% CI 0.967-2.595, p = 0.0677, n = 417, BAP1-mutant 8.63%
#   (36 of 417 cases).
#
# ALSO MEASURED, and this settles what an earlier version of this block called
# "DERIVED, NOT RECORDED". Run 31375702141 (verify-module2.yml LEVEL 3),
# transcribed at docs/results/module4-run-31375702141.txt:
#   n_events = 147, n_mutant = 36, n_events_mutant = 18,
#   events_required = 470.4, underpowered = TRUE.
#
#   The reconstruction this block used to carry — 138 events, ~12 in the mutant
#   arm — is superseded and must not be quoted. The recorded 18 is also what the
#   recorded interval implied: SE(log HR) = log(2.5946 / 0.9670) / (2 * 1.96)
#   = 0.2518 reproduces p = 0.0677 exactly, and for a two-arm Cox
#   SE ~= sqrt(1/d1 + 1/d0) gives d1 ~ 18, not 12. The tension is resolved in
#   favour of the measurement, not of the arithmetic.
#
# Schoenfeld's requirement for a two-group Cox comparison is
#   d = (z_{1-alpha/2} + z_{power})^2 / (p * (1 - p) * log(HR)^2)
# with p the exposed fraction. At alpha 0.05 two-sided, 80% power, p = 0.0863:
#
#   HR 1.584 (the OBSERVED effect)          -> d = 470 events
#   HR 1.2   (the weakest published effect) -> d = 2993 events
#   HR 3.0   (the strongest)                -> d = 82 events
#
# The cohort supplies 147, 18 of them in the mutant arm (RECORDED, run
# 31375702141). A 417-case ccRCC series with an 8.6%-prevalence marker of this
# effect size is STRUCTURALLY incapable of reaching alpha 0.05: it would need
# ~3.2x its event count. Requiring `p < 0.05` and `ci_low > 1` of
# it was therefore a requirement no correct pipeline running on this snapshot
# could ever satisfy — it tested the size of TCGA KIRC, not the correctness of
# this code. The DIRECTION and MAGNITUDE are what the literature anchors and
# what this cohort can speak to, so those are the hard requirements; p, CI and n
# are REPORTED and asserted to be well-formed.
#
# This is the ONE re-specification in this correction. Nothing else in the
# Module 3 suite was re-thresholded: SANITY_MAX_PLATFORM_ARI, SANITY_MIN_SILHOUETTE
# and every published range are UNCHANGED, and the m1-m4 anchor stays RED.
#
# The under-power is itself ASSERTED (see the BAP1 anchor in test-sanity.R), so
# the limitation is tested rather than merely written down here. Note what that
# assertion is and is not: it is a statement about the DESIGN — how many events
# an effect of this size needs — and it is NOT a post-hoc-power defence of the
# non-significant p-value. It exists so that if a future cohort ever does supply
# enough events, this comment stops being true and the test says so.

# Maximum adjusted Rand index between the MOFA subtype assignment and the assay
# platform. Distinct from SANITY_MAX_PLATFORM_ARI (0.25), which governs the
# m1-m4 methylation k-means: that partition IS platform-driven on this snapshot
# (measured ARI 0.583) and its anchor is red. The subtype assignment is the one
# integration output the platform diagnostic cleared, and this constant PINS
# that so it cannot silently break.
#
# MEASURED (GitHub Actions run 30911448546, frozen curatedTCGAData KIRC
# snapshot, 524 cases, HM27 214 / HM450 310): ARI(subtypes_mofa, platform) =
# 0.0058, cross-tab S1 11/9, S2 120/186, S3 33/43, S4 50/72. Under independence
# the expected ARI is 0, so 0.0058 is indistinguishable from chance.
#
# 0.05 is an order of magnitude below the m1-m4 ceiling and ~9x above the
# measured value, i.e. tight enough that a genuine platform-driven subtype
# assignment (individual factors reach AUC 0.888 against platform, so the raw
# material for one is present) trips it, loose enough that resampling noise
# around zero does not. Like SANITY_MAX_PLATFORM_ARI this term can only turn a
# green verdict red; never raise it to make an observation pass.
SUBTYPE_MAX_PLATFORM_ARI <- 0.05

# The MOFA factors Phase 4 may use as survival predictors.
#
# SELECTION RULE — DETERMINISTIC AND OUTCOME-BLIND. Re-running the diagnostic
# and applying the three clauses below mechanically yields this vector; if it
# does not, one of them is wrong and must be fixed, not worked around.
#
#   (1) ELIGIBILITY. A factor is eligible iff its association with the
#       HM27/HM450 split is non-significant after BH correction across ALL
#       factors: q > SANITY_MAX_P in the `factor_platform` target.
#   (2) RANKING. Eligible factors are ranked by TOTAL variance explained,
#       SUMMED across the three views (RNA + Methylation + CNV) in
#       `mofa_varexp`. Summed, not max and not CNV-only: the three views
#       disagree in scale and picking one of them after the fact is a free
#       parameter.
#   (3) COUNT. Take the top N_SURVIVAL_MOFA_FACTORS = 2. The count is fixed IN
#       ADVANCE by the predictor budget (2 factors + age + stage + platform =
#       5 terms, against an EPV-10 cap of 17, 14 at the 5-year horizon), not
#       chosen after looking at the ranking.
#
# SURVIVAL ASSOCIATION IS NEVER CONSULTED. Picking factors by how well they
# predict the outcome and then reporting that model's discrimination on the same
# cohort is selection bias; a held-out split does not repair it.
#
# APPLIED TO run 30911448546 + the Module 2 variance table (run 30718392588):
#   eligible (q > 0.05): Factor1 (AUC 0.500, q 0.993), Factor4 (0.535, q 0.217),
#     Factor7 (0.532, q 0.243), Factor13 (0.523, q 0.389), Factor15 (0.552,
#     q 0.058).
#   ranked by summed variance explained (RNA/Methyl/CNV %):
#     Factor1  2.98 + 3.20 + 12.93 = 19.11
#     Factor4  3.18 + 0.97 +  8.16 = 12.31
#     Factor7  0.72 + 0.34 +  3.65 =  4.71
#     Factor13 0.29 + 0.28 +  2.04 =  2.61
#     Factor15 0.19 + 0.02 +  2.17 =  2.38
#   top 2 -> Factor1, Factor4. Factor15's marginal q (0.058) needs no special
#   rule: it is eligible and simply ranks fifth.
#
# INELIGIBLE, for the record: Factor2 (AUC 0.888, q 2.9e-50 corrected -- the
# transcript's 1.4e-50 was depressed by the impossible Factor6 p = 0, see
# docs/results/platform-diagnosis-run-30911448546.txt) and Factor3 (0.658,
# q 2.7e-09) are the two that a pre-diagnostic predictor set contained; Factor5,
# Factor6, Factor9, Factor8, Factor11, Factor14, Factor10 and Factor12 are also
# significant after BH and likewise excluded.
#
# PINNED, NOT ASSUMED: the `factor_platform` target recomputes clause (1) every
# run and an ANCHOR asserts every name below still satisfies it. That anchor can
# only turn a green verdict RED.
N_SURVIVAL_MOFA_FACTORS <- 2L
PLATFORM_CLEAN_MOFA_FACTORS <- c("Factor1", "Factor4")

# Vital-status values that encode an overall-survival EVENT (death), lower-cased
# before matching. This is the decode the repo has ALREADY VERIFIED on the real
# snapshot: GitHub Actions run 30708943504 read colData(mae)$vital_status with
# `tolower(as.character(v)) %in% c("dead", "deceased", "1")` and obtained 177
# events over 536 cases / 173 over the 524-case main cohort. The design spec and
# this file's EPV block both state the same set.
#
# "1" is LOAD-BEARING and must not be dropped: curatedTCGAData stores
# vital_status numerically on some cohorts, and a {dead, deceased} test on a
# 0/1 column matches NOTHING. That yields zero events, and a survival positive
# control fitted on zero events cannot fail — the exact silent-green failure
# mode the Module 3 suite exists to prevent.
VITAL_STATUS_DEAD_VALUES <- c("dead", "deceased", "1")

# Guard for fn_check_methyl_strata: fraction of CpGs that must be complete.
# stats::kmeans errors on NA/NaN/Inf and the real methyl_mat is only 91.4%
# complete (432/5000 CpGs carry a non-finite value), so incomplete CpGs are
# dropped before clustering — the same resolution already applied to SNF
# (MIN_SNF_COMPLETE_FRAC). This is DATA HYGIENE, not an anchor: it never
# touches METHYL_N_STRATA, SANITY_MIN_SILHOUETTE or SANITY_MAX_P, so it cannot
# make a failing stratum check pass. A matrix gutted by missingness stops
# loudly rather than clustering its remnants into a silent green.
SANITY_MIN_COMPLETE_FRAC <- 0.5

# --- Module 4: survival + classifier model parameters ---
HELDOUT_FRACTION      <- 0.30      # held-out test fraction (spec §6c)
CV_FOLDS              <- 5L        # cv.glmnet / StratifiedKFold folds
MODEL_SEED            <- 20160128L # snapshot date, reused as RNG seed
RSF_NTREE             <- 1000L     # randomForestSRC trees
CALIBRATION_BINS      <- 5L        # risk-group bins for grouped calibration
SURVIVAL_HORIZON_DAYS <- 5 * 365   # 5-year OS calibration horizon
BAP1_LABEL_COL        <- "BAP1"    # column of mut_annot holding BAP1 status

# --- Live GDC API panel (Module 5) -------------------------------------------
# The ONLY live-updating source. The research core stays on the frozen
# curatedTCGAData snapshot (SNAPSHOT_DATE); this panel is refreshed weekly.
#
# THESE COUNTS ARE NOT THE COHORT. Anything this panel returns is a CURRENT
# GDC-portal figure for the whole TCGA-KIRC project; the analysis cohort is
# COHORT_SIZES$rna_methyl_cnv cases from the frozen 20160128 snapshot. The two
# disagree in both directions (the 2016 legacy MAF covers MORE cases than
# today's masked-somatic MAF; CNV and HM450 cover fewer), so a number from here
# must never be quoted as a snapshot number, or vice versa.
GDC_API_BASE   <- "https://api.gdc.cancer.gov"
GDC_PROJECT_ID <- "TCGA-KIRC"
# VERIFIED against https://api.gdc.cancer.gov/cases/_mapping on 2026-08-10.
# The plan's third field was `demographic.gender`, which the current GDC cases
# index does NOT define -- a facet request for it comes back with
# `warnings.facets = "unrecognized values: [demographic.gender]"` and no
# aggregation, which would take the whole live panel down. The field is now
# `demographic.sex_at_birth`.
GDC_FACET_FIELDS <- c(
  "demographic.vital_status",
  "demographic.sex_at_birth",
  "diagnoses.ajcc_pathologic_stage"
)

# --- Module 5: dashboard presentation -----------------------------------------
# Significant digits for the evidence strings on the landing-page sanity table.
# Presentation only: rounding here can never change a pass/fail verdict, which
# is computed in Module 3 and merely reported by Module 5.
SANITY_DETAIL_DIGITS <- 3L
# Attached to any sanity check that flags itself `underpowered`. A green tick
# on such a check means "the direction matches the literature", never "the
# result is confirmed", and the dashboard must say so in the same cell.
SANITY_UNDERPOWERED_NOTE <- paste(
  "UNDERPOWERED: direction consistent with the literature,",
  "the test itself is not conclusive at this event count"
)

# Where a rendering .qmd finds the frozen store. `_quarto.yml` sets
# `execute-dir: file`, so a page's working directory is `dashboard/` and the
# store restored from the release asset sits one level up. Every page reads it
# through fn_dashboard_read(); no page constructs this path itself.
DASHBOARD_STORE_PATH <- "../_targets"

# House accent, kept identical to --kirc-accent in dashboard/styles.css so the
# Plotly traces and the CSS do not drift apart.
DASHBOARD_ACCENT <- "#2c6e91"

# The two methylation assays are the project's largest confound, so they get
# fixed colours: whenever a plot splits by platform it splits the same way on
# every page. HM450 borrows the "live"/warning amber deliberately -- a reader
# should register the split as a caveat, not as decoration.
DASHBOARD_PLATFORM_COLOURS <- c(HM27 = "#2c6e91", HM450 = "#d08700")

# Verdict vocabulary for the per-factor platform annotation on factors.qmd.
# `confounded` is assigned on the SAME q <= SANITY_MAX_P rule that Phase 4 used
# to choose its predictors, so the page cannot describe a factor as clean that
# the model treated as dirty (or the reverse).
FACTOR_VERDICT_CONFOUNDED <- "CONFOUNDED by assay"
FACTOR_VERDICT_CLEAN      <- "platform-clean"
FACTOR_VERDICT_DEGENERATE <- "AUC UNUSABLE (p underflowed to 0)"
FACTOR_VERDICT_UNSCORED   <- "not scored against platform"

# A value box whose underlying target is absent from THIS render. The main
# dashboard is the most quoted surface in the project, so a value box is the
# single easiest place for an unmeasured quantity to enter a conversation as a
# result. An absent number is therefore printed as an absent number, in the
# warning colour, and never as a plausible-looking figure.
DASHBOARD_PENDING_VALUE  <- "not in this render"
DASHBOARD_PENDING_COLOUR <- "warning"

# Rendered-site directory, relative to the repo root. `dashboard/_quarto.yml`
# sets `output-dir: ../_site`, the `dashboard_site` target returns this path as
# its file output, and `.github/workflows/pages.yml` uploads it. One constant so
# those three cannot drift into publishing an empty directory.
DASHBOARD_SITE_DIR <- "_site"
