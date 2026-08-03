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
