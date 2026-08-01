# All magic values for the pipeline live here (coding-style: no magic numbers).

# --- Frozen data snapshot (curatedTCGAData 1.34.0) ---
SNAPSHOT_DATE <- "20160128"
GENOME_BUILD <- "hg19"

# --- ccRCC somatic driver genes (mutation used as annotation only, spec 6a) ---
DRIVER_GENES <- c("VHL", "PBRM1", "SETD2", "BAP1", "MTOR", "KDM5C")

# --- Published ccRCC mutation-frequency ranges (fraction of tumours) ---
# Anchors for fn_check_mutation_freq (TCGA KIRC, Nature 2013; COSMIC).
MUTATION_FREQ_RANGES <- list(
  VHL   = c(0.40, 0.60),
  PBRM1 = c(0.30, 0.45),
  SETD2 = c(0.08, 0.15),
  BAP1  = c(0.08, 0.15)
)

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
