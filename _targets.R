# Top-level targets DAG for the KIRC multi-omics pipeline.
# Phase 0 seeds an environment-verification target; each later module appends
# its own targets here (mae_raw -> ... -> dashboard_site). See docs/architecture.md.

library(targets)

# Load run parameters / feature flags and expose the globals downstream modules
# assume exist at parse time: Module 1's mae_raw gates its tar_cue on
# HEAVY_PULL and Module 6 reads config$singlecell. Defined here once so every
# appended target can reference them without re-reading the YAML. An explicit
# HEAVY_PULL=true env var overrides the config default (used by fixture
# regeneration / local full runs).
config <- yaml::read_yaml("config/params.yml")
HEAVY_PULL <- isTRUE(config$heavy_pull) ||
  identical(tolower(Sys.getenv("HEAVY_PULL", "false")), "true")

# Source project constants and every pure-function module.
source("R/constants.R")
invisible(lapply(
  list.files("R", pattern = "^functions_.*\\.R$", full.names = TRUE),
  source
))

tar_option_set(
  # Module 1 adds MultiAssayExperiment: its S4 `[` / `[[` methods are only
  # dispatchable when the package is ATTACHED (loading the namespace is not
  # enough), and mae_qc / the *_raw assay targets subset the MAE directly.
  packages = "MultiAssayExperiment",
  format = "rds"
)

list(
  tar_target(
    scaffold_env_check,
    list(
      r_version = R.version.string,
      reticulate_python = Sys.getenv("RETICULATE_PYTHON"),
      snapshot_date = SNAPSHOT_DATE,
      heavy_pull = HEAVY_PULL
    )
  ),

  # --- Module 1: ingest + preprocess ---------------------------------------
  # tar_cue(mode = "never") only suppresses RE-runs: a target with no metadata
  # record still builds on a fresh store, so the command itself must be gated
  # or `tar_make()` on a clean clone starts the 15-30 min ExperimentHub pull.
  tar_target(
    mae_raw,
    if (HEAVY_PULL) {
      fn_load_mae()
    } else {
      stop("mae_raw requires HEAVY_PULL=true (or a restored _targets store); ",
           "rerun with HEAVY_PULL=true Rscript -e 'targets::tar_make()'")
    },
    cue = tar_cue(mode = if (HEAVY_PULL) "thorough" else "never")
  ),
  tar_target(mae_qc, fn_qc_mae(mae_raw)),

  # fn_experiment resolves the real curatedTCGAData experiment name
  # ("KIRC_<stub>-20160128"); a bare mae_qc[["<stub>"]] returns NULL silently
  # and only blows up later inside assay().
  tar_target(rna_raw,     SummarizedExperiment::assay(fn_experiment(mae_qc, "RNASeq2GeneNorm"))),
  tar_target(meth27_raw,  SummarizedExperiment::assay(fn_experiment(mae_qc, "Methylation_methyl27"))),
  tar_target(meth450_raw, SummarizedExperiment::assay(fn_experiment(mae_qc, "Methylation_methyl450"))),
  tar_target(cnv_raw,     SummarizedExperiment::assay(fn_experiment(mae_qc, "GISTIC_ThresholdedByGene"))),

  # Same attach-not-namespace rule as MultiAssayExperiment above:
  # minfi::getAnnotation() resolves its tables off the search list, so the
  # annotation package must be ATTACHED. Declared per-target rather than in
  # tar_option_set() because only this target needs it, and a global entry
  # would make every light target (e.g. scaffold_env_check) unrunnable on a
  # machine without the ~200 MB Bioconductor manifest installed.
  # fn_load_methyl_annotation() also attaches it defensively, for callers
  # outside this pipeline (testthat, bare Rscript).
  tar_target(
    methyl_anno,
    fn_load_methyl_annotation(),
    packages = c(tar_option_get("packages"), METHYL_ANNO_PKG)
  ),

  tar_target(common_ids, fn_intersect_cases(list(
    fn_harmonise_ids(colnames(rna_raw)),
    union(fn_harmonise_ids(colnames(meth27_raw)),
          fn_harmonise_ids(colnames(meth450_raw))),
    fn_harmonise_ids(colnames(cnv_raw))
  ))),

  # "correct n recorded" — enforce the cohort contract (spec section 2). The
  # verified snapshot value is 524 (COHORT_SIZES$rna_methyl_cnv); [520, 535] is
  # the drift band around it.
  tar_target(cohort_n, {
    n <- length(common_ids)
    stopifnot(n >= COHORT_MIN, n <= COHORT_MAX)
    n
  }),

  tar_target(rna_mat, fn_top_variable(
    fn_align_samples(fn_log2_normalise_rna(rna_raw), common_ids),
    N_TOP_GENES
  )),

  # Restrict to the CpGs the merge would keep anyway and drop SNP/sex probes
  # BEFORE the M-value transform. fn_beta_to_mvalue realises the HDF5-backed
  # HM450 DelayedMatrix (~485k x 320) and then allocates several full-size
  # copies (pmax/pmin/log2), ~95% of which the merge immediately discards.
  # Both steps are row-wise, so the result is identical (see test-preprocess.R
  # "probe filtering ... commute with the M-value transform").
  tar_target(methyl_common_cpgs,
             intersect(rownames(meth27_raw), rownames(meth450_raw))),
  tar_target(methyl_merged, fn_merge_methyl_platforms(
    fn_beta_to_mvalue(fn_drop_bad_probes(
      meth27_raw[methyl_common_cpgs, , drop = FALSE], methyl_anno)),
    fn_beta_to_mvalue(fn_drop_bad_probes(
      meth450_raw[methyl_common_cpgs, , drop = FALSE], methyl_anno))
  )),
  tar_target(methyl_mat, fn_top_variable(
    fn_align_samples(methyl_merged, common_ids),
    N_TOP_CPGS
  )),

  tar_target(cnv_mat, fn_align_samples(fn_prep_cnv(cnv_raw), common_ids)),

  tar_target(mut_annot, fn_extract_mutation_status(mae_qc, DRIVER_GENES))
)
