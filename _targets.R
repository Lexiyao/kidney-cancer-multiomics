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
  # Module 2's MOFA2 / reticulate / SNFtool are deliberately NOT listed here:
  # a global entry makes EVERY target — including scaffold_env_check and all of
  # Module 1 — unbuildable on a machine without them ("could not find packages
  # MOFA2, reticulate in library paths"), which is the same rule the methyl_anno
  # target documents below. They are declared per-target instead.
  packages = c("MultiAssayExperiment"),
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

  # rna_full is the log2-normalised, cohort-aligned expression matrix BEFORE the
  # top-variable filter (~20500 genes x 524). Split out of rna_mat because a
  # PUBLISHED marker panel must not be subject to a data-driven feature filter:
  # Module 3's ccA/ccB check scores the Brannon 2010 proxy panels, and any
  # low-variance member (EPAS1, KDR and FLT1 are the likely casualties) would
  # otherwise be silently dropped by fn_top_variable, degrading a 6-vs-6
  # comparison to a 2-vs-2 one with nothing in the result to say so.
  # rna_mat is unchanged — the same expression, now derived from rna_full.
  tar_target(rna_full,
             fn_align_samples(fn_log2_normalise_rna(rna_raw), common_ids)),
  tar_target(rna_mat, fn_top_variable(rna_full, N_TOP_GENES)),

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

  # Per-sample assay platform, in the column order of methyl_mat. Module 3's
  # m1-m4 check needs it because methyl_merged is cbind(HM27, HM450) with NO
  # batch correction, so the strongest axis in the merged matrix is the assay
  # rather than the biology, and a k-means partition that merely reproduces the
  # platform would otherwise report a green m1-m4 verdict.
  #
  # The rule below MIRRORS fn_align_samples exactly: it keeps the FIRST column
  # per patient (`!duplicated`) and methyl_merged puts HM27 first, so a case
  # assayed on both platforms is carried by its HM27 column. colnames(methyl_mat)
  # is `common_ids` by construction, hence the names.
  tar_target(methyl_platform, {
    hm27_ids <- fn_harmonise_ids(colnames(meth27_raw))
    p <- factor(ifelse(common_ids %in% hm27_ids, "HM27", "HM450"),
                levels = METHYL_PLATFORMS)
    names(p) <- common_ids
    stopifnot(length(unique(p)) == length(METHYL_PLATFORMS))
    p
  }),

  # How many cohort cases were assayed on BOTH platforms. This is the SOLE
  # stated justification for refusing a batch correction ("only 3 cases overlap,
  # so ComBat could not be validated"), and until this target existed it was
  # recorded nowhere and computed nowhere — asserted in six places, derivable
  # from none of them. It CANNOT be corroborated by the 214/310 split either:
  # methyl_platform's own rule labels every dual-assayed case HM27, so the split
  # is identical for any overlap count.
  #
  # Cheap (two colnames vectors and an intersect), so it costs nothing to make
  # the claim executable. Restricted to `common_ids` because that is the cohort
  # the refusal is about.
  tar_target(methyl_platform_overlap, {
    both <- intersect(fn_harmonise_ids(colnames(meth27_raw)),
                      fn_harmonise_ids(colnames(meth450_raw)))
    length(intersect(both, common_ids))
  }),

  tar_target(cnv_mat, fn_align_samples(fn_prep_cnv(cnv_raw), common_ids)),

  tar_target(mut_annot, fn_extract_mutation_status(mae_qc, DRIVER_GENES)),

  # --- Module 3 prerequisite: shared clinical survival frame ----------------
  # Derived once from colData(mae_qc); Module 4's survival model MUST reuse
  # this target instead of building its own survival_df, so the survival
  # columns (sample_id / os_time / os_event) stay consistent across the DAG.
  # os_event = 1 for deceased patients; os_time uses days_to_death for events
  # and days_to_last_followup for censored cases.
  #
  # TWO DEPARTURES from the plan's literal block, both required for the BAP1
  # positive control to be capable of failing at all:
  #
  #  1. The death set is VITAL_STATUS_DEAD_VALUES (= dead / deceased / "1"),
  #     not a literal c("dead", "deceased"). The narrower set is contradicted
  #     by this repo's own MEASURED census (run 30708943504) and by the design
  #     spec, and on a snapshot storing vital_status as 0/1 it matches nothing
  #     -> zero events -> a survival anchor that cannot fail.
  #  2. The required colData columns are checked up front. `cd$days_to_death`
  #     on a snapshot that spells a column differently returns NULL, and the
  #     ifelse() below then silently turns a whole arm into NA — censored cases
  #     would vanish and the anchor would be fitted on deaths only. Fail loudly.
  #
  # IDs are harmonised to patient barcodes because mut_annot, common_ids and
  # every aligned matrix in this DAG are keyed that way; an unharmonised
  # rowname would join to nothing and merge() would silently return 0 rows.
  #
  # SCOPE: this frame covers every case in colData (536), NOT just the 524-case
  # main cohort. fn_check_bap1_survival inner-joins it to mut_annot (417), so
  # the BAP1 anchor is evaluated on the mutation subset by construction. Module
  # 4 must restrict to `common_ids` itself before fitting the survival model.
  #
  # PLATFORM COVARIATE (fn_attach_platform, R/functions_clinical.R). The frame
  # gains a FOURTH column, `platform`, a factor over METHYL_PLATFORMS carrying
  # the HM27/HM450 assay of each case; sample_id / os_time / os_event are
  # untouched in name, order and value, so Module 4's survival_df contract and
  # fn_check_bap1_survival (which selects those three columns explicitly) are
  # unaffected. It is here rather than in a parallel frame so that exactly one
  # clinical table exists in the DAG and the covariate cannot drift away from
  # the outcome it is fitted beside.
  #
  # It is a COVARIATE, not a correction: too few cases are assayed on both
  # platforms (the figure "3" has been stated but is NOT in any committed
  # transcript — the methyl_platform_overlap target above now computes it, and
  # the next container run records it) and the probe sets differ, so ComBat
  # could not be validated on this snapshot; the uncorrected batch stays in the
  # data and Phase 4 adjusts for it instead
  # (measured platform structure: run 30911448546, see fn_attach_platform).
  # NOTHING downstream of this changes a threshold — the m1-m4 anchor is still
  # red at platform_ari 0.583 against the unchanged 0.25 ceiling.
  #
  # `platform` is NA outside the 524-case main cohort, because methyl_platform
  # is only defined there and inventing an assay for a case with no methylation
  # would be a fabrication; see fn_attach_platform for the full argument.
  tar_target(
    clinical,
    {
      cd <- as.data.frame(MultiAssayExperiment::colData(mae_qc))
      required_cols <- c("vital_status", "days_to_death", "days_to_last_followup")
      absent <- setdiff(required_cols, colnames(cd))
      if (length(absent) > 0L) {
        stop("colData(mae_qc) lacks required survival columns: ",
             paste(absent, collapse = ", "))
      }
      days_to_death    <- suppressWarnings(as.numeric(cd$days_to_death))
      days_to_followup <- suppressWarnings(as.numeric(cd$days_to_last_followup))
      os_event <- as.integer(
        tolower(as.character(cd$vital_status)) %in% VITAL_STATUS_DEAD_VALUES
      )
      os_time  <- ifelse(os_event == 1L, days_to_death, days_to_followup)
      sample_id <- fn_harmonise_ids(rownames(cd))
      stopifnot(!anyDuplicated(sample_id))
      fn_attach_platform(
        data.frame(
          sample_id = sample_id,
          os_time   = os_time,
          os_event  = os_event,
          stringsAsFactors = FALSE
        ),
        methyl_platform
      )
    }
  ),

  # --- Module 2: integrate --------------------------------------------------
  # MOFA2 is the MAIN integration; SNF is a cheap second opinion and the two
  # partitions are compared by ARI (`concordance`) instead of the old
  # "consensus clustering" framing. mut_annot is an EXTERNAL LABEL for factor
  # interpretation only — it is never one of the three MOFA views.
  #
  # Same attach-not-namespace, declare-per-target rule as methyl_anno above.
  # mofa_model needs MOFA2 + reticulate to TRAIN; mofa_factors / mofa_varexp
  # need the MOFA2 namespace so the RDS-stored S4 object deserialises and
  # methods::is(x, "MOFA") resolves. Everything downstream of them works on
  # plain matrices and factors, so it stays light.
  tar_target(
    mofa_model,
    fn_run_mofa(list(RNA = rna_mat, Methylation = methyl_mat, CNV = cnv_mat)),
    packages = c(tar_option_get("packages"), "MOFA2", "reticulate")
  ),
  tar_target(mofa_factors, fn_extract_factors(mofa_model),
             packages = c(tar_option_get("packages"), "MOFA2")),
  tar_target(mofa_varexp,  fn_variance_explained(mofa_model),
             packages = c(tar_option_get("packages"), "MOFA2")),
  tar_target(subtypes_mofa, fn_assign_subtypes(mofa_factors)),
  tar_target(
    snf_clusters,
    fn_run_snf(list(RNA = rna_mat, Methylation = methyl_mat, CNV = cnv_mat)),
    packages = c(tar_option_get("packages"), "SNFtool")
  ),
  # Per-factor platform association, in the DAG rather than in a manually
  # triggered workflow. This is the BASIS on which Phase 4 selects its
  # predictors (PLATFORM_CLEAN_MOFA_FACTORS, R/constants.R), and until this
  # target existed that basis was hard-coded numbers in comments citing run
  # 30911448546 — nothing executable could tell you whether Factor1 and Factor4
  # are STILL clean after a MOFA seed / version / upstream-matrix change, and
  # re-adding Factor2 to the predictor vector would have broken no check.
  # Contrast subtypes_mofa, whose platform-cleanliness has been falsifiable at
  # run time since fn_check_subtype_platform.
  tar_target(factor_platform,
             fn_factor_platform_association(mofa_factors, methyl_platform)),

  tar_target(concordance, fn_cluster_concordance(subtypes_mofa, snf_clusters)),
  tar_target(mutation_factor_annot, fn_annotate_mutation(mofa_factors, mut_annot)),

  # --- Module 3: sanity-check positive controls (credibility anchor) --------
  # Literature-anchored ccRCC checks, each returning a structured pass/fail
  # object (spec section 7), plus one platform-cleanliness guard on the subtype
  # assignment. Computed ONCE on the frozen real data and cached;
  # tests/testthat/test-sanity.R reads this target and asserts it against the
  # published literature. Light on packages: every fn_check_* works on plain
  # matrices/data.frames, so no per-target `packages` entry is needed.
  tar_target(
    sanity_results,
    list(
      mutation_freq  = fn_check_mutation_freq(mut_annot),
      bap1_survival  = fn_check_bap1_survival(clinical, mut_annot),
      # methyl_platform is REQUIRED here: without it the m1-m4 verdict cannot
      # distinguish the published biology from the uncorrected HM27/HM450 batch.
      methyl_strata  = fn_check_methyl_strata(methyl_mat,
                                              platform = methyl_platform),
      # rna_full, NOT rna_mat: the published ccA/ccB panels must not pass
      # through the top-5000-variable filter (see the rna_full target above).
      ccab_signature = fn_check_ccab_signature(rna_full),
      # The counterpart to methyl_strata, and the one platform result that came
      # back CLEAN (run 30911448546: ARI 0.0058). Individual MOFA factors carry
      # heavy platform information (Factor2 separates HM27 from HM450 at
      # AUC 0.888), so a subtype assignment independent of the assay is a real,
      # falsifiable property rather than a foregone conclusion — pinned here so
      # a future change to fn_assign_subtypes or to the factor set cannot break
      # it silently.
      subtype_platform = fn_check_subtype_platform(subtypes_mofa,
                                                   methyl_platform)
    )
  )
)
