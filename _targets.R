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
  ),

  # --- Module 4: survival ---------------------------------------------------
  # SELECTION RULE (load-bearing — read before changing this vector).
  # It is DETERMINISTIC and OUTCOME-BLIND, and it lives in R/constants.R:
  #   (1) ELIGIBILITY: q > SANITY_MAX_P against the HM27/HM450 split in the
  #       `factor_platform` target (the in-DAG recomputation of run
  #       30911448546).
  #   (2) RANKING: eligible factors by TOTAL variance explained SUMMED across
  #       RNA + Methylation + CNV in `mofa_varexp` — summed, not max, not one
  #       chosen view, because the views differ in scale and picking one after
  #       the fact is a free parameter.
  #   (3) COUNT: the top N_SURVIVAL_MOFA_FACTORS = 2, fixed IN ADVANCE by the
  #       predictor budget (2 factors + age + stage + platform = 5 terms
  #       against an EPV-10 cap of 17, 14 at the 5-year horizon), not chosen
  #       after seeing the ranking. Factor15's marginal q (0.058) needs no
  #       special rule: it is eligible and simply ranks fifth.
  # Applied: Factor1 (19.11% summed) and Factor4 (12.31%) — the full ranking,
  # the ineligible factors (Factor2 AUC 0.888, Factor3 0.658, ...) and the
  # provenance are recorded once, on PLATFORM_CLEAN_MOFA_FACTORS.
  #
  # It is NEVER chosen on survival association. Picking factors by how well
  # they predict the outcome and then reporting that model's discrimination on
  # the same cohort is selection bias — the C-index would be optimistic by
  # construction and the held-out split would not repair it. The choice was
  # fixed from the platform diagnostic and the Module-2 variance table BEFORE
  # any survival model was fitted.
  #
  # DEVIATION FROM THE PLAN, RECORDED: the plan's literal block writes
  # `c("Factor1", "Factor4", "age_years", "stage_num", "platform")`. That is a
  # SECOND copy of PLATFORM_CLEAN_MOFA_FACTORS, and the anchor that keeps the
  # selection falsifiable (test-clinical.R, "every factor wired as a survival
  # predictor is platform-clean") asserts against the CONSTANT, not against
  # this target — so an edit to the constant would leave the wired vector
  # stale and the anchor would still pass. Referencing the constant is the same
  # one-source-of-truth rule this Task applies to the OS decode below.
  #
  # `platform` is an ADJUSTMENT COVARIATE, not a finding: the merged
  # HM27 214 / HM450 310 matrix received NO batch correction (too few
  # overlapping cases to validate one — see `methyl_platform_overlap` — and the
  # probe sets differ), so the residual assay effect is modelled explicitly
  # instead of being assumed away. It costs 1 df.
  tar_target(
    survival_predictors,
    c(PLATFORM_CLEAN_MOFA_FACTORS, "age_years", "stage_num", "platform")
  ),

  tar_target(survival_df, {
    cd <- as.data.frame(MultiAssayExperiment::colData(mae_qc))

    # SURVIVAL COMES FROM `clinical` (Module 1/3), NOT from a second decode
    # here. This target used to re-derive vital_status / days_to_death /
    # days_to_last_followup independently and key on RAW `rownames(cd)`, while
    # `clinical` decodes with VITAL_STATUS_DEAD_VALUES and keys on harmonised
    # patient barcodes. That put two divergent OS derivations, with
    # incompatible keys, into the same DAG — and the failure mode where a
    # {dead, deceased} test against a 0/1-coded vital_status yields zero events
    # and an unfalsifiable survival model.
    #
    # Only the NON-survival covariates are read from colData here.
    req_cols <- c("years_to_birth", "pathologic_stage")
    missing_cols <- setdiff(req_cols, colnames(cd))
    stopifnot(length(missing_cols) == 0)

    covars <- data.frame(
      sample_id = fn_harmonise_ids(rownames(cd)),
      age_years = as.numeric(cd$years_to_birth),
      stage_num = as.integer(factor(tolower(trimws(cd$pathologic_stage)),
                                    levels = c("stage i", "stage ii",
                                               "stage iii", "stage iv"))),
      stringsAsFactors = FALSE
    )

    # Rename at THIS boundary and nowhere else: `clinical` is os_time/os_event
    # across the DAG, while the fitters (fn_fit_cox and friends all
    # `stopifnot(all(c("time","status") %in% names(surv_df)))`) and Module 5's
    # km_subtype_df contract on time/status. One rename, one source of truth.
    base <- merge(clinical, covars, by = "sample_id")
    names(base)[names(base) == "os_time"]  <- "time"
    names(base)[names(base) == "os_event"] <- "status"
    stopifnot(sum(base$status == 1L) > 0)  # non-degenerate: at least one event

    fac <- as.data.frame(mofa_factors)
    fac$sample_id <- rownames(mofa_factors)
    merged <- merge(base, fac, by = "sample_id")

    # Platform adjustment covariate — ALREADY PRESENT, inherited from
    # `clinical`. DO NOT RE-JOIN IT HERE. `clinical` carries it "so that exactly
    # one clinical table exists in the DAG and the covariate cannot drift away
    # from the outcome it is fitted beside", and the merge above already brings
    # it in. An earlier draft re-derived it as
    # `merged$platform <- methyl_platform[merged$sample_id]`, silently
    # OVERWRITING the inherited column: the values agreed (same upstream
    # target) so nothing failed, but the survival frame's platform then no
    # longer came from the canonical clinical table and the assertions
    # validated the re-join rather than the inherited column. It was also
    # ORDER-DEPENDENT — move the overwrite below the complete.cases filter and
    # the filter uses a different column from the model.
    #
    # `clinical$platform` is NA outside the 524-case main cohort; the merge
    # with `mofa_factors` above restricts to that cohort, which is why the
    # no-NA assertion can be unconditional here.
    stopifnot("platform" %in% names(merged),
              !any(is.na(merged$platform)),
              nlevels(droplevels(merged$platform)) == 2L)

    merged <- merged[stats::complete.cases(
      merged[, c("time", "status", survival_predictors)]), , drop = FALSE]
    merged <- merged[merged$time > 0, , drop = FALSE]
    stopifnot(sum(merged$status == 1L) > 0)  # events survive the row filtering
    merged
  }),

  # All three arms take the SAME held-out split (same helper, same default
  # seed), so their C-indices are comparable and none of them ever sees the
  # test rows during fitting: cv.glmnet picks lambda inside TRAIN, the forest
  # is grown on TRAIN, and the test part is only ever scored.
  tar_target(cox_fit,
             fn_fit_cox(survival_df, survival_predictors)),
  tar_target(penalised_cox_fit,
             fn_fit_penalised_cox(survival_df, survival_predictors),
             packages = c(tar_option_get("packages"), "glmnet")),
  tar_target(rsf_fit,
             fn_fit_rsf(survival_df, survival_predictors),
             packages = c(tar_option_get("packages"), "randomForestSRC")),

  tar_target(survival_metrics, {
    c_cox <- fn_cindex(cox_fit$test$time, cox_fit$test$status, cox_fit$risk_test)
    c_pen <- fn_cindex(penalised_cox_fit$test$time,
                       penalised_cox_fit$test$status, penalised_cox_fit$risk_test)
    c_rsf <- fn_cindex(rsf_fit$test$time, rsf_fit$test$status, rsf_fit$risk_test)
    # Optimism = apparent (train) minus validated (held-out) C-index for Cox.
    # REPORTED, never subtracted away: subtypes and factors were derived from
    # the same omics, so the apparent figure is optimistic by construction and
    # the gap is the finding, not an embarrassment to be hidden.
    #
    # THE HELD-OUT FIGURE IS NOT FULLY OUT-OF-SAMPLE EITHER, and saying so here
    # is the point: `mofa_factors` is fitted by fn_run_mofa on ALL 524 cases —
    # the same rows fn_split_train_test later calls the TEST set — so the
    # latent axes the model is built on were defined with the held-out rows
    # contributing. MOFA never sees an outcome, so this is NOT a label leak and
    # the non-circularity claim stands; it is UNSUPERVISED-TRANSDUCTIVE
    # optimism. The number reported below therefore bounds the SUPERVISED
    # component ONLY. A fully out-of-sample estimate would mean refitting MOFA
    # (and the top-variable-gene filter) inside every training split, which is
    # out of scope at this cohort size. Mirrored in the README limitations.
    c_cox_train <- fn_cindex(cox_fit$train$time, cox_fit$train$status,
                             cox_fit$risk_train)
    # Calibration for Cox at the 5-year horizon
    # (survival prob = exp(-baseline * exp(lp))). basehaz(centered = TRUE) and
    # predict(type = "lp") share the SAME centring at the training means, so
    # the two halves of this expression are on one scale.
    bh <- survival::basehaz(cox_fit$model, centered = TRUE)
    h0 <- approx(bh$time, bh$hazard, xout = SURVIVAL_HORIZON_DAYS,
                 rule = 2)$y
    pred_surv <- exp(-h0 * exp(cox_fit$risk_test))
    cal <- fn_calibration(cox_fit$test$time, cox_fit$test$status,
                          pred_surv, SURVIVAL_HORIZON_DAYS)
    list(
      cindex = list(cox = c_cox, penalised = c_pen, rsf = c_rsf),
      optimism = list(cox = c_cox_train - c_cox),
      calibration = cal
    )
  }),

  # ---- Module 4: non-circular BAP1 classifier (R+Python via reticulate) ----
  # The label (BAP1 somatic-mutation status) is EXTERNAL to the features
  # (log-transformed normalised expression), so this is the non-tautological
  # counterpart the spec demands in section 6b -- unlike predicting subtypes,
  # which were themselves derived from this expression matrix.
  #
  # `rna_mat` is the top-variable-gene subset. That filter is UNSUPERVISED
  # (variance only, no reference to BAP1), so it selects features without
  # seeing the label and does not leak it. It is, however, computed over ALL
  # samples — including the rows the classifier is later scored on — so
  # `heldout_auroc` is transductive in the same sense the survival C-index is:
  # outcome-blind, hence non-circular, but not fully out-of-sample. Within the
  # Python function the standardiser sits inside the pipeline, so it is
  # refitted per CV fold and on the training split only.
  #
  # reticulate is declared per-target, never in tar_option_set: a global entry
  # makes every unrelated target unbuildable on a machine without it (the same
  # rule mofa_model and methyl_anno document above).
  tar_target(bap1_auroc, {
    reticulate::source_python("python/bap1_classifier.py")
    # Restrict to the n=417 mutation subset and align samples to expression
    # (413 of those cases fall inside the 524-case main cohort).
    common <- intersect(colnames(rna_mat), rownames(mut_annot))
    # The plan's assertion here was `length(common) > 0`, which passes on an
    # overlap of one. Reuse the guard fn_annotate_mutation already applies to
    # the same join instead: a divergent ID space (aliquot barcodes vs
    # harmonised patient IDs) collapses the overlap, and a handful of samples
    # would otherwise reach scikit-learn and produce an AUROC on noise.
    stopifnot(BAP1_LABEL_COL %in% names(mut_annot),
              length(common) >= MIN_MUT_ANNOT_SAMPLES)
    expr_df <- as.data.frame(t(rna_mat[, common, drop = FALSE]))  # samples × genes
    labels  <- as.integer(mut_annot[common, BAP1_LABEL_COL])
    res <- train_bap1_classifier(
      reticulate::r_to_py(expr_df),
      reticulate::r_to_py(labels)
    )
    # source_python(convert = TRUE) -- the default -- already hands back an R
    # list, and reticulate::py_to_r() STOPS with "Object to convert is not a
    # Python object" when it is given one. The guard keeps the plan's explicit
    # conversion for a convert = FALSE session without making the default
    # session error on a value that is already converted.
    if (inherits(res, "python.builtin.object")) reticulate::py_to_r(res) else res
  }, packages = c(tar_option_get("packages"), "reticulate")),

  # --- Module 5: live GDC panel (independent of frozen research core) ---------
  # THE ONLY TARGET IN THIS DAG THAT IS NOT FROZEN. Everything above descends
  # from the curatedTCGAData 20160128 snapshot and is byte-stable across runs;
  # this one re-queries api.gdc.cancer.gov and legitimately changes. It shares
  # NO edge with any research target, which is the point: the "regularly updated
  # tool" claim (spec section 9) is scoped to this panel and nothing else. Any
  # page rendering it must label it live, with `retrieved_at`, and must never
  # let a reader take `total_cases` for the frozen cohort's n.
  #
  # `packages` is declared here rather than inherited, for the same reason
  # MOFA2/reticulate are declared per-target above: the global default attaches
  # MultiAssayExperiment, and the weekly cron (.github/workflows/cron.yml)
  # builds THIS TARGET AND NOTHING ELSE. Inheriting the default would force a
  # Bioconductor install into a job whose entire work is one HTTPS request, and
  # would make the panel refresh fail for a reason that has nothing to do with
  # the panel. httr is all fn_query_gdc() touches.
  tar_target(
    gdc_live_panel,
    fn_query_gdc(),
    packages = "httr",
    cue = tar_cue(mode = "always")   # always re-query so the panel truly updates
  ),

  # --- Module 5: tidy consumption targets for the dashboard -------------------
  # The .qmd pages need rectangular, join-ready shapes. Module 2 emits
  # `subtypes_mofa` as a NAMED FACTOR, Module 3 emits `sanity_results` as a
  # NESTED NAMED LIST, and Module 4's Cox model has NO `subtype` term — so these
  # helper targets adapt those outputs once, here, rather than in every page.
  tar_target(
    subtypes_df,
    data.frame(
      sample_id = names(subtypes_mofa),
      subtype   = unname(as.character(subtypes_mofa)),
      stringsAsFactors = FALSE
    )
  ),
  tar_target(
    km_subtype_df,
    {
      merged <- merge(
        survival_df[, c("sample_id", "time", "status")],
        subtypes_df,
        by = "sample_id"
      )
      merged[!is.na(merged$subtype), , drop = FALSE]
    }
  ),
  tar_target(
    sanity_table,
    do.call(rbind, lapply(names(sanity_results), function(nm) {
      el <- sanity_results[[nm]]
      data.frame(
        check  = if (!is.null(el$label))  el$label  else nm,
        passed = isTRUE(el$pass),
        # DEVIATION FROM THE PLAN, deliberately. The plan takes the element's
        # own `detail` field and falls back to an empty string.
        # NO element of `sanity_results` has a `detail` field, so every row
        # of the landing-page table would render blank. Two rows cannot survive
        # that: BAP1 survival passes its `hr > 1` criterion while being
        # UNDERPOWERED, and m1-m4 fails for a reason `fn_check_methyl_strata`
        # already wrote out in `$message`. A blank cell next to a green tick is
        # how an underpowered result gets read as a confirmation.
        detail = fn_sanity_detail(el),
        stringsAsFactors = FALSE
      )
    }))
  ),

  # --- Module 5: render the full Quarto site from the frozen store ------------
  # The .qmd pages do NOT receive these objects as arguments; each reads the
  # store itself through fn_dashboard_read() (`../_targets`, because
  # `execute-dir: file` puts a page's working directory in dashboard/). Naming
  # the symbols here is what creates the DAG edges, so an upstream change
  # re-renders the site instead of leaving a stale page beside a fresh one.
  #
  # LONGER THAN THE PLAN'S LIST, deliberately. The plan names ten targets,
  # which was the set the pages read when it was written. `factors.qmd` also
  # reads `factor_platform`, `methyl_platform` and `mutation_factor_annot`, and
  # `survival.qmd` reads `cox_fit`; without their edges the platform verdicts
  # and the EPV/split panel could go stale silently -- and those are precisely
  # the numbers that qualify the headline results.
  tar_target(
    dashboard_site,
    {
      invisible(list(
        mofa_factors, mofa_varexp, subtypes_df, km_subtype_df, concordance,
        survival_metrics, bap1_auroc, sanity_table, rna_mat, cox_fit,
        factor_platform, methyl_platform, mutation_factor_annot,
        gdc_live_panel
      ))
      quarto::quarto_render(input = "dashboard", quiet = TRUE)
      DASHBOARD_SITE_DIR
    },
    format = "file",
    packages = c(tar_option_get("packages"), "quarto")
  )
)
