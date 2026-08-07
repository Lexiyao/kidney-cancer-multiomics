#' Attach the per-sample assay platform to the shared clinical frame.
#'
#' WHY THIS EXISTS. The HM27/HM450 split is real, large and NOT corrected
#' anywhere in this DAG (`fn_merge_methyl_platforms` only column-binds on the
#' common CpGs). MEASURED on the frozen curatedTCGAData KIRC snapshot, GitHub
#' Actions run 30911448546: the main cohort is 214 HM27 / 310 HM450; MOFA
#' Factor2 / Factor5 / Factor6 separate the two platforms at AUC 0.888 / 0.818 /
#' 0.735 (BH q = 2.9e-50 / 2.0e-34 / not reliably recorded — the transcript
#' printed p = q = 0 for Factor6, which is impossible at these group sizes and
#' which also depressed the other two q values; see the post-hoc annotation in
#' docs/results/platform-diagnosis-run-30911448546.txt); and the merged 4-means silhouette
#' (0.1197) is HIGHER than either arm computed alone (HM27 0.0858, HM450
#' 0.0489), i.e. a substantial part of the "structure" in the merged methylation
#' matrix IS the assay. ComBat and friends are deliberately NOT used: too few
#' cases are assayed on both platforms and the probe sets differ, so a batch
#' correction could not be validated here and might erase real signal along with
#' the batch. (The overlap has been stated as "3 cases" but is NOT in any
#' committed transcript; the `methyl_platform_overlap` target computes it, and
#' the next container run records it. The 214/310 split cannot corroborate it —
#' `methyl_platform` labels every dual-assayed case HM27, so the split is the
#' same for any overlap count.) The honest alternative is to ADJUST for platform
#' in the survival model — which first requires platform to travel with the
#' clinical frame.
#'
#' This is a plumbing function. It changes no threshold and no verdict: the
#' m1-m4 anchor stays red (platform_ari 0.583 vs a 0.25 ceiling), and nothing
#' here touches SANITY_MAX_PLATFORM_ARI or any published range.
#'
#' SCOPE AND THE DELIBERATE NAs. `clinical` covers every case in
#' `colData(mae_qc)` (536 on this snapshot), whereas `methyl_platform` is
#' defined on the 524-case main cohort (`common_ids`) only. Cases outside that
#' cohort therefore receive NA. The alternative — defaulting them to HM450, as
#' `methyl_platform`'s own `ifelse` does INSIDE the cohort, where every case is
#' known to carry methylation — would invent an assay for cases that may have no
#' methylation data at all. Module 4 restricts to `common_ids` before fitting,
#' so those NA rows never reach a model; the test suite asserts that every
#' cohort case is non-NA, which is the property the model actually depends on.
#'
#' Pure: returns a NEW data.frame, mutates neither argument, and preserves the
#' existing columns, their names, their order and their values (Module 4's
#' survival_df is built from sample_id / os_time / os_event).
#'
#' @param clinical data.frame with at least sample_id / os_time / os_event and
#'   unique sample_id values.
#' @param platform named factor with levels METHYL_PLATFORMS, one entry per
#'   cohort case, names being harmonised patient IDs (the `methyl_platform`
#'   target).
#' @return new data.frame: `clinical` plus a trailing `platform` factor column
#'   with levels METHYL_PLATFORMS, NA outside `names(platform)`.
fn_attach_platform <- function(clinical, platform) {
  if (!is.data.frame(clinical)) {
    stop("fn_attach_platform needs a clinical data.frame")
  }
  # The same trio fn_check_bap1_survival asserts; kept as literals here for the
  # same reason it does, rather than introducing a third name for one contract.
  absent <- setdiff(c("sample_id", "os_time", "os_event"), names(clinical))
  if (length(absent) > 0L) {
    stop("clinical lacks required survival columns: ",
         paste(absent, collapse = ", "))
  }
  if ("platform" %in% names(clinical)) {
    stop("clinical already carries a `platform` column; refusing to overwrite it")
  }
  if (anyDuplicated(clinical$sample_id) > 0L) {
    stop("clinical$sample_id is not unique; the platform join would fan out")
  }
  if (!is.factor(platform) || !identical(levels(platform), METHYL_PLATFORMS)) {
    stop("platform must be a factor with levels ",
         paste(METHYL_PLATFORMS, collapse = "/"))
  }
  if (anyNA(platform)) {
    stop("platform is NA for ", sum(is.na(platform)), " of ", length(platform),
         " cases it claims to cover; the assay split must be known for each")
  }
  ids <- names(platform)
  if (is.null(ids) || anyDuplicated(ids) > 0L || anyNA(ids)) {
    stop("platform must be named by unique, non-missing sample IDs")
  }
  uncovered <- setdiff(ids, clinical$sample_id)
  if (length(uncovered) > 0L) {
    stop(length(uncovered), " case(s) carry a platform but are absent from ",
         "clinical (e.g. ", paste(head(uncovered, 3L), collapse = ", "),
         "): the join would silently drop them")
  }
  attached <- factor(
    as.character(platform)[match(clinical$sample_id, ids)],
    levels = METHYL_PLATFORMS
  )
  out <- data.frame(clinical, stringsAsFactors = FALSE)
  out[["platform"]] <- attached
  out
}

#' Per-factor association between a MOFA factor and the assay platform.
#'
#' WHY THIS IS A DAG TARGET AND NOT A WORKFLOW SCRIPT. Phase 4 chooses its
#' survival predictors on platform-cleanliness (BH q > SANITY_MAX_P) and
#' variance explained, and NOTHING else. Until this function existed, that basis
#' lived only as hard-coded numbers in comments citing run 30911448546, computed
#' by a manually-triggered workflow (.github/workflows/diagnose-platform.yml)
#' that wrote a text artifact and touched no target. Nothing in the pipeline
#' could tell you that Factor1 and Factor4 are still clean.
#'
#' The asymmetry that motivated it: the SUBTYPE claim is falsifiable at run time
#' (fn_check_subtype_platform + SUBTYPE_MAX_PLATFORM_ARI + its anchor), while
#' the FACTOR claim — the actual selection criterion — was not. If MOFA output
#' shifts (seed, MOFA2 version, upstream matrix change) Factor1/Factor4 could
#' become platform-loaded and the pipeline would report nothing; symmetrically,
#' re-adding "Factor2" to the predictor vector would break no check.
#'
#' The statistic is the one the diagnostic used, so the target is comparable
#' with the recorded transcript: a two-sample Wilcoxon of the factor scores
#' across the two assays, with the rank-biserial AUC as the readable effect size
#' (0.5 = the factor carries no platform information, 1.0 = it separates the
#' assays perfectly). AUC is reported as max(auc, 1 - auc) so it does not depend
#' on which platform happens to be the first level. Full precision is kept here;
#' the transcript rounds for display.
#'
#' DIRECTION OF THE GUARD: the anchor built on this can only turn a green
#' verdict RED. It asserts that the factors already wired as predictors are
#' platform-clean; it never selects them, and it never looks at survival.
#'
#' @param factors numeric matrix, samples x factors, rownames = sample IDs.
#' @param platform named factor of assay platform, names = sample IDs. Must
#'   cover every row of `factors`.
#' @return data.frame(factor, n_hm27-style per-level counts collapsed into `n`,
#'   auc, p_value, q_value), one row per factor, in `colnames(factors)` order.
#'   `q_value` is BH across the factors tested.
fn_factor_platform_association <- function(factors, platform) {
  if (!is.matrix(factors) || !is.numeric(factors)) {
    stop("factors must be a numeric samples x factors matrix")
  }
  ids <- rownames(factors)
  if (is.null(ids) || anyDuplicated(ids) > 0L) {
    stop("factors must have unique rownames (sample IDs); an unnamed matrix ",
         "would be compared against platform in arbitrary order")
  }
  if (is.null(names(platform))) {
    stop("platform must be NAMED by sample id")
  }
  absent <- setdiff(ids, names(platform))
  if (length(absent) > 0L) {
    stop(length(absent), " sample(s) in `factors` carry no platform label ",
         "(e.g. ", paste(utils::head(absent, 3L), collapse = ", "), ")")
  }
  p <- factor(as.character(platform[ids]))
  if (nlevels(p) != 2L) {
    stop("the platform comparison needs exactly two assay levels, found ",
         nlevels(p))
  }
  levs <- levels(p)
  rows <- lapply(colnames(factors), function(f) {
    x <- factors[, f]
    if (!all(is.finite(x))) {
      stop("factor ", f, " carries non-finite scores; the platform comparison ",
           "would silently drop samples and report an AUC over a denominator ",
           "that never saw them")
    }
    a <- x[p == levs[1L]]
    b <- x[p == levs[2L]]
    w <- suppressWarnings(stats::wilcox.test(a, b))
    raw_auc <- as.numeric(w$statistic) / (length(a) * length(b))
    data.frame(factor = f, n = length(x),
               auc = max(raw_auc, 1 - raw_auc),
               p_value = w$p.value, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$q_value <- stats::p.adjust(out$p_value, "BH")
  # A p-value of EXACTLY zero is not a very small p-value, it is a broken one.
  # Above n = 50 per group wilcox.test uses the normal approximation, whose
  # p-value has a strictly positive floor — MEASURED at the real group sizes
  # (214 vs 310) it is 1.98e-84 even at perfect separation. Zero therefore means
  # the approximation's variance collapsed (heavy ties send SIGMA to 0 and z to
  # -Inf), which also makes that row's AUC meaningless.
  #
  # This is not hypothetical: run 30911448546 recorded Factor6 with AUC 0.735
  # and p = q = 0, ~20 orders of magnitude from what that AUC implies, and
  # because BH ranks a zero first it depressed the reported Factor2 and Factor5
  # q values by factors of 2 and 5/3. See the post-hoc annotation in
  # docs/results/platform-diagnosis-run-30911448546.txt.
  #
  # FLAGGED rather than fatal, for the same reason fn_within_platform_silhouette
  # reports NA for an unclusterable arm: this is a diagnostic target and one bad
  # factor must not take out the DAG. The anchor asserts the WIRED factors are
  # not degenerate, which is the property predictor selection depends on.
  out$degenerate <- !(out$p_value > 0)
  rownames(out) <- NULL
  out
}
