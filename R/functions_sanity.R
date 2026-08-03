# Module 3: literature positive controls as structured pass/fail objects.
# Each fn_check_* is pure (returns a new list, never mutates its inputs).
#
# These are the project's CREDIBILITY ANCHOR (design spec section 7). Every
# check compares an OBSERVED quantity against an INDEPENDENT published anchor
# held in R/constants.R, and returns a structured object with a `pass` flag —
# never just a plot or a printed number. A control that cannot fail proves
# nothing, so the published ranges, thresholds and directions here must never
# be tuned to make an observed value pass.

#' Capture the global RNG stream and return a function that restores it.
#'
#' The two seeded checks below call `set.seed()` so their k-means partitions are
#' reproducible, but `set.seed()` mutates the GLOBAL stream. Left unrestored,
#' the session sits at the seed-42 state once `sanity_results` has built, and any
#' later randomised code in the same process — other targets, or subsequent
#' `test_that` blocks that draw randoms without re-seeding — is silently coupled
#' to Module 3 and to whether Module 3 ran at all. That also makes test outcomes
#' order-dependent. MEASURED before this was added: `.Random.seed` differed
#' before vs after each call, and after `set.seed(7)` a caller's `rnorm(3)[1]`
#' was 2.287 without the check but 0.581 with it.
#'
#' Determinism of the checks themselves is unaffected: `set.seed()` still runs
#' before `kmeans`; only the caller's stream is put back afterwards. This is what
#' the file header's "each fn_check_* is pure" claim requires.
#'
#' Usage: `restore <- fn_capture_rng(); on.exit(restore(), add = TRUE)` — the
#' `on.exit` must be registered by the caller so it fires on that frame's exit.
#'
#' @return a zero-argument function that restores the captured stream.
fn_capture_rng <- function() {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  function() {
    if (is.null(old_seed)) {
      suppressWarnings(rm(".Random.seed", envir = .GlobalEnv))
    } else {
      # `.Random.seed` is base R's reserved name, not ours to rename (same
      # nolint as R/functions_integrate.R:130).
      assign(".Random.seed", old_seed, envir = .GlobalEnv) # nolint: object_name_linter.
    }
    invisible(NULL)
  }
}

#' Check ccRCC driver mutation frequencies against published ranges.
#' @param mut_annot data.frame with sample_id + one 0/1 or logical column per
#'   gene. Module 1's `fn_extract_mutation_status` emits INTEGER 0/1 columns
#'   (plus rownames == sample_id); the plan's interface note says "logical".
#'   `as.logical()` below accepts both, and both shapes are unit-tested.
#' @param gene_panel driver genes eligible to be scored. The scored set is the
#'   three-way intersection of this panel, `names(ranges)` and the columns
#'   actually present. Under the default this changes nothing — DRIVER_GENE_PANEL
#'   is a superset of `names(PUBLISHED_MUT_FREQ_RANGES)` — but a caller passing
#'   `gene_panel = "BAP1"` now gets a verdict about BAP1 alone. (It previously
#'   was NOT read at all and silently scored all four ranged genes.)
#' @param ranges named list of c(low=, high=) published frequency bounds.
#'   INDEPENDENT literature anchor — never widen it to make an observation pass.
#' @return list(label, per_gene, n, pass).
fn_check_mutation_freq <- function(mut_annot,
                                   gene_panel = DRIVER_GENE_PANEL,
                                   ranges = PUBLISHED_MUT_FREQ_RANGES) {
  stopifnot(is.data.frame(mut_annot), "sample_id" %in% names(mut_annot))
  genes <- intersect(intersect(names(ranges), gene_panel), colnames(mut_annot))
  if (length(genes) == 0L) {
    stop("mut_annot contains none of the ranged driver genes")
  }
  rows <- lapply(genes, function(g) {
    # VALIDATE THE COERCION. as.logical() on a character or factor 0/1 column
    # returns all NA, and na.rm = TRUE then averages an empty vector.
    # MEASURED with a character VHL column: observed = NaN, that gene's pass =
    # NA, and the OVERALL verdict = NA — a frequency anchor that scored nothing
    # and reported no failure.
    col <- mut_annot[[g]]
    if (!(is.logical(col) || is.numeric(col))) {
      stop(sprintf("mut_annot$%s must be logical or 0/1 numeric, got %s",
                   g, class(col)[1]))
    }
    observed <- mean(as.logical(col), na.rm = TRUE)
    if (!is.finite(observed)) {
      stop(sprintf("observed frequency for %s is not finite (no scorable values)",
                   g))
    }
    rng <- ranges[[g]]
    data.frame(
      gene     = g,
      observed = observed,
      low      = unname(rng["low"]),
      high     = unname(rng["high"]),
      pass     = observed >= rng["low"] & observed <= rng["high"],
      stringsAsFactors = FALSE
    )
  })
  per_gene <- do.call(rbind, rows)
  rownames(per_gene) <- NULL
  # `all()` over an NA turns the overall verdict into NA, which is not a verdict.
  stopifnot(!anyNA(per_gene$pass))
  list(
    label    = "ccRCC driver mutation frequencies within published ranges",
    per_gene = per_gene,
    # The DENOMINATOR the proportions came from. Without it, `mean()` over 20
    # rows and over 417 is indistinguishable in the output, so a regression that
    # collapsed the cohort — while leaving the frequencies inside the published
    # bands, which small samples easily do — would satisfy every anchor. (Per
    # gene the effective denominator can be smaller if that column carries NAs;
    # this is the cohort size, which is what the anchor bounds.)
    n        = nrow(mut_annot),
    pass     = all(per_gene$pass)
  )
}

#' Check that BAP1-mutant ccRCC tumours have worse overall survival.
#'
#' DIRECTIONAL anchor (design spec section 7): the published finding is that
#' BAP1 loss is associated with WORSE survival, so `pass` requires HR > 1, not
#' merely "different from 1". A two-sided "is there any effect" test would be
#' satisfied by a protective BAP1 effect, which would contradict the literature
#' while still reporting green.
#'
#' NOTE: `pass` is directional only and carries NO significance requirement.
#' Consumers that need a significant result must read `p_value` alongside it
#' (the Task 3.7 credibility anchor asserts both).
#'
#' REFUSES degenerate designs (see the guards in the body): fewer than
#' MIN_OS_EVENTS observed deaths, or a BAP1 column with only one level, both
#' stop() rather than returning an NA verdict.
#'
#' @param clinical data.frame with sample_id / os_time (days) / os_event (0/1).
#' @param mut_annot data.frame with sample_id and a 0/1 or logical BAP1 column.
#' @return list(label, hr, ci_low, ci_high, p_value, n, n_events, pass).
fn_check_bap1_survival <- function(clinical, mut_annot) {
  stopifnot(all(c("sample_id", "os_time", "os_event") %in% names(clinical)))
  stopifnot("BAP1" %in% names(mut_annot))
  merged <- merge(
    clinical[, c("sample_id", "os_time", "os_event")],
    mut_annot[, c("sample_id", "BAP1")],
    by = "sample_id"
  )
  merged <- merged[stats::complete.cases(merged) & merged$os_time > 0, ]
  merged$bap1_mut <- as.integer(as.logical(merged$BAP1))
  # REFUSE rather than report NA. survival::coxph does not error on a degenerate
  # design: with zero events, or with BAP1 constant, it returns quietly with NA
  # coefficients and this function emitted `pass = NA` — not a verdict, and a
  # violation of the suite's own shape anchor. With a single event it returned
  # hr = 3.7e+09, ci = [0, Inf], p = 0.999 and pass = TRUE, i.e. GREEN from a fit
  # carrying no information. Both are checked before anything is estimated.
  stopifnot(all(merged$os_event %in% c(0L, 1L)))
  n_events <- sum(merged$os_event == 1L)
  if (n_events < MIN_OS_EVENTS) {
    stop(sprintf(
      paste0("BAP1 survival control has only %d OS events among %d cases ",
             "(floor %d): check the vital_status decode ",
             "(VITAL_STATUS_DEAD_VALUES). A Cox model fitted on too few ",
             "events cannot fail and proves nothing."),
      n_events, nrow(merged), MIN_OS_EVENTS
    ))
  }
  if (length(unique(merged$bap1_mut)) < 2L) {
    stop("BAP1 is constant across the merged cohort; the control cannot be ",
         "fitted (no wild-type/mutant contrast exists)")
  }
  fit <- survival::coxph(
    survival::Surv(os_time, os_event) ~ bap1_mut,
    data = merged
  )
  s <- summary(fit)
  hr <- unname(s$conf.int["bap1_mut", "exp(coef)"])
  p_value <- unname(s$coefficients["bap1_mut", "Pr(>|z|)"])
  if (!is.finite(hr) || !is.finite(p_value)) {
    stop("coxph returned a degenerate BAP1 fit (non-finite HR or p-value)")
  }
  list(
    label    = "BAP1-mutant tumours show worse overall survival (HR > 1)",
    hr       = hr,
    ci_low   = unname(s$conf.int["bap1_mut", "lower .95"]),
    ci_high  = unname(s$conf.int["bap1_mut", "upper .95"]),
    p_value  = p_value,
    n        = nrow(merged),
    n_events = n_events,
    pass     = hr > 1
  )
}

#' Keep only CpGs that are finite in every sample.
#'
#' `stats::kmeans` errors outright on NA/NaN/Inf, and the real Module 1
#' `methyl_mat` is 91.4% complete (432/5000 CpGs carry a non-finite value from
#' failed probes), so without this the m1-m4 anchor cannot run on the data it
#' exists to check. Identical failure mode and identical resolution to
#' `fn_complete_features` for SNF; imputation is deliberately NOT used, because
#' invented values in a positive control would be invented evidence.
#'
#' Dropping is per-CpG, never per-sample: every case stays in the stratum
#' assignment. The floor is the surviving FRACTION, so a matrix gutted by
#' missingness (e.g. one bad sample turning most CpGs incomplete) stops loudly
#' instead of clustering its remnants into a silent green.
#'
#' @param methyl_mat CpGs x samples numeric matrix.
#' @param min_frac minimum fraction of CpGs that must be complete.
#' @return new matrix restricted to complete CpGs (inputs unchanged).
fn_complete_cpgs <- function(methyl_mat, min_frac = SANITY_MIN_COMPLETE_FRAC) {
  keep <- apply(is.finite(methyl_mat), 1L, all)
  n_keep <- sum(keep)
  if (n_keep < 2L || n_keep < min_frac * nrow(methyl_mat)) {
    stop(sprintf(
      "methyl_mat retains %d of %d complete CpGs (%.1f%%, floor %.0f%%)",
      n_keep, nrow(methyl_mat), 100 * n_keep / nrow(methyl_mat), 100 * min_frac
    ))
  }
  methyl_mat[keep, , drop = FALSE]
}

#' Recover the four TCGA KIRC DNA-methylation strata (m1-m4).
#'
#' FALSIFIABILITY, measured on constructed data, so the limits are on record:
#' only the SILHOUETTE term discriminates. On 1000 iid-noise CpGs x 80 samples
#' (no strata whatsoever) this returns silhouette 0.005 -> pass FALSE, while
#' `n_strata == k` is TAUTOLOGICAL (Hartigan-Wong k-means with centers = 4
#' errors rather than returning fewer non-empty groups) and the Kruskal test is
#' CIRCULAR — it compares column means across clusters derived from those same
#' columns, and fires at p = 0.009 on that pure noise. Both weak conjuncts are
#' kept because the plan specifies them and neither can mask a low silhouette,
#' but neither is evidence on its own. See the negative-control test.
#'
#' PLATFORM CONFOUND, made explicit and falsifiable: `methyl_mat` is
#' cbind(HM27, HM450) with NO batch correction, and platform is the strongest
#' single axis in merged 27k/450k M-values. MEASURED on constructed data with
#' iid noise, NO biological strata and only a per-platform mean offset, this
#' check returned pass = TRUE from an offset of 1.5 SD upward (silhouette 0.122
#' / 0.164 / 0.222 at 1.5 / 2.0 / 3.0 SD, kruskal p = 2.7e-80), i.e. a green
#' light produced entirely by the assay. `platform_ari` is therefore part of
#' `pass`: a partition that merely reproduces the assay cannot report green.
#'
#' @param methyl_mat CpGs x samples M-value matrix (top-variable).
#' @param platform per-sample assay factor, in the column order of methyl_mat
#'   (named, it is checked against colnames). NULL is only legitimate for
#'   single-platform constructed data: the platform term is then reported as NA
#'   and dropped from `pass`, and the Task 3.7 anchor REQUIRES a finite
#'   platform_ari, so the real run cannot silently lose this guard.
#' @param k published number of strata; an ANCHOR, never tuned to the data.
#' @param max_platform_ari refusal ceiling on cluster-vs-platform agreement.
#' @param seed fixed so the k-means partition is reproducible.
#' @return list(label, n_strata, n_cpg_used, silhouette, kw_p_value,
#'   platform_ari, platform_p, cluster, pass).
fn_check_methyl_strata <- function(methyl_mat, platform = NULL,
                                   k = METHYL_N_STRATA,
                                   max_platform_ari = SANITY_MAX_PLATFORM_ARI,
                                   seed = SANITY_SEED) {
  stopifnot(is.matrix(methyl_mat), ncol(methyl_mat) > k)
  if (!is.null(platform)) {
    if (length(platform) != ncol(methyl_mat)) {
      stop(sprintf("platform has %d entries but methyl_mat has %d samples",
                   length(platform), ncol(methyl_mat)))
    }
    if (!is.null(names(platform)) &&
          !identical(names(platform), colnames(methyl_mat))) {
      stop("platform names do not match colnames(methyl_mat): the assay ",
           "labels would be misaligned with the samples they describe")
    }
    if (length(unique(platform)) < 2L) {
      stop("platform has a single level; pass NULL for single-platform data")
    }
  }
  complete <- fn_complete_cpgs(methyl_mat)
  restore <- fn_capture_rng()
  on.exit(restore(), add = TRUE)
  set.seed(seed)
  feat <- t(complete)                        # samples x CpGs
  km <- stats::kmeans(feat, centers = k, nstart = 25L, iter.max = 100L)
  d <- stats::dist(feat)
  sil <- cluster::silhouette(km$cluster, d)
  mean_sil <- mean(sil[, "sil_width"])
  # Global methylation differs across strata (CIMP-like stratum) -> Kruskal.
  # Taken over the SAME complete CpGs the clusters came from, so the "global"
  # mean is not itself a function of which probes failed in which sample.
  sample_mean_m <- colMeans(complete, na.rm = TRUE)
  kw <- stats::kruskal.test(sample_mean_m, factor(km$cluster))
  n_strata <- length(unique(km$cluster))
  # Is the partition just the assay? See the roxygen note and the calibration in
  # R/constants.R. Reported as NA when no platform vector was supplied, which
  # the Task 3.7 anchor refuses.
  platform_ari <- NA_real_
  platform_p <- NA_real_
  if (!is.null(platform)) {
    platform_ari <- mclust::adjustedRandIndex(km$cluster, platform)
    platform_p <- suppressWarnings(
      stats::chisq.test(table(km$cluster, platform))$p.value
    )
  }
  list(
    label        = "TCGA KIRC methylation resolves into m1-m4 strata",
    n_strata     = n_strata,
    n_cpg_used   = nrow(complete),
    silhouette   = mean_sil,
    kw_p_value   = kw$p.value,
    platform_ari = platform_ari,
    platform_p   = platform_p,
    cluster      = km$cluster,
    pass         = n_strata == k &&
      mean_sil > SANITY_MIN_SILHOUETTE &&
      kw$p.value < SANITY_MAX_P &&
      (is.na(platform_ari) || platform_ari < max_platform_ari)
  )
}

#' Check ccA/ccB expression-signature separation (Brannon 2010 proxy).
#' ccB = proliferative axis-high, ccA = angiogenic axis-low.
#'
#' DELIBERATE STRENGTHENING of the plan's body, which is otherwise kept
#' verbatim. As specified, `pass` required only (a) mean silhouette >
#' SANITY_MIN_SILHOUETTE and (b) a Wilcoxon p on the ccB-ccA axis across the
#' k-means groups. BOTH are circular: the groups are k-means on the very pair
#' of scores that then get clustered and compared, so 2-D k-means always splits
#' the cloud (silhouette 0.37-0.44 even on pure noise) and the Wilcoxon can
#' fire on an arbitrary split direction. MEASURED on five structureless
#' matrices carrying every published marker and no ccA/ccB structure at all,
#' the specified pair returned pass = TRUE on three of five. A positive control
#' that greenlights noise is worse than none.
#'
#' The added conjunct is the published claim itself and is NON-CIRCULAR: the
#' proliferation and angiogenesis programmes must be OPPOSED across tumours
#' (Spearman rho < 0, p < SANITY_MAX_P), computed from the marker scores alone
#' with no reference to the clustering. It uses no new threshold — the
#' direction is the null (0) and the alpha is the existing SANITY_MAX_P — and
#' it makes the check STRICTER, never more permissive, so it cannot turn a real
#' failure green. Measured behaviour: rejects all five structureless matrices,
#' rejects markers that all move together (rho ~ +0.99), accepts genuine
#' opposition (rho = -0.71, p = 9.9e-07).
#'
#' The panels are simplified proxies, NOT ClearCode34; that limitation is
#' stated in R/constants.R and belongs in the dashboard too.
#'
#' @param rna_mat genes x samples log2-normalised expression, symbol rownames.
#' @param ccb_markers,cca_markers published marker panels — ANCHORS, never
#'   edited to make an observation pass.
#' @param min_markers refusal floor on surviving markers PER PANEL. Not an
#'   anchor: it touches no threshold and no direction, so it cannot turn a
#'   failing check green — it only stops a gutted panel scoring a published axis
#'   from a couple of genes.
#' @param seed fixed so the k-means partition is reproducible.
#' @return list(label, n_ccb_used, n_cca_used, markers_used, silhouette,
#'   separation_p_value, anticorr_rho, anticorr_p_value, group, axis_score,
#'   pass).
fn_check_ccab_signature <- function(rna_mat,
                                    ccb_markers = CCB_PROLIFERATION_MARKERS,
                                    cca_markers = CCA_ANGIOGENESIS_MARKERS,
                                    min_markers = SANITY_MIN_MARKERS_PER_PANEL,
                                    seed = SANITY_SEED) {
  stopifnot(is.matrix(rna_mat), !is.null(rownames(rna_mat)))
  ccb <- intersect(ccb_markers, rownames(rna_mat))
  cca <- intersect(cca_markers, rownames(rna_mat))
  if (length(ccb) < min_markers || length(cca) < min_markers) {
    stop(sprintf(
      paste0("insufficient ccA/ccB marker genes present in rna_mat: ",
             "%d/%d ccB and %d/%d ccA survived (floor %d per panel)"),
      length(ccb), length(ccb_markers), length(cca), length(cca_markers),
      min_markers
    ))
  }
  ccb_score <- colMeans(rna_mat[ccb, , drop = FALSE])
  cca_score <- colMeans(rna_mat[cca, , drop = FALSE])
  # A high axis value marks a ccB-like (proliferative) tumour.
  axis <- ccb_score - cca_score
  restore <- fn_capture_rng()
  on.exit(restore(), add = TRUE)
  set.seed(seed)
  # Cluster AND measure in the SAME space. This previously clustered on the
  # standardised scores but took the silhouette from a raw-space distance, so
  # the reported number was not the separation quality of the partition being
  # reported. MEASURED on a matrix whose two panels differ ~37x in spread:
  # reported silhouette -0.020 versus 0.538 in the space k-means actually used.
  # fn_check_methyl_strata already does this correctly.
  feat <- cbind(scale(ccb_score), scale(cca_score))
  km <- stats::kmeans(feat, centers = 2L, nstart = 25L)
  sil <- cluster::silhouette(km$cluster, stats::dist(feat))
  wt <- stats::wilcox.test(axis ~ factor(km$cluster))
  mean_sil <- mean(sil[, "sil_width"])
  # The non-circular term: proliferation vs angiogenesis must be OPPOSED.
  # exact = FALSE pins the normal approximation so tied expression values
  # cannot change the test used (and cannot raise a warning) run to run.
  rho_test <- stats::cor.test(ccb_score, cca_score, method = "spearman",
                              exact = FALSE)
  rho <- unname(rho_test$estimate)
  list(
    label              = "ccA/ccB expression signatures separate into two groups",
    # The panel ACTUALLY used. Without this the returned object was identical
    # whether the full 6-vs-6 published panel scored the axis or a silently
    # truncated 2-vs-2 remnant did, so no anchor could detect attrition.
    n_ccb_used         = length(ccb),
    n_cca_used         = length(cca),
    markers_used       = list(ccb = ccb, cca = cca),
    silhouette         = mean_sil,
    separation_p_value = wt$p.value,
    anticorr_rho       = rho,
    anticorr_p_value   = rho_test$p.value,
    group              = km$cluster,
    axis_score         = axis,
    # SANITY_MIN_SILHOUETTE_2D, not SANITY_MIN_SILHOUETTE: this is a 2-D 2-means
    # whose null distribution is nothing like the methylation check's. See the
    # calibration recorded in R/constants.R (null ceiling 0.466 over 600 runs).
    pass               = mean_sil > SANITY_MIN_SILHOUETTE_2D &&
      wt$p.value < SANITY_MAX_P &&
      rho < 0 &&
      rho_test$p.value < SANITY_MAX_P
  )
}
