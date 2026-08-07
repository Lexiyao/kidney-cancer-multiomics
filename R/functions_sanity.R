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

#' Events required to detect a hazard ratio, Schoenfeld's formula.
#'
#' d = (z_{1-alpha/2} + z_{power})^2 / (p * (1 - p) * log(HR)^2), with `p` the
#' fraction of the cohort in the exposed group. This is DESIGN arithmetic: it
#' answers "how many deaths would a study need to observe before an effect of
#' this size were detectable at this alpha", and it is deliberately implemented
#' as a function rather than left as a number in a comment, so the claim that
#' the BAP1 anchor is under-powered is EXECUTED and TESTED rather than asserted
#' in prose that can rot.
#'
#' It is NOT a post-hoc power calculation dressed up: nothing here interprets a
#' non-significant p-value as evidence of no effect. It bounds what the cohort
#' could ever have shown.
#'
#' Returns Inf for hr == 1 (a null effect is never detectable), which is the
#' mathematically correct limit and keeps callers from having to special-case it.
#'
#' @param hr hazard ratio to be detected (> 0, != 1 for a finite answer).
#' @param p_exposed fraction of the cohort carrying the exposure, in (0, 1).
#' @param alpha two-sided significance level.
#' @param power target power.
#' @return numeric(1) required number of EVENTS (not cases).
fn_schoenfeld_events <- function(hr, p_exposed,
                                 alpha = SANITY_MAX_P,
                                 power = SURVIVAL_TARGET_POWER) {
  stopifnot(is.numeric(hr), length(hr) == 1L, is.finite(hr), hr > 0)
  stopifnot(is.numeric(p_exposed), length(p_exposed) == 1L,
            is.finite(p_exposed), p_exposed > 0, p_exposed < 1)
  stopifnot(alpha > 0, alpha < 1, power > 0, power < 1)
  z_sum <- stats::qnorm(1 - alpha / 2) + stats::qnorm(power)
  z_sum^2 / (p_exposed * (1 - p_exposed) * log(hr)^2)
}

#' Check that BAP1-mutant ccRCC tumours have worse overall survival.
#'
#' DIRECTIONAL anchor (design spec section 7): the published finding is that
#' BAP1 loss is associated with WORSE survival, so `pass` requires HR > 1, not
#' merely "different from 1". A two-sided "is there any effect" test would be
#' satisfied by a protective BAP1 effect, which would contradict the literature
#' while still reporting green.
#'
#' NOTE: `pass` is directional only and carries NO significance requirement, and
#' on this cohort it CANNOT acquire one. The requirement was re-specified once,
#' on arithmetic, and the arithmetic is recorded in full in R/constants.R
#' (PUBLISHED_BAP1_HR_RANGE / SURVIVAL_TARGET_POWER block). In short: the
#' measured effect (HR 1.584 at 8.63% mutant prevalence) needs ~470 OS events
#' for 80% power at two-sided 0.05, and the 417-case mutation subset supplies
#' 138. A `p < 0.05` demand on this snapshot tested the size of TCGA KIRC, not
#' the correctness of this pipeline.
#'
#' CAVEAT ON THE 138: it is DERIVED, not recorded. No committed transcript
#' prints an event count for this fit — see the block in R/constants.R, which
#' also shows that the companion "~12 in the mutant arm" is in tension with the
#' recorded confidence interval (which implies ~18). The re-specification rests
#' on it, so it is labelled as derived until the next container run prints it;
#' `n_events` and friends are returned below precisely so it can be.
#'
#' What replaced it is STRICTER, not looser — but only because the SECOND thing
#' the old `ci_low > 1` guaranteed was replaced too. It carried two claims at
#' once: (a) the effect is significant, and (b) the fit is INFORMATIVE. (a) is
#' the mis-specified one; (b) is not, and dropping the line silently dropped it.
#' The anchor therefore now asserts, none of it a significance demand:
#'   - the point estimate lies inside PUBLISHED_BAP1_HR_RANGE, not merely > 1;
#'   - `mutant_frac` lies inside the SAME published BAP1 frequency band the
#'     mutation_freq check is anchored to, so the survival fit's exposed arm is
#'     the arm the rest of the suite certifies;
#'   - `n_events_mutant >= EPV_CAP`, a floor on the arm that carries the
#'     contrast (the existing `< n_events / 2` term is an upper bound only);
#'   - `ci_high / ci_low < BAP1_MAX_CI_RATIO`, a PRECISION bound — it is
#'     invariant to where the interval sits relative to 1, so it cannot turn a
#'     red check green;
#'   - the study is under-powered (`n_events < events_required`), so the
#'     limitation is itself a tested claim.
#' `p_value`, `ci_low`, `ci_high` and `n` are REPORTED and checked for
#' well-formedness, never for significance.
#'
#' REFUSES degenerate designs (see the guards in the body): fewer than
#' MIN_OS_EVENTS observed deaths, or a BAP1 column with only one level, both
#' stop() rather than returning an NA verdict.
#'
#' @param clinical data.frame with sample_id / os_time (days) / os_event (0/1).
#' @param mut_annot data.frame with sample_id and a 0/1 or logical BAP1 column.
#' @return list(label, hr, ci_low, ci_high, p_value, n, n_events, n_mutant,
#'   n_events_mutant, mutant_frac, events_required, underpowered, pass).
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
  # Design adequacy, reported alongside the estimate so no consumer can read the
  # p-value without also seeing what this cohort would have needed to produce a
  # significant one. See the arithmetic block in R/constants.R.
  n_mutant <- sum(merged$bap1_mut == 1L)
  n_events_mutant <- sum(merged$os_event == 1L & merged$bap1_mut == 1L)
  mutant_frac <- n_mutant / nrow(merged)
  events_required <- if (mutant_frac <= 0 || mutant_frac >= 1) {
    Inf
  } else {
    fn_schoenfeld_events(hr, mutant_frac)
  }
  list(
    label    = "BAP1-mutant tumours show worse overall survival (HR > 1)",
    hr       = hr,
    ci_low   = unname(s$conf.int["bap1_mut", "lower .95"]),
    ci_high  = unname(s$conf.int["bap1_mut", "upper .95"]),
    p_value  = p_value,
    n        = nrow(merged),
    n_events = n_events,
    n_mutant = n_mutant,
    n_events_mutant = n_events_mutant,
    mutant_frac = mutant_frac,
    # Events needed for SURVIVAL_TARGET_POWER at SANITY_MAX_P against an effect
    # the size of the one observed, at the observed exposure prevalence.
    events_required = events_required,
    underpowered    = n_events < events_required,
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

#' Re-run the k-means stratum search INSIDE each assay platform separately.
#'
#' The discriminating question when `platform_ari` fires is not "is the merged
#' partition contaminated" — that is already known — but "is there any m1-m4
#' structure UNDERNEATH the contamination". A real biological stratification
#' present in both arms would survive being looked for within each arm; a
#' partition that is only the assay split collapses, and the MERGED silhouette
#' then EXCEEDS both single-platform values, because the between-platform gap is
#' doing the separating.
#'
#' Each arm gets its OWN complete-CpG set (probes fail per platform), which is
#' why the per-arm CpG counts exceed the merged one on the real data.
#'
#' MEASURED on the frozen snapshot (GitHub Actions run 30911448546):
#'   HM27   n = 214  CpGs = 4658  silhouette 0.0858
#'   HM450  n = 310  CpGs = 4905  silhouette 0.0489
#'   merged n = 524  CpGs = 4568  silhouette 0.1197
#' The merged value is HIGHER than either arm alone. That ordering is the
#' signature of an assay-driven partition, and it is the evidence the m1-m4
#' verdict now carries with it.
#'
#' @param methyl_mat CpGs x samples M-value matrix (pre-completeness-filter, so
#'   each platform can be filtered on its own probes).
#' @param platform per-sample assay factor in the column order of methyl_mat.
#' @param k published number of strata.
#' @param seed fixed so each per-platform partition is reproducible.
#' @return data.frame(platform, n, n_cpg, silhouette); silhouette is NA for an
#'   arm too small or too incomplete to cluster (recorded, never imputed).
fn_within_platform_silhouette <- function(methyl_mat, platform,
                                          k = METHYL_N_STRATA,
                                          seed = SANITY_SEED) {
  stopifnot(is.matrix(methyl_mat),
            length(platform) == ncol(methyl_mat))
  restore <- fn_capture_rng()
  on.exit(restore(), add = TRUE)
  levs <- levels(factor(platform))
  rows <- lapply(levs, function(lv) {
    idx <- which(as.character(platform) == lv)
    sub <- methyl_mat[, idx, drop = FALSE]
    keep <- apply(is.finite(sub), 1L, all)
    n_cpg <- sum(keep)
    row <- data.frame(platform = lv, n = length(idx), n_cpg = n_cpg,
                      silhouette = NA_real_, stringsAsFactors = FALSE)
    # Too few samples for k groups, or nothing left to cluster on: report NA
    # rather than inventing a number. This is diagnostic reporting, not a
    # verdict term, so a missing arm must not stop the whole check.
    if (length(idx) <= k || n_cpg < 2L) {
      return(row)
    }
    set.seed(seed)
    feat <- t(sub[keep, , drop = FALSE])
    km <- stats::kmeans(feat, centers = k, nstart = 25L, iter.max = 100L)
    sil <- cluster::silhouette(km$cluster, stats::dist(feat))
    row$silhouette <- mean(sil[, "sil_width"])
    row
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
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
#' MEASURED OUTCOME ON THE REAL SNAPSHOT — this check is RED, and that is a
#' finding, not an unexplained failure. Run 30840373033: silhouette 0.1197,
#' Kruskal p 1.3e-82, but platform_ari 0.583 against the 0.25 ceiling
#' (platform_p 3.0e-113). Run 30911448546 then asked whether anything survives
#' within a single assay: HM27 0.0858 (n = 214), HM450 0.0489 (n = 310) — BOTH
#' BELOW the merged 0.1197. A merged silhouette that exceeds both single-platform
#' values is the signature of a partition separating assays rather than
#' biology, so `within_platform` and `merged_exceeds_within` are returned with
#' the verdict and `message` states the finding in words. The thresholds are
#' UNCHANGED; the failure is now informative instead of bare.
#'
#' @param k published number of strata; an ANCHOR, never tuned to the data.
#' @param max_platform_ari refusal ceiling on cluster-vs-platform agreement.
#' @param seed fixed so the k-means partition is reproducible.
#' @return list(label, message, n_strata, n_cpg_used, silhouette, kw_p_value,
#'   platform_ari, platform_p, within_platform, merged_exceeds_within, cluster,
#'   pass).
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
  within_platform <- NULL
  merged_exceeds_within <- NA
  if (!is.null(platform)) {
    platform_ari <- mclust::adjustedRandIndex(km$cluster, platform)
    platform_p <- suppressWarnings(
      stats::chisq.test(table(km$cluster, platform))$p.value
    )
    # The follow-up question the red verdict raises: is there m1-m4 structure
    # underneath the assay split? Computed on the RAW matrix so each arm is
    # filtered on its own probes (see fn_within_platform_silhouette).
    within_platform <- fn_within_platform_silhouette(methyl_mat, platform,
                                                     k = k, seed = seed)
    w <- within_platform$silhouette
    merged_exceeds_within <- if (all(is.na(w))) {
      NA
    } else {
      mean_sil > max(w, na.rm = TRUE)
    }
  }
  pass <- n_strata == k &&
    mean_sil > SANITY_MIN_SILHOUETTE &&
    kw$p.value < SANITY_MAX_P &&
    (is.na(platform_ari) || platform_ari < max_platform_ari) &&
    # A merged silhouette that BEATS both single-platform arms means the
    # between-assay gap is doing the separating. This conjunct introduces no new
    # threshold — it compares the check's own statistic to itself computed
    # within each arm — and it can only turn a green verdict RED. MEASURED on
    # the constructed controls: assay-only offsets 1.5/3.0 SD give merged 0.122
    # / 0.222 against within-arm 0.006 / 0.005 (TRUE, vetoed); four real strata
    # spread over both platforms give merged 0.293 against within-arm 0.294 /
    # 0.292 (FALSE, still green), with or without a 1.0 SD platform offset.
    !isTRUE(merged_exceeds_within)
  list(
    label        = "TCGA KIRC methylation resolves into m1-m4 strata",
    # State the verdict as a finding. A red m1-m4 light on this snapshot means
    # something specific and defensible, and a bare FALSE invites a future
    # reader to "fix" it by moving a threshold — which is exactly the failure
    # mode this suite exists to prevent.
    message      = fn_methyl_strata_message(pass, mean_sil, platform_ari,
                                            max_platform_ari, within_platform,
                                            merged_exceeds_within),
    n_strata     = n_strata,
    n_cpg_used   = nrow(complete),
    silhouette   = mean_sil,
    kw_p_value   = kw$p.value,
    platform_ari = platform_ari,
    platform_p   = platform_p,
    # Evidence that distinguishes "contaminated" from "contaminated AND empty".
    within_platform       = within_platform,
    merged_exceeds_within = merged_exceeds_within,
    cluster      = km$cluster,
    pass         = pass
  )
}

#' Put the m1-m4 verdict into words, including WHY it is red.
#'
#' Separated from fn_check_methyl_strata so the reporting logic stays small and
#' testable and the check body keeps to one job.
#'
#' @param pass the verdict.
#' @param mean_sil merged silhouette.
#' @param platform_ari cluster-vs-platform ARI (NA when no platform supplied).
#' @param max_platform_ari the ceiling it was judged against.
#' @param within_platform per-arm data.frame from fn_within_platform_silhouette.
#' @param merged_exceeds_within TRUE when the merged silhouette beats both arms.
#' @return character(1).
fn_methyl_strata_message <- function(pass, mean_sil, platform_ari,
                                     max_platform_ari, within_platform,
                                     merged_exceeds_within) {
  within_txt <- ""
  if (!is.null(within_platform)) {
    within_txt <- paste0(
      "; within-platform silhouettes ",
      paste(sprintf("%s %.4f (n=%d, %d CpGs)", within_platform$platform,
                    within_platform$silhouette, within_platform$n,
                    within_platform$n_cpg),
            collapse = ", "),
      if (isTRUE(merged_exceeds_within)) {
        " -- the MERGED silhouette EXCEEDS both, the signature of an assay-driven partition"
      } else if (isFALSE(merged_exceeds_within)) {
        paste0(" -- at least one arm matches or beats the merged value, so ",
               "structure survives within a platform")
      } else {
        ""
      }
    )
  }
  if (isTRUE(pass)) {
    return(sprintf("m1-m4 strata recovered (silhouette %.4f)%s",
                   mean_sil, within_txt))
  }
  if (!is.na(platform_ari) && platform_ari >= max_platform_ari) {
    fmt <- paste0(
      "FINDING, not an unexplained failure: the merged methylation ",
      "partition tracks ASSAY PLATFORM (ARI %.3f vs ceiling %.2f), so ",
      "the m1-m4 strata are NOT recovered on this snapshot ",
      "(merged silhouette %.4f)%s. The cohort is deliberately NOT ",
      "restricted to one platform and NOT batch-corrected (too few ",
      "cases are assayed on both platforms for a correction to be ",
      "validated, and the probe sets differ -- see the ",
      "methyl_platform_overlap target for the count); platform is ",
      "instead adjusted for downstream. This anchor stays RED as a ",
      "recorded negative result."
    )
    return(sprintf(fmt, platform_ari, max_platform_ari, mean_sil, within_txt))
  }
  # THE VETO-ONLY CASE, and it needs its own branch. When the ARI sits under the
  # ceiling but the merged silhouette still beats both single-platform arms,
  # `merged_exceeds_within` is the sole failing term — and falling through to the
  # floor message below stated a comparison the partition had CLEARED.
  # VERIFIED by direct call before this branch existed: platform_ari 0.10 with
  # merged silhouette 0.1197 (ABOVE the 0.10 floor) rendered as "m1-m4 strata
  # NOT recovered (silhouette 0.1197, floor 0.10)", i.e. it named the floor as
  # the cause. That is precisely the reading that invites a future maintainer to
  # "fix" the red by lowering SANITY_MIN_SILHOUETTE — the failure mode this
  # phase exists to prevent. Not reachable on the frozen snapshot (ARI 0.583
  # fires the branch above), but reachable the moment platform handling changes
  # downstream, which is the declared direction of Phase 4.
  if (isTRUE(merged_exceeds_within)) {
    fmt <- paste0(
      "FINDING, not an unexplained failure: the merged methylation partition ",
      "separates the ASSAYS rather than the tumours — its silhouette (%.4f) ",
      "EXCEEDS every within-platform arm, which is what an assay-driven ",
      "partition looks like%s. The silhouette floor (%.2f) is NOT the reason ",
      "for this red and lowering it would not address the finding; the ",
      "cluster-vs-platform ARI (%.3f) merely happens to sit under its %.2f ",
      "ceiling. This anchor stays RED as a recorded negative result."
    )
    return(sprintf(fmt, mean_sil, within_txt, SANITY_MIN_SILHOUETTE,
                   platform_ari, max_platform_ari))
  }
  if (mean_sil <= SANITY_MIN_SILHOUETTE) {
    return(sprintf("m1-m4 strata NOT recovered (silhouette %.4f, floor %.2f)%s",
                   mean_sil, SANITY_MIN_SILHOUETTE, within_txt))
  }
  # Everything else: the silhouette cleared its floor and the platform terms are
  # satisfied, so some OTHER conjunct (stratum count, Kruskal p) failed. Say so
  # rather than blaming a threshold that was met.
  sprintf(paste0("m1-m4 strata NOT recovered; the silhouette (%.4f) cleared its ",
                 "floor (%.2f) and the platform terms are satisfied, so the ",
                 "failing conjunct is the stratum count or the Kruskal test — ",
                 "read the full verdict object%s"),
          mean_sil, SANITY_MIN_SILHOUETTE, within_txt)
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

#' Check that the MOFA subtype assignment is INDEPENDENT of the assay platform.
#'
#' The one piece of good news in the platform diagnostic, pinned so a future
#' change cannot silently lose it.
#'
#' Individual MOFA factors carry a great deal of platform information on this
#' snapshot — Factor2 separates HM27 from HM450 at AUC 0.888, Factor5 at 0.818,
#' Factor6 at 0.735 (run 30911448546; that run's Factor6 p/q of 0 is a recorded
#' defect, see the annotation on the transcript — the AUC ordering is unaffected)
#' — so the raw material for a
#' platform-driven subtype assignment is unambiguously present. The 4-means
#' partition over the 15-factor space nevertheless does NOT recover it:
#' ARI(subtypes_mofa, platform) = 0.0058, cross-tab S1 11/9, S2 120/186,
#' S3 33/43, S4 50/72 over HM27/HM450. Under independence the expected ARI is 0.
#'
#' That is exactly the contrast with the m1-m4 methylation k-means, which IS
#' platform-driven (ARI 0.583) and whose anchor is red. Both facts are true at
#' once and both belong in the credibility suite.
#'
#' Direction of the guard: like SANITY_MAX_PLATFORM_ARI this term can only turn
#' a green verdict RED. Never raise `max_ari` to accommodate an observation.
#'
#' @param subtypes named factor/character of subtype labels, names = sample IDs
#'   (the shape `fn_assign_subtypes` returns).
#' @param platform named factor of assay platform, names = sample IDs.
#' @param max_ari refusal ceiling on subtype-vs-platform agreement.
#' @return list(label, ari, p_value, n, n_subtypes, cross_tab, pass).
fn_check_subtype_platform <- function(subtypes, platform,
                                      max_ari = SUBTYPE_MAX_PLATFORM_ARI) {
  stopifnot(length(subtypes) > 0L, length(platform) > 0L)
  if (is.null(names(subtypes)) || is.null(names(platform))) {
    stop("subtypes and platform must both be NAMED by sample id; unnamed ",
         "vectors would be compared in whatever order they happen to be in, ",
         "which can only produce a meaningless ARI")
  }
  # A partial overlap would score the guard on a silently reduced cohort, and a
  # near-zero ARI computed from a handful of cases is not evidence of anything.
  # Require the same sample set on both sides.
  if (!setequal(names(subtypes), names(platform))) {
    stop(sprintf(
      paste0("subtypes and platform cover different samples (%d vs %d, %d ",
             "shared): the platform-cleanliness guard must be scored on the ",
             "whole cohort, not on an accidental intersection"),
      length(subtypes), length(platform),
      length(intersect(names(subtypes), names(platform)))
    ))
  }
  ids <- names(subtypes)
  sub_v <- factor(as.character(subtypes[ids]))
  plat_v <- factor(as.character(platform[ids]))
  if (nlevels(sub_v) < 2L || nlevels(plat_v) < 2L) {
    stop("subtype-vs-platform independence needs at least two levels on both ",
         "sides; a constant vector gives ARI 0 for free and would report a ",
         "green light carrying no information")
  }
  ari <- mclust::adjustedRandIndex(sub_v, plat_v)
  cross_tab <- table(subtype = sub_v, platform = plat_v)
  p_value <- suppressWarnings(stats::chisq.test(cross_tab)$p.value)
  list(
    label      = "MOFA subtypes are independent of the HM27/HM450 assay platform",
    ari        = ari,
    p_value    = p_value,
    n          = length(ids),
    n_subtypes = nlevels(sub_v),
    cross_tab  = cross_tab,
    pass       = ari < max_ari
  )
}
