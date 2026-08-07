# Module 3 (sanity): the credibility anchor (spec section 7). Literature
# positive controls as real assertions, not figures. Task 3.1 covers the
# published constants themselves; Tasks 3.2-3.5 append the fn_check_* unit
# tests and Task 3.7 appends the anchors on the frozen sanity_results target.

test_that("sanity constants are defined with the published ccRCC ranges", {
  # Arrange / Act — constants sourced via helper-fixtures.R

  # Assert
  expect_setequal(DRIVER_GENE_PANEL,
                  c("VHL", "PBRM1", "SETD2", "BAP1", "MTOR", "KDM5C"))
  expect_true(all(c("VHL", "PBRM1", "SETD2", "BAP1") %in%
                    names(PUBLISHED_MUT_FREQ_RANGES)))
  expect_equal(unname(PUBLISHED_MUT_FREQ_RANGES$VHL), c(0.40, 0.60))
  expect_true(PUBLISHED_MUT_FREQ_RANGES$BAP1["low"] < 0.18)
  expect_identical(METHYL_N_STRATA, 4L)
  expect_true(SANITY_MAX_P == 0.05)
  expect_length(CCB_PROLIFERATION_MARKERS, 6L)
  expect_length(CCA_ANGIOGENESIS_MARKERS, 6L)

  # --- Additive strengthening (beyond the plan's Task 3.1 block) -------------
  # The assertions above still pass if a published range is silently widened to
  # c(low = 0, high = 1) — which would make fn_check_mutation_freq unfalsifiable
  # and is exactly the drift this suite exists to catch. Pin every anchor and
  # enforce the structural invariants the Task 3.2 check relies on.
  expect_identical(
    PUBLISHED_MUT_FREQ_RANGES,
    list(
      VHL   = c(low = 0.40, high = 0.60),
      PBRM1 = c(low = 0.28, high = 0.45),
      SETD2 = c(low = 0.08, high = 0.18),
      BAP1  = c(low = 0.06, high = 0.18)
    )
  )
  for (rng in PUBLISHED_MUT_FREQ_RANGES) {
    expect_named(rng, c("low", "high"))
    expect_true(rng[["low"]] < rng[["high"]])
    expect_true(rng[["low"]] >= 0 && rng[["high"]] <= 1)
  }

  # Collision resolution: DRIVER_GENES is the canonical panel (Modules 1-2 and
  # _targets.R consume it); DRIVER_GENE_PANEL is a single-source alias, so the
  # two literals can never drift apart. MUTATION_FREQ_RANGES was the duplicate
  # of PUBLISHED_MUT_FREQ_RANGES and must stay deleted.
  expect_identical(DRIVER_GENE_PANEL, DRIVER_GENES)
  expect_false(exists("MUTATION_FREQ_RANGES", inherits = TRUE))

  # Marker panels must be disjoint, non-empty gene symbols: an overlapping
  # panel would make the ccA/ccB axis score partly self-cancelling.
  expect_length(intersect(CCB_PROLIFERATION_MARKERS,
                          CCA_ANGIOGENESIS_MARKERS), 0L)
  expect_true(all(nzchar(c(CCB_PROLIFERATION_MARKERS,
                           CCA_ANGIOGENESIS_MARKERS))))

  expect_identical(SANITY_MIN_SILHOUETTE, 0.10)
  expect_identical(SANITY_SEED, 42L)
})

test_that("the platform-confound thresholds are PINNED at their calibrated values", {
  # These are the constants a future reader would be tempted to move in order to
  # turn the m1-m4 anchor green. Pin them, so doing so is a deliberate, visible
  # edit to this test rather than a one-character change in constants.R.
  #
  # MEASURED against them on the frozen snapshot: methylation k-means
  # platform_ari 0.583 (run 30840373033) — RED, and correctly so; MOFA subtypes
  # vs platform ARI 0.0058 (run 30911448546) — clean.
  expect_identical(SANITY_MAX_PLATFORM_ARI, 0.25)
  expect_identical(SUBTYPE_MAX_PLATFORM_ARI, 0.05)
  # The subtype guard is the tighter of the two: an ARI expected to be 0 under
  # independence has no business being judged against a ceiling calibrated for a
  # k-means that might legitimately correlate with anything.
  expect_lt(SUBTYPE_MAX_PLATFORM_ARI, SANITY_MAX_PLATFORM_ARI)
  # Both are strictly inside (0, 1): an ARI ceiling of 1 could never fire.
  expect_true(SANITY_MAX_PLATFORM_ARI > 0 && SANITY_MAX_PLATFORM_ARI < 1)
  expect_true(SUBTYPE_MAX_PLATFORM_ARI > 0 && SUBTYPE_MAX_PLATFORM_ARI < 1)
  expect_identical(METHYL_PLATFORMS, c("HM27", "HM450"))
})

test_that("the BAP1 effect band is a published anchor, not a fitted one", {
  # PUBLISHED_BAP1_HR_RANGE replaced the `p < 0.05` / `ci_low > 1` demand the
  # cohort cannot meet (see the BAP1 anchor for the Schoenfeld arithmetic). It
  # is only a legitimate replacement if it is STRICTER than the `hr > 1` it sits
  # beside and if it excludes the values a broken fit produces — so assert that,
  # not merely that the constant exists.
  expect_named(PUBLISHED_BAP1_HR_RANGE, c("low", "high"))
  expect_identical(unname(PUBLISHED_BAP1_HR_RANGE), c(1.2, 3.0))
  # Strictly inside the harmful direction, so the band can never admit a
  # protective or null effect.
  expect_gt(PUBLISHED_BAP1_HR_RANGE[["low"]], 1)
  expect_lt(PUBLISHED_BAP1_HR_RANGE[["low"]], PUBLISHED_BAP1_HR_RANGE[["high"]])
  # The MEASURED point estimate must sit inside it — if a future snapshot moves
  # outside, that is a finding to report, not a bound to widen.
  expect_gte(1.583991, PUBLISHED_BAP1_HR_RANGE[["low"]])
  expect_lte(1.583991, PUBLISHED_BAP1_HR_RANGE[["high"]])
  # The failure modes the old significance demand was standing in for.
  expect_lt(1.0001, PUBLISHED_BAP1_HR_RANGE[["low"]])   # a null "pass"
  expect_gt(40, PUBLISHED_BAP1_HR_RANGE[["high"]])      # a merge blow-up

  expect_identical(SURVIVAL_TARGET_POWER, 0.80)
  expect_true(SURVIVAL_TARGET_POWER > 0 && SURVIVAL_TARGET_POWER < 1)
})

test_that("fn_check_mutation_freq passes when all driver freqs are in range", {
  # Arrange — VHL 0.50, PBRM1 0.35, SETD2 0.12, BAP1 0.10
  n <- 100L
  mut_annot <- data.frame(
    sample_id = paste0("P", seq_len(n)),
    VHL   = c(rep(TRUE, 50), rep(FALSE, 50)),
    PBRM1 = c(rep(TRUE, 35), rep(FALSE, 65)),
    SETD2 = c(rep(TRUE, 12), rep(FALSE, 88)),
    BAP1  = c(rep(TRUE, 10), rep(FALSE, 90)),
    stringsAsFactors = FALSE
  )

  # Act
  res <- fn_check_mutation_freq(mut_annot)

  # Assert
  expect_true(res$pass)
  expect_equal(res$per_gene$observed[res$per_gene$gene == "VHL"], 0.50)
  expect_true(all(res$per_gene$pass))

  # --- Additive strengthening (beyond the plan's Task 3.2 block) -------------
  # `expect_true(res$pass)` alone would still pass if `pass` were hard-coded
  # TRUE, or if the ranges were never consulted. Pin the whole per_gene frame:
  # every observed frequency, and the low/high actually used, must be the
  # published anchor for that gene.
  expect_setequal(res$per_gene$gene, names(PUBLISHED_MUT_FREQ_RANGES))
  observed <- setNames(res$per_gene$observed, res$per_gene$gene)
  expect_equal(observed[["PBRM1"]], 0.35)
  expect_equal(observed[["SETD2"]], 0.12)
  expect_equal(observed[["BAP1"]], 0.10)
  for (i in seq_len(nrow(res$per_gene))) {
    rng <- PUBLISHED_MUT_FREQ_RANGES[[res$per_gene$gene[i]]]
    expect_identical(res$per_gene$low[i], unname(rng[["low"]]))
    expect_identical(res$per_gene$high[i], unname(rng[["high"]]))
  }
  expect_type(res$label, "character")

  # The DENOMINATOR must be reported. Without it, `mean()` over 20 rows and over
  # 417 rows are indistinguishable in the output, so a regression that collapsed
  # the cohort — while leaving the proportions inside the published bands, which
  # small samples easily do — would pass every anchor in the suite.
  expect_identical(res$n, n)

  # Immutability (coding-style): the input frame must be untouched.
  expect_identical(ncol(mut_annot), 5L)
  expect_identical(mut_annot$VHL, c(rep(TRUE, 50), rep(FALSE, 50)))
})

test_that("fn_check_mutation_freq flags a driver outside its published range", {
  # Arrange — VHL 0.20, below the 0.40-0.60 range
  n <- 100L
  mut_annot <- data.frame(
    sample_id = paste0("P", seq_len(n)),
    VHL   = c(rep(TRUE, 20), rep(FALSE, 80)),
    PBRM1 = c(rep(TRUE, 35), rep(FALSE, 65)),
    SETD2 = c(rep(TRUE, 12), rep(FALSE, 88)),
    BAP1  = c(rep(TRUE, 10), rep(FALSE, 90)),
    stringsAsFactors = FALSE
  )

  # Act
  res <- fn_check_mutation_freq(mut_annot)

  # Assert
  expect_false(res$pass)
  expect_false(res$per_gene$pass[res$per_gene$gene == "VHL"])

  # --- Additive strengthening ------------------------------------------------
  # Falsifiability is the whole point: exactly ONE gene must be flagged, and
  # only the out-of-range one. A check that fails everything on any deviation
  # is as uninformative as one that never fails.
  expect_equal(res$per_gene$observed[res$per_gene$gene == "VHL"], 0.20)
  failed <- res$per_gene$gene[!res$per_gene$pass]
  expect_identical(failed, "VHL")

  # A frequency ABOVE the published high bound must also fail (the check is
  # two-sided, not a floor).
  high_vhl <- mut_annot
  high_vhl$VHL <- c(rep(TRUE, 90), rep(FALSE, 10))
  res_high <- fn_check_mutation_freq(high_vhl)
  expect_false(res_high$pass)
  expect_identical(res_high$per_gene$gene[!res_high$per_gene$pass], "VHL")
})

test_that("fn_check_mutation_freq handles the real integer-0/1 mut_annot shape", {
  # Arrange — Module 1's fn_extract_mutation_status emits INTEGER 0/1 columns
  # plus rownames == sample_id, not the logical columns the plan's interface
  # note describes. The check must agree with the object the DAG actually
  # feeds it, or the anchor would be verified on a shape that never occurs.
  n <- 100L
  ids <- paste0("TCGA-B0-", sprintf("%04d", seq_len(n)))
  mut_annot <- data.frame(
    sample_id = ids,
    VHL   = c(rep(1L, 50), rep(0L, 50)),
    PBRM1 = c(rep(1L, 35), rep(0L, 65)),
    SETD2 = c(rep(1L, 12), rep(0L, 88)),
    BAP1  = c(rep(1L, 10), rep(0L, 90)),
    MTOR  = c(rep(1L,  6), rep(0L, 94)),   # in the panel, no published range
    stringsAsFactors = FALSE
  )
  rownames(mut_annot) <- ids

  # Act
  res <- fn_check_mutation_freq(mut_annot)

  # Assert
  expect_true(res$pass)
  expect_equal(res$per_gene$observed[res$per_gene$gene == "VHL"], 0.50)
  # MTOR has no published range, so it must not be scored (silently scoring it
  # against a missing anchor would produce NA and an unfalsifiable result).
  expect_false("MTOR" %in% res$per_gene$gene)
  expect_false(anyNA(res$per_gene$pass))
})

test_that("fn_check_mutation_freq is inclusive at the published bounds", {
  # Arrange — VHL exactly 0.40 (low bound) and exactly 0.60 (high bound).
  n <- 100L
  base <- data.frame(
    sample_id = paste0("P", seq_len(n)),
    PBRM1 = c(rep(1L, 35), rep(0L, 65)),
    SETD2 = c(rep(1L, 12), rep(0L, 88)),
    BAP1  = c(rep(1L, 10), rep(0L, 90)),
    stringsAsFactors = FALSE
  )
  at_low  <- cbind(base, VHL = c(rep(1L, 40), rep(0L, 60)))
  at_high <- cbind(base, VHL = c(rep(1L, 60), rep(0L, 40)))

  # Act / Assert
  expect_true(fn_check_mutation_freq(at_low)$pass)
  expect_true(fn_check_mutation_freq(at_high)$pass)
})

test_that("fn_check_mutation_freq refuses inputs it cannot score", {
  # A silent pass on an annotation carrying none of the ranged drivers would
  # be the worst failure mode: a green anchor computed from nothing.
  bad <- data.frame(sample_id = c("P1", "P2"), FOO = c(1L, 0L),
                    stringsAsFactors = FALSE)
  expect_error(fn_check_mutation_freq(bad), "none of the ranged driver genes")
  expect_error(fn_check_mutation_freq(list(sample_id = "P1")))
  expect_error(
    fn_check_mutation_freq(data.frame(VHL = 1L)),
    "sample_id"
  )
})

test_that("fn_check_mutation_freq honours the gene_panel it is given", {
  # `gene_panel` used to be declared and never read: MEASURED, a caller writing
  # gene_panel = "BAP1" still got a verdict over all four ranged genes. It was
  # documented as a KNOWN DEFECT, but documenting a trap does not disarm it, and
  # .lintr disables object_usage_linter for this project so nothing else would
  # ever flag it. Pin the behaviour in BOTH directions so a future change cannot
  # silently alter what the anchor scores.
  n <- 100L
  mut_annot <- data.frame(
    sample_id = paste0("P", seq_len(n)),
    VHL   = c(rep(1L, 50), rep(0L, 50)),
    PBRM1 = c(rep(1L, 35), rep(0L, 65)),
    SETD2 = c(rep(1L, 12), rep(0L, 88)),
    BAP1  = c(rep(1L, 10), rep(0L, 90)),
    stringsAsFactors = FALSE
  )

  narrowed <- fn_check_mutation_freq(mut_annot, gene_panel = c("VHL", "BAP1"))
  expect_identical(narrowed$per_gene$gene, c("VHL", "BAP1"))

  # The DEFAULT path is unchanged: DRIVER_GENE_PANEL is a superset of the ranged
  # genes, so the anchor still scores all four.
  expect_true(all(names(PUBLISHED_MUT_FREQ_RANGES) %in% DRIVER_GENE_PANEL))
  expect_setequal(fn_check_mutation_freq(mut_annot)$per_gene$gene,
                  names(PUBLISHED_MUT_FREQ_RANGES))

  # A panel that intersects nothing must stop, not silently score everything.
  expect_error(fn_check_mutation_freq(mut_annot, gene_panel = "MTOR"),
               "none of the ranged driver genes")
})

test_that("fn_check_mutation_freq refuses a column it cannot coerce", {
  # as.logical() on a character or factor 0/1 column returns ALL NA, and
  # na.rm = TRUE then averages an empty vector. MEASURED before this guard, with
  # a character VHL column: observed = NaN, VHL pass = NA, and the OVERALL
  # verdict = NA — a mutation-frequency anchor that scored nothing and reported
  # no failure. "A control that cannot fail proves nothing" (file header).
  n <- 100L
  base <- data.frame(
    sample_id = paste0("P", seq_len(n)),
    PBRM1 = c(rep(1L, 35), rep(0L, 65)),
    SETD2 = c(rep(1L, 12), rep(0L, 88)),
    BAP1  = c(rep(1L, 10), rep(0L, 90)),
    stringsAsFactors = FALSE
  )
  chr_col <- cbind(base, VHL = as.character(c(rep(1L, 50), rep(0L, 50))))
  expect_error(fn_check_mutation_freq(chr_col), "logical or 0/1 numeric")
  fct_col <- cbind(base, VHL = factor(c(rep(1L, 50), rep(0L, 50))))
  expect_error(fn_check_mutation_freq(fct_col), "logical or 0/1 numeric")

  # A column that IS numeric but carries no scorable value must also stop
  # rather than produce a NaN frequency.
  all_na <- cbind(base, VHL = rep(NA_integer_, n))
  expect_error(fn_check_mutation_freq(all_na), "not finite")

  # The verdict column can never contain NA on an input the function accepts.
  ok <- cbind(base, VHL = c(rep(1L, 50), rep(0L, 50)))
  expect_false(anyNA(fn_check_mutation_freq(ok)$per_gene$pass))
  expect_type(fn_check_mutation_freq(ok)$pass, "logical")
  expect_false(is.na(fn_check_mutation_freq(ok)$pass))
})

test_that("fn_check_bap1_survival detects worse OS in BAP1-mutant tumours", {
  # Arrange — mutants get a higher hazard, so shorter times
  set.seed(3)
  n <- 200L
  bap1 <- rbinom(n, 1L, 0.15)
  os_time  <- rexp(n, rate = 1 / 600 + bap1 * (1 / 250))
  os_event <- rbinom(n, 1L, 0.7)
  clinical <- data.frame(
    sample_id = paste0("P", seq_len(n)),
    os_time = os_time, os_event = os_event, stringsAsFactors = FALSE
  )
  mut_annot <- data.frame(
    sample_id = paste0("P", seq_len(n)),
    BAP1 = as.logical(bap1), stringsAsFactors = FALSE
  )

  # Act
  res <- fn_check_bap1_survival(clinical, mut_annot)

  # Assert
  expect_gt(res$hr, 1)
  expect_true(res$pass)
  expect_lt(res$p_value, 0.05)

  # --- Additive strengthening (beyond the plan's Task 3.3 block) -------------
  # The simulated hazard ratio is (1/600 + 1/250) / (1/600) = 3.4, so a
  # correctly specified Cox model must land near it. Asserting only "hr > 1"
  # would still pass if the model were fitted on the wrong column or the
  # coefficient read on the log scale.
  expect_gt(res$hr, 1.5)
  expect_lt(res$hr, 8)
  # The reported CI must actually bracket the point estimate and exclude the
  # null, consistent with p < 0.05.
  expect_lt(res$ci_low, res$hr)
  expect_gt(res$ci_high, res$hr)
  expect_gt(res$ci_low, 1)
  expect_identical(res$n, n)
  expect_type(res$label, "character")

  # Immutability (coding-style): neither input frame may be modified. In
  # particular no `bap1_mut` column may leak back into the caller's mut_annot.
  expect_named(clinical, c("sample_id", "os_time", "os_event"))
  expect_named(mut_annot, c("sample_id", "BAP1"))
})

test_that("fn_check_bap1_survival FAILS when BAP1-mutants survive LONGER", {
  # This is the assertion that makes the anchor worth having. A positive
  # control that returns pass = TRUE for any BAP1 effect — or for none — would
  # manufacture false confidence. Here the mutants are PROTECTED (lower
  # hazard), which contradicts the published direction, so pass MUST be FALSE.
  set.seed(11)
  n <- 200L
  bap1 <- rbinom(n, 1L, 0.30)
  # Mutant hazard 1/2000 vs wild-type 1/300 -> protective, HR well below 1.
  os_time  <- rexp(n, rate = ifelse(bap1 == 1L, 1 / 2000, 1 / 300))
  clinical <- data.frame(
    sample_id = paste0("P", seq_len(n)),
    os_time = os_time, os_event = rep(1L, n), stringsAsFactors = FALSE
  )
  mut_annot <- data.frame(
    sample_id = paste0("P", seq_len(n)),
    BAP1 = as.logical(bap1), stringsAsFactors = FALSE
  )

  res <- fn_check_bap1_survival(clinical, mut_annot)

  expect_lt(res$hr, 1)
  expect_false(res$pass)

  # DOCUMENTED LIMITATION, not an endorsement: `pass` is DIRECTIONAL ONLY
  # (pass = hr > 1) and carries no significance requirement, so a null result
  # with hr = 1.0001 would still set pass = TRUE. Statistical significance is
  # asserted separately by the Task 3.7 anchor (expect_lt(bs$p_value, 0.05)).
  # Anything consuming sanity_results$bap1_survival$pass as a standalone
  # verdict must check p_value alongside it.
  expect_true(is.numeric(res$p_value) && res$p_value >= 0 && res$p_value <= 1)
})

test_that("fn_check_bap1_survival drops unusable rows and honours the merge", {
  # Arrange — one non-overlapping sample on each side, one NA time, one
  # non-positive time. Only the clean intersection may reach the model.
  set.seed(5)
  n <- 120L
  ids <- paste0("P", seq_len(n))
  bap1 <- rbinom(n, 1L, 0.25)
  clinical <- data.frame(
    sample_id = c(ids, "ONLY_CLINICAL"),
    os_time   = c(rexp(n, rate = 1 / 500 + bap1 * (1 / 200)), 400),
    os_event  = c(rep(1L, n), 1L),
    stringsAsFactors = FALSE
  )
  clinical$os_time[1] <- NA_real_    # unusable: missing time
  clinical$os_time[2] <- 0           # unusable: non-positive time
  mut_annot <- data.frame(
    sample_id = c(ids, "ONLY_MUT"),
    BAP1 = as.logical(c(bap1, 1L)),
    stringsAsFactors = FALSE
  )

  # Act
  res <- fn_check_bap1_survival(clinical, mut_annot)

  # Assert — n = 120 shared ids, minus the NA row and the zero-time row.
  expect_identical(res$n, n - 2L)
  expect_gt(res$hr, 1)
})

test_that("fn_check_bap1_survival REFUSES to fit a degenerate Cox model", {
  # coxph does NOT error on a degenerate design — it returns quietly with NA
  # coefficients — so without a guard the anchor emits `pass = NA`, which is not
  # a verdict at all: `if (res$pass)` errors on it, and it violates the suite's
  # own shape anchor (`expect_false(is.na(sr[[nm]]$pass))`).
  #
  # MEASURED before this guard, on 300 cases:
  #   zero OS events    -> hr NA, p NA, pass NA (no error, no warning)
  #   BAP1 all wild-type-> hr NA, p NA, pass NA
  #   BAP1 all mutant   -> hr NA, p NA, pass NA
  #   a single OS event -> hr 3.7e+09, ci [0, Inf], p 0.999, pass TRUE
  # The last is the worst: a fit carrying no information reporting GREEN.
  #
  # R/constants.R:155-160 names the zero-event case as "the exact silent-green
  # failure mode the Module 3 suite exists to prevent", but the guard lived only
  # in the vital-status decode constant, never in the check itself.
  set.seed(21)
  n <- 300L
  ids <- paste0("P", seq_len(n))
  base_clin <- data.frame(
    sample_id = ids, os_time = rexp(n, rate = 1 / 500),
    os_event = rbinom(n, 1L, 0.5), stringsAsFactors = FALSE
  )
  base_mut <- data.frame(sample_id = ids, BAP1 = rbinom(n, 1L, 0.2),
                         stringsAsFactors = FALSE)

  # (i) A vital_status decode that matched nothing -> zero events.
  no_events <- base_clin
  no_events$os_event <- rep(0L, n)
  expect_error(fn_check_bap1_survival(no_events, base_mut), "0 OS events")
  expect_error(fn_check_bap1_survival(no_events, base_mut),
               "VITAL_STATUS_DEAD_VALUES")

  # (ii) Too few events to license any estimate, even though some exist.
  few_events <- base_clin
  few_events$os_event <- c(rep(1L, MIN_OS_EVENTS - 1L),
                           rep(0L, n - MIN_OS_EVENTS + 1L))
  expect_error(fn_check_bap1_survival(few_events, base_mut), "OS events")

  # (iii) BAP1 constant: no contrast exists, so nothing can be estimated.
  all_wt <- base_mut
  all_wt$BAP1 <- rep(0L, n)
  expect_error(fn_check_bap1_survival(base_clin, all_wt), "constant")
  all_mut <- base_mut
  all_mut$BAP1 <- rep(1L, n)
  expect_error(fn_check_bap1_survival(base_clin, all_mut), "constant")

  # The floor is far below the MEASURED cohort (173 OS events among 522 usable
  # main-cohort cases), so it can never mask a real result.
  expect_lt(MIN_OS_EVENTS, 173L)
  # Exactly at the floor the check must still run, not refuse.
  at_floor <- base_clin
  at_floor$os_event <- c(rep(1L, MIN_OS_EVENTS), rep(0L, n - MIN_OS_EVENTS))
  res <- fn_check_bap1_survival(at_floor, base_mut)
  expect_type(res$pass, "logical")
  expect_false(is.na(res$pass))
})

test_that("fn_check_bap1_survival never returns an NA verdict", {
  # `pass` must be a real logical(1) on every input the function ACCEPTS —
  # anything else is refused above rather than reported as NA.
  set.seed(23)
  n <- 250L
  bap1 <- rbinom(n, 1L, 0.2)
  clinical <- data.frame(
    sample_id = paste0("P", seq_len(n)),
    os_time = rexp(n, rate = 1 / 500 + bap1 * (1 / 300)),
    os_event = rbinom(n, 1L, 0.6), stringsAsFactors = FALSE
  )
  mut_annot <- data.frame(sample_id = paste0("P", seq_len(n)), BAP1 = bap1,
                          stringsAsFactors = FALSE)

  res <- fn_check_bap1_survival(clinical, mut_annot)

  expect_false(is.na(res$pass))
  expect_true(is.finite(res$hr))
  expect_true(is.finite(res$p_value))
  expect_gte(res$n_events, MIN_OS_EVENTS)
  expect_lte(res$n_events, res$n)
})

test_that("fn_check_bap1_survival refuses inputs missing required columns", {
  ok_clin <- data.frame(sample_id = "P1", os_time = 10, os_event = 1L,
                        stringsAsFactors = FALSE)
  ok_mut  <- data.frame(sample_id = "P1", BAP1 = TRUE, stringsAsFactors = FALSE)
  expect_error(fn_check_bap1_survival(ok_clin[, c("sample_id", "os_time")], ok_mut))
  expect_error(fn_check_bap1_survival(ok_clin, ok_mut[, "sample_id", drop = FALSE]),
               "BAP1")
})

test_that("fn_check_methyl_strata recovers four separated methylation strata", {
  # Arrange — four blocks with distinct global M-value offsets
  set.seed(1)
  n_per <- 15L; k <- 4L; n_cpg <- 30L
  offsets <- c(-3, -1, 1, 3)
  blocks <- lapply(seq_len(k), function(i) {
    matrix(rnorm(n_cpg * n_per, mean = offsets[i], sd = 0.3), nrow = n_cpg)
  })
  methyl_mat <- do.call(cbind, blocks)
  rownames(methyl_mat) <- paste0("cg", seq_len(n_cpg))
  colnames(methyl_mat) <- paste0("S", seq_len(k * n_per))

  # Act
  res <- fn_check_methyl_strata(methyl_mat)

  # Assert
  expect_identical(res$n_strata, 4L)
  expect_gt(res$silhouette, SANITY_MIN_SILHOUETTE)
  expect_lt(res$kw_p_value, SANITY_MAX_P)
  expect_true(res$pass)
})

test_that("fn_check_methyl_strata FAILS on methylation with no strata structure", {
  # THE test that makes this anchor worth having. If a check returns pass = TRUE
  # on structureless data it proves nothing. Here 1000 iid-noise CpGs over 80
  # samples contain no strata at all, so k-means still returns 4 non-empty
  # groups and the Kruskal test is still significant (both conjuncts are weak),
  # but the silhouette collapses and `pass` MUST be FALSE.
  set.seed(99)
  n_cpg <- 1000L; n_samp <- 80L
  methyl_mat <- matrix(rnorm(n_cpg * n_samp), nrow = n_cpg)
  rownames(methyl_mat) <- paste0("cg", seq_len(n_cpg))
  colnames(methyl_mat) <- paste0("S", seq_len(n_samp))

  res <- fn_check_methyl_strata(methyl_mat)

  expect_false(res$pass)
  expect_lt(res$silhouette, SANITY_MIN_SILHOUETTE)
  # DOCUMENTED WEAKNESS, asserted so it cannot be forgotten: neither of the
  # other two conjuncts discriminates. `n_strata` is tautological (k-means with
  # centers = 4 cannot return other than 4 non-empty groups; Hartigan-Wong
  # errors rather than emptying one), and the Kruskal test is CIRCULAR — the
  # clusters are derived from the same matrix whose column means it compares,
  # so it fires on pure noise too. The silhouette is the only falsifiable term.
  expect_identical(res$n_strata, 4L)
  expect_lt(res$kw_p_value, SANITY_MAX_P)
})

helper_platform_matrix <- function(seed, plat_offset, strata_sd = 0,
                                   n27 = 200L, n450 = 324L, n_cpg = 500L) {
  # A merged HM27+HM450 matrix. `plat_offset` is a per-platform mean shift (the
  # batch effect fn_merge_methyl_platforms does NOT correct); `strata_sd` adds
  # four TRUE biological strata spread evenly across BOTH platforms.
  set.seed(seed)
  n <- n27 + n450
  m <- matrix(stats::rnorm(n_cpg * n), nrow = n_cpg)
  g <- rep_len(1:4, n)
  if (strata_sd > 0) {
    bio <- matrix(stats::rnorm(n_cpg * 4L, 0, strata_sd), nrow = n_cpg)
    m <- m + bio[, g]
  }
  m[, seq_len(n27)] <- m[, seq_len(n27)] + plat_offset
  rownames(m) <- paste0("cg", seq_len(n_cpg))
  colnames(m) <- paste0("S", seq_len(n))
  platform <- factor(c(rep("HM27", n27), rep("HM450", n450)),
                     levels = METHYL_PLATFORMS)
  names(platform) <- colnames(m)
  list(m = m, platform = platform, strata = g)
}

test_that("fn_check_methyl_strata FAILS when the partition is just the assay", {
  # THE confound this check could not see. methyl_mat is cbind(HM27, HM450) with
  # NO batch correction (fn_merge_methyl_platforms only column-binds), and
  # platform is the strongest single axis in merged 27k/450k M-values. So the
  # top-line green light "TCGA KIRC methylation resolves into m1-m4 strata"
  # could be produced ENTIRELY by the assay rather than by the published m1-m4
  # biology, and nothing in the returned object, in `pass`, or in the Task 3.7
  # anchor could distinguish the two.
  #
  # MEASURED on constructed data — 500 CpGs, 200 HM27-like + 324 HM450-like
  # samples, iid noise, NO biological strata, only a per-platform mean offset:
  #   offset 1.0 SD -> pass FALSE, silhouette 0.072 (fails only narrowly)
  #   offset 1.5 SD -> pass TRUE,  silhouette 0.122
  #   offset 2.0 SD -> pass TRUE,  silhouette 0.164
  #   offset 3.0 SD -> pass TRUE,  silhouette 0.222
  # with kruskal p = 2.7e-80 and ARI(cluster, platform) = 0.504 throughout.
  for (off in c(1.5, 2.0, 3.0)) {
    d <- helper_platform_matrix(99, plat_offset = off)
    res <- fn_check_methyl_strata(d$m, platform = d$platform)

    # The silhouette still clears its floor — that term cannot see the confound.
    expect_gt(res$silhouette, SANITY_MIN_SILHOUETTE)
    # The platform term must, and must veto the verdict.
    expect_gt(res$platform_ari, SANITY_MAX_PLATFORM_ARI)
    expect_lt(res$platform_p, SANITY_MAX_P)
    expect_false(res$pass)
  }
})

test_that("fn_check_methyl_strata passes strata that CROSS the two platforms", {
  # The mirror image: four true strata distributed evenly over both platforms.
  # The partition then tracks biology, not the assay, so the platform term must
  # not fire and the check must report green.
  #
  # CALIBRATION behind SANITY_MAX_PLATFORM_ARI: when k-means recovered the true
  # strata (ARI vs truth 1.000) the platform ARI was -0.003; when it locked onto
  # the assay instead (ARI vs truth 0.299) the platform ARI was 0.532. The two
  # regimes are cleanly bimodal and 0.25 sits in the gap.
  d <- helper_platform_matrix(7, plat_offset = 0, strata_sd = 1.0)
  res <- fn_check_methyl_strata(d$m, platform = d$platform)

  expect_lt(res$platform_ari, SANITY_MAX_PLATFORM_ARI)
  expect_true(res$pass)
  expect_gt(mclust::adjustedRandIndex(res$cluster, d$strata), 0.9)

  # Even with a real platform offset present, biology strong enough to dominate
  # still reports green — the term rejects assay-driven partitions, not merged
  # data as such.
  d2 <- helper_platform_matrix(7, plat_offset = 1.0, strata_sd = 1.0)
  res2 <- fn_check_methyl_strata(d2$m, platform = d2$platform)
  expect_lt(res2$platform_ari, SANITY_MAX_PLATFORM_ARI)
  expect_true(res2$pass)
})

test_that("fn_check_methyl_strata validates the platform vector it is given", {
  d <- helper_platform_matrix(3, plat_offset = 1.0)
  expect_error(fn_check_methyl_strata(d$m, platform = d$platform[1:10]),
               "platform")
  misnamed <- d$platform
  names(misnamed) <- rev(names(misnamed))
  expect_error(fn_check_methyl_strata(d$m, platform = misnamed), "platform")

  # Without a platform vector the term is NOT silently satisfied — it is
  # reported as NA so the Task 3.7 anchor can refuse the result.
  res <- fn_check_methyl_strata(d$m)
  expect_true(is.na(res$platform_ari))
  expect_true(is.na(res$platform_p))
})

test_that("fn_check_methyl_strata tolerates the non-finite probes real data carries", {
  # The real Module 1 `methyl_mat` is 91.4% complete (432/5000 CpGs carry a
  # non-finite value from failed probes). stats::kmeans ERRORS on NA/NaN/Inf, so
  # without a completeness guard this anchor cannot run on the data it exists to
  # check. Same failure mode already fixed once for SNF (fn_complete_features).
  set.seed(1)
  n_per <- 15L; k <- 4L; n_cpg <- 30L
  offsets <- c(-3, -1, 1, 3)
  blocks <- lapply(seq_len(k), function(i) {
    matrix(rnorm(n_cpg * n_per, mean = offsets[i], sd = 0.3), nrow = n_cpg)
  })
  methyl_mat <- do.call(cbind, blocks)
  rownames(methyl_mat) <- paste0("cg", seq_len(n_cpg))
  colnames(methyl_mat) <- paste0("S", seq_len(k * n_per))
  dirty <- methyl_mat
  dirty[1L, 3L]  <- NA_real_
  dirty[2L, 10L] <- Inf
  dirty[3L, 1L]  <- NaN

  res <- fn_check_methyl_strata(dirty)

  # The three damaged CpGs are dropped; the other 27 still resolve the strata.
  expect_identical(res$n_cpg_used, 27L)
  expect_true(res$pass)
  expect_gt(res$silhouette, SANITY_MIN_SILHOUETTE)
  expect_identical(res$n_strata, 4L)
  # Every sample is retained: dropping is per-CpG, never per-sample.
  expect_length(res$cluster, ncol(dirty))
  # Immutability: the caller's matrix is untouched.
  expect_true(is.na(dirty[1L, 3L]))
  expect_identical(dim(dirty), c(n_cpg, k * n_per))
})

test_that("fn_check_methyl_strata refuses input it cannot honestly cluster", {
  # A silent pass computed from almost nothing is the worst failure mode, so a
  # matrix gutted by missingness must stop rather than cluster its remnants.
  set.seed(7)
  m <- matrix(rnorm(20L * 30L), nrow = 20L)
  rownames(m) <- paste0("cg", seq_len(20L))
  colnames(m) <- paste0("S", seq_len(30L))
  gutted <- m
  gutted[1:18, 1L] <- NA_real_          # one bad sample kills 18/20 CpGs
  expect_error(fn_check_methyl_strata(gutted), "complete CpGs")

  # Structural guards from the plan's own stopifnot.
  expect_error(fn_check_methyl_strata(as.data.frame(m)))
  expect_error(fn_check_methyl_strata(m[, 1:3, drop = FALSE]))
})

test_that("fn_check_methyl_strata is deterministic and returns a labelled object", {
  # A seeded check that drifts between runs cannot anchor anything.
  set.seed(1)
  n_per <- 15L; k <- 4L; n_cpg <- 30L
  offsets <- c(-3, -1, 1, 3)
  blocks <- lapply(seq_len(k), function(i) {
    matrix(rnorm(n_cpg * n_per, mean = offsets[i], sd = 0.3), nrow = n_cpg)
  })
  methyl_mat <- do.call(cbind, blocks)
  rownames(methyl_mat) <- paste0("cg", seq_len(n_cpg))
  colnames(methyl_mat) <- paste0("S", seq_len(k * n_per))

  a <- fn_check_methyl_strata(methyl_mat)
  b <- fn_check_methyl_strata(methyl_mat)

  expect_identical(a$cluster, b$cluster)
  expect_identical(a$silhouette, b$silhouette)
  expect_type(a$label, "character")
  expect_named(a$cluster, colnames(methyl_mat))
  # k = 4 is the published stratum count, taken from the constant, not inlined.
  expect_identical(a$n_strata, METHYL_N_STRATA)
})

test_that("fn_check_ccab_signature separates ccA-like and ccB-like groups", {
  # Arrange — cols 1:20 ccB-high (proliferative), cols 21:40 ccA-high
  set.seed(2)
  n <- 40L; half <- n / 2L
  ccb <- CCB_PROLIFERATION_MARKERS
  cca <- CCA_ANGIOGENESIS_MARKERS
  mk <- function(hi_cols) {
    v <- rnorm(n, mean = 4, sd = 0.3)
    v[hi_cols] <- rnorm(length(hi_cols), mean = 8, sd = 0.3)
    v
  }
  rna_mat <- rbind(
    t(sapply(ccb, function(g) mk(seq_len(half)))),
    t(sapply(cca, function(g) mk((half + 1L):n)))
  )
  colnames(rna_mat) <- paste0("T", seq_len(n))

  # Act
  res <- fn_check_ccab_signature(rna_mat)

  # Assert
  expect_gt(res$silhouette, SANITY_MIN_SILHOUETTE_2D)
  expect_lt(res$separation_p_value, SANITY_MAX_P)
  expect_true(res$pass)
})

# Helper for the ccA/ccB negative controls: a genes x samples matrix carrying
# every published marker but NO ccA/ccB structure.
make_structureless_rna <- function(seed, n = 40L, shared_axis = FALSE) {
  set.seed(seed)
  genes <- c(CCB_PROLIFERATION_MARKERS, CCA_ANGIOGENESIS_MARKERS)
  m <- if (shared_axis) {
    # One common axis: proliferation and angiogenesis move TOGETHER, the
    # opposite of the published ccA/ccB opposition.
    shared <- rnorm(n, mean = 4, sd = 1.5)
    t(vapply(seq_along(genes), function(i) shared + rnorm(n, 0, 0.2),
             numeric(n)))
  } else {
    matrix(rnorm(length(genes) * n, mean = 4, sd = 0.3), nrow = length(genes))
  }
  rownames(m) <- genes
  colnames(m) <- paste0("T", seq_len(n))
  m
}

test_that("fn_check_ccab_signature FAILS on expression with no ccA/ccB structure", {
  # THE falsifiability test. MEASURED: the plan's silhouette + Wilcoxon pair
  # alone returns pass = TRUE on 3 of these 5 structureless matrices (seeds
  # 2, 4, 5), because 2-D k-means always splits a blob (silhouette 0.37-0.44)
  # and the Wilcoxon is CIRCULAR — it compares the ccB-ccA axis across clusters
  # built from that same pair of scores. A control that greenlights pure noise
  # manufactures false confidence, so `pass` must also require the published
  # ccA/ccB OPPOSITION, which noise does not produce.
  for (seed in 1:5) {
    res <- fn_check_ccab_signature(make_structureless_rna(seed))
    expect_false(res$pass)
  }
})

test_that("the 2-D ccA/ccB silhouette is judged against a CALIBRATED threshold", {
  # SANITY_MIN_SILHOUETTE = 0.10 was applied to two clusterings with completely
  # different null distributions: a ~4568-dimensional 4-means (methylation,
  # where iid noise gives ~0.005) and this 2-dimensional 2-means. In 2-D,
  # k-means always splits a blob, so the statistic is large under the NULL and
  # a 0.10 floor cannot fail on any input the function accepts — it was carried
  # on the dashboard as if it were evidence.
  #
  # CALIBRATION (this file, R 4.6.0): 200 structureless matrices at each of
  # n = 40 / 100 / 524, all 12 published markers, iid noise, NO ccA/ccB
  # structure -> silhouette min 0.284, median 0.318-0.350, MAX 0.466.
  # Genuine opposition at the weakest effect tested (a 1-SD-unit shift between
  # the two halves) -> 0.490-0.562. SANITY_MIN_SILHOUETTE_2D = 0.50 sits above
  # the measured null ceiling and below the weakest real structure.
  noise_sil <- vapply(1:20, function(s) {
    fn_check_ccab_signature(make_structureless_rna(s, n = 200L))$silhouette
  }, numeric(1))

  # The old constant is NON-DISCRIMINATING here: pure noise clears it every time.
  expect_true(all(noise_sil > SANITY_MIN_SILHOUETTE))
  # The calibrated one is not.
  expect_true(all(noise_sil < SANITY_MIN_SILHOUETTE_2D))
  expect_gt(SANITY_MIN_SILHOUETTE_2D, SANITY_MIN_SILHOUETTE)

  # Genuine ccA/ccB opposition must still clear it, or the threshold would be
  # unreachable rather than strict.
  set.seed(31)
  n <- 200L; half <- n / 2L
  mk <- function(hi) {
    v <- stats::rnorm(n, 4, 0.3)
    v[hi] <- stats::rnorm(length(hi), 4.3, 0.3)   # a 1-SD-unit shift
    v
  }
  real <- rbind(
    t(vapply(CCB_PROLIFERATION_MARKERS, function(g) mk(seq_len(half)), numeric(n))),
    t(vapply(CCA_ANGIOGENESIS_MARKERS, function(g) mk((half + 1L):n), numeric(n)))
  )
  colnames(real) <- paste0("T", seq_len(n))
  res <- fn_check_ccab_signature(real)
  expect_gt(res$silhouette, SANITY_MIN_SILHOUETTE_2D)
  expect_true(res$pass)
})

test_that("fn_check_ccab_signature FAILS when the two programmes move together", {
  # Brannon 2010's claim is an OPPOSITION: ccB tumours are proliferative and
  # angiogenesis-low, ccA the reverse. Markers that all rise and fall together
  # are the clearest contradiction of it, so this must never report green.
  for (seed in 1:3) {
    res <- fn_check_ccab_signature(make_structureless_rna(seed, shared_axis = TRUE))
    expect_false(res$pass)
    expect_gt(res$anticorr_rho, 0)
  }
})

test_that("fn_check_ccab_signature reports the ccA/ccB opposition it tested", {
  # The pass flag must be traceable to published quantities, not opaque.
  set.seed(2)
  n <- 40L; half <- n / 2L
  mk <- function(hi_cols) {
    v <- rnorm(n, mean = 4, sd = 0.3)
    v[hi_cols] <- rnorm(length(hi_cols), mean = 8, sd = 0.3)
    v
  }
  rna_mat <- rbind(
    t(sapply(CCB_PROLIFERATION_MARKERS, function(g) mk(seq_len(half)))),
    t(sapply(CCA_ANGIOGENESIS_MARKERS, function(g) mk((half + 1L):n)))
  )
  colnames(rna_mat) <- paste0("T", seq_len(n))

  res <- fn_check_ccab_signature(rna_mat)

  # ccB-high tumours are cols 1:20 by construction, so the axis must be high
  # there and low in cols 21:40 — not merely "different between two clusters".
  expect_gt(mean(res$axis_score[seq_len(half)]),
            mean(res$axis_score[(half + 1L):n]))
  expect_lt(res$anticorr_rho, 0)
  expect_lt(res$anticorr_p_value, SANITY_MAX_P)
  expect_length(res$group, n)
  expect_named(res$axis_score, colnames(rna_mat))
  expect_type(res$label, "character")
  # Immutability: the caller's expression matrix is untouched.
  expect_identical(dim(rna_mat), c(12L, n))
  expect_identical(rownames(rna_mat),
                   c(CCB_PROLIFERATION_MARKERS, CCA_ANGIOGENESIS_MARKERS))
})

test_that("fn_check_ccab_signature reports and floors the panel it used", {
  # The check intersected the published panels with rownames(rna_mat) and then
  # required only 2 survivors per side, recording NOWHERE how many actually
  # survived. The returned object was indistinguishable between a full 6-vs-6
  # run and a gutted 2-vs-2 one, so no anchor could detect attrition — while
  # MEASURED false-green rates rise as the panel shrinks (on 20 structureless
  # matrices the 6+6 panel returned pass = TRUE 0/20, the 2+2 panel 1/20).
  #
  # This is inconsistent with the rest of the suite: fn_check_methyl_strata
  # reports n_cpg_used and the mutation anchor asserts the exact gene set.
  full <- make_structureless_rna(3)

  res <- fn_check_ccab_signature(full)

  expect_identical(res$n_ccb_used, length(CCB_PROLIFERATION_MARKERS))
  expect_identical(res$n_cca_used, length(CCA_ANGIOGENESIS_MARKERS))
  expect_setequal(res$markers_used$ccb, CCB_PROLIFERATION_MARKERS)
  expect_setequal(res$markers_used$cca, CCA_ANGIOGENESIS_MARKERS)

  # A gutted panel must stop LOUDLY rather than score a published axis from a
  # couple of genes. Three per side is below the floor.
  three_each <- full[c(CCB_PROLIFERATION_MARKERS[1:3],
                       CCA_ANGIOGENESIS_MARKERS[1:3]), , drop = FALSE]
  expect_error(fn_check_ccab_signature(three_each),
               "insufficient ccA/ccB marker genes")

  # At the floor it still runs, and reports the reduced panel honestly.
  at_floor <- full[c(CCB_PROLIFERATION_MARKERS[1:SANITY_MIN_MARKERS_PER_PANEL],
                     CCA_ANGIOGENESIS_MARKERS[1:SANITY_MIN_MARKERS_PER_PANEL]), ,
                   drop = FALSE]
  res_floor <- fn_check_ccab_signature(at_floor)
  expect_identical(res_floor$n_ccb_used, SANITY_MIN_MARKERS_PER_PANEL)
  expect_identical(res_floor$n_cca_used, SANITY_MIN_MARKERS_PER_PANEL)
  expect_gte(SANITY_MIN_MARKERS_PER_PANEL, 4L)
})

test_that("fn_check_ccab_signature scores the partition it actually built", {
  # The check clustered in STANDARDISED space but measured the silhouette in RAW
  # space, so the reported number was not the separation quality of the
  # clustering being reported. When the two marker scores have different spreads
  # the two quantities diverge badly.
  #
  # MEASURED on this matrix before the fix — ccB carries a 37x larger spread but
  # no structure, ccA carries the split — reported silhouette -0.020 versus
  # 0.538 recomputed in the space k-means used. fn_check_methyl_strata already
  # gets this right (it clusters and measures on the same `feat`).
  set.seed(11)
  n <- 60L
  ccb <- t(vapply(CCB_PROLIFERATION_MARKERS,
                  function(g) stats::rnorm(n, 100, 50), numeric(n)))
  cca <- t(vapply(CCA_ANGIOGENESIS_MARKERS,
                  function(g) c(stats::rnorm(n / 2, -0.5, 0.05),
                                stats::rnorm(n / 2, 0.5, 0.05)), numeric(n)))
  rna_mat <- rbind(ccb, cca)
  colnames(rna_mat) <- paste0("T", seq_len(n))

  res <- fn_check_ccab_signature(rna_mat)

  # Independently rebuild the clustering and score it in ITS OWN space.
  b <- colMeans(rna_mat[CCB_PROLIFERATION_MARKERS, , drop = FALSE])
  a <- colMeans(rna_mat[CCA_ANGIOGENESIS_MARKERS, , drop = FALSE])
  feat <- cbind(scale(b), scale(a))
  set.seed(SANITY_SEED)
  km <- stats::kmeans(feat, centers = 2L, nstart = 25L)
  expected_sil <- mean(
    cluster::silhouette(km$cluster, stats::dist(feat))[, "sil_width"]
  )

  # The spreads really do differ by more than an order of magnitude here, so
  # this is a discriminating comparison, not a coincidence.
  expect_gt(stats::sd(b) / stats::sd(a), 10)
  expect_equal(res$silhouette, expected_sil)
  expect_identical(unname(res$group), unname(km$cluster))
})

test_that("fn_check_ccab_signature refuses inputs it cannot score", {
  # Scoring an axis from one marker, or from none, would be a green computed
  # from nothing — the failure mode this whole suite exists to prevent.
  set.seed(4)
  full <- make_structureless_rna(4)
  expect_error(fn_check_ccab_signature(full[1:7, , drop = FALSE]),
               "insufficient ccA/ccB marker genes")
  no_names <- full
  rownames(no_names) <- NULL
  expect_error(fn_check_ccab_signature(no_names))
  expect_error(fn_check_ccab_signature(as.data.frame(full)))
})

test_that("the seeded checks leave the caller's RNG stream intact", {
  # The file header claims "Each fn_check_* is pure". Both seeded checks called
  # set.seed() on the GLOBAL stream and never restored it, so after
  # sanity_results built, the session was left at the seed-42 state and ANY
  # later randomised code — other targets in the same process, or subsequent
  # test_that blocks that draw randoms without re-seeding — was silently coupled
  # to Module 3 and to whether Module 3 ran at all. That also makes test
  # outcomes order-dependent.
  #
  # MEASURED before the fix: .Random.seed differed before vs after each call,
  # and after set.seed(7) a caller's rnorm(3)[1] was 2.287 without the check but
  # 0.581 with it.
  set.seed(1)
  m <- matrix(rnorm(30L * 60L), nrow = 30L)
  rownames(m) <- paste0("cg", seq_len(30L))
  colnames(m) <- paste0("S", seq_len(60L))
  rna <- make_structureless_rna(1)

  set.seed(101)
  before <- .Random.seed
  invisible(fn_check_methyl_strata(m))
  expect_identical(.Random.seed, before)
  invisible(fn_check_ccab_signature(rna))
  expect_identical(.Random.seed, before)

  # The caller's DRAWS must be identical whether or not the check ran.
  set.seed(7)
  without <- stats::rnorm(3)
  set.seed(7)
  invisible(fn_check_ccab_signature(rna))
  with_check <- stats::rnorm(3)
  expect_identical(without, with_check)

  set.seed(7)
  invisible(fn_check_methyl_strata(m))
  expect_identical(stats::rnorm(3), without)

  # Determinism of the checks themselves is unaffected: set.seed still runs
  # before kmeans, only the caller's stream is restored afterwards.
  expect_identical(fn_check_methyl_strata(m)$cluster,
                   fn_check_methyl_strata(m)$cluster)
  expect_identical(fn_check_ccab_signature(rna)$silhouette,
                   fn_check_ccab_signature(rna)$silhouette)
})

# --- Power arithmetic behind the BAP1 re-specification ----------------------

test_that("fn_schoenfeld_events reproduces the arithmetic the BAP1 anchor rests on", {
  # The number that licensed re-specifying the BAP1 significance demand. If this
  # drifts, the justification recorded in R/constants.R and in the BAP1 anchor
  # is no longer true and both must be rewritten.
  #
  # d = (z_0.975 + z_0.80)^2 / (p (1-p) log(HR)^2)
  #   = (1.959964 + 0.8416212)^2 / (0.0863309 * 0.9136691 * 0.4599817^2)
  expect_equal(fn_schoenfeld_events(1.583991, 36 / 417), 470.3656,
               tolerance = 1e-4)

  # Hand-computed independently of the implementation.
  z <- stats::qnorm(0.975) + stats::qnorm(0.8)
  p <- 36 / 417
  expect_equal(fn_schoenfeld_events(1.583991, p),
               z^2 / (p * (1 - p) * log(1.583991)^2))

  # Monotonic in the right directions: weaker effects and rarer exposures both
  # cost events. At the low edge of the published band the requirement is ~3000
  # events, i.e. seven times the whole TCGA KIRC series.
  expect_gt(fn_schoenfeld_events(PUBLISHED_BAP1_HR_RANGE[["low"]], p),
            fn_schoenfeld_events(1.583991, p))
  expect_lt(fn_schoenfeld_events(PUBLISHED_BAP1_HR_RANGE[["high"]], p),
            fn_schoenfeld_events(1.583991, p))
  expect_gt(fn_schoenfeld_events(1.583991, 0.02),
            fn_schoenfeld_events(1.583991, 0.20))

  # A null effect is never detectable, and that is the correct limit, not an
  # error to be papered over.
  expect_identical(fn_schoenfeld_events(1, 0.2), Inf)

  # Refuses nonsense rather than returning a plausible-looking number.
  expect_error(fn_schoenfeld_events(0, 0.2))
  expect_error(fn_schoenfeld_events(1.5, 0))
  expect_error(fn_schoenfeld_events(1.5, 1))
  expect_error(fn_schoenfeld_events(c(1.5, 2), 0.2))
})

test_that("fn_check_bap1_survival reports the design adequacy of its own fit", {
  # The under-power must be readable off the verdict object, not reconstructed
  # by whoever happens to look at the p-value.
  set.seed(31)
  n <- 400L
  bap1 <- c(rep(1L, 40L), rep(0L, n - 40L))
  clinical <- data.frame(
    sample_id = paste0("P", seq_len(n)),
    os_time = rexp(n, rate = 1 / 600 + bap1 * (1 / 1200)),
    os_event = rbinom(n, 1L, 0.35), stringsAsFactors = FALSE
  )
  mut_annot <- data.frame(sample_id = paste0("P", seq_len(n)), BAP1 = bap1,
                          stringsAsFactors = FALSE)

  res <- fn_check_bap1_survival(clinical, mut_annot)

  expect_identical(res$n_mutant, 40L)
  expect_equal(res$mutant_frac, res$n_mutant / res$n)
  expect_lte(res$n_events_mutant, res$n_events)
  expect_equal(res$events_required,
               fn_schoenfeld_events(res$hr, res$mutant_frac))
  expect_identical(res$underpowered, res$n_events < res$events_required)
  # A thin mutant arm against a modest effect: this construction is underpowered
  # by the same arithmetic as the real cohort.
  expect_true(res$underpowered)

  # `pass` is still the DIRECTION and nothing else — the power fields are
  # reported, never folded into the verdict.
  expect_identical(res$pass, res$hr > 1)
})

# --- Within-platform evidence for the m1-m4 verdict -------------------------

test_that("fn_within_platform_silhouette separates assay-driven from real structure", {
  # The discriminator that makes the m1-m4 red state informative. On an
  # assay-only matrix the merged silhouette must EXCEED both arms (the gap
  # between platforms is doing the separating); with true strata spread over
  # both platforms it must not.
  assay_only <- helper_platform_matrix(99, plat_offset = 2.0)
  res_assay <- fn_check_methyl_strata(assay_only$m,
                                      platform = assay_only$platform)

  expect_s3_class(res_assay$within_platform, "data.frame")
  expect_setequal(res_assay$within_platform$platform, METHYL_PLATFORMS)
  expect_identical(sum(res_assay$within_platform$n), ncol(assay_only$m))
  # Neither arm contains any structure at all: the whole silhouette is the gap.
  expect_lt(max(res_assay$within_platform$silhouette), SANITY_MIN_SILHOUETTE)
  expect_true(res_assay$merged_exceeds_within)
  expect_false(res_assay$pass)

  real_strata <- helper_platform_matrix(7, plat_offset = 0, strata_sd = 1.0)
  res_bio <- fn_check_methyl_strata(real_strata$m,
                                    platform = real_strata$platform)
  expect_false(res_bio$merged_exceeds_within)
  expect_gt(min(res_bio$within_platform$silhouette), SANITY_MIN_SILHOUETTE)
  expect_true(res_bio$pass)

  # ... and it still holds with a genuine platform offset ON TOP of genuine
  # biology, which is the case the merged check must not reject.
  both <- helper_platform_matrix(7, plat_offset = 1.0, strata_sd = 1.0)
  res_both <- fn_check_methyl_strata(both$m, platform = both$platform)
  expect_false(res_both$merged_exceeds_within)
  expect_true(res_both$pass)
})

test_that("the within-platform term can only turn a verdict RED", {
  # It introduces no new threshold — it compares the check's own statistic to
  # itself computed inside each arm — so it must never rescue a partition that
  # the existing floors reject. Structureless data stays FALSE either way.
  set.seed(77)
  n_cpg <- 400L
  n27 <- 60L
  n450 <- 80L
  m <- matrix(rnorm(n_cpg * (n27 + n450)), nrow = n_cpg)
  rownames(m) <- paste0("cg", seq_len(n_cpg))
  colnames(m) <- paste0("S", seq_len(n27 + n450))
  plat <- factor(c(rep("HM27", n27), rep("HM450", n450)),
                 levels = METHYL_PLATFORMS)
  names(plat) <- colnames(m)

  res <- fn_check_methyl_strata(m, platform = plat)

  expect_false(res$pass)
  expect_lt(res$silhouette, SANITY_MIN_SILHOUETTE)
  expect_type(res$message, "character")
  expect_true(grepl("NOT recovered", res$message, fixed = TRUE))

  # With no platform vector the evidence is absent rather than fabricated, and
  # the term drops out of `pass` exactly as platform_ari does.
  bare <- fn_check_methyl_strata(m)
  expect_null(bare$within_platform)
  expect_true(is.na(bare$merged_exceeds_within))
})

test_that("the m1-m4 verdict states its finding in words", {
  # A bare FALSE invites a future reader to move a threshold. The object must
  # say what it found and why.
  d <- helper_platform_matrix(99, plat_offset = 2.0)
  red <- fn_check_methyl_strata(d$m, platform = d$platform)

  expect_false(red$pass)
  expect_true(grepl("tracks ASSAY PLATFORM", red$message, fixed = TRUE))
  expect_true(grepl("m1-m4 strata are NOT recovered", red$message, fixed = TRUE))
  expect_true(grepl("MERGED silhouette EXCEEDS both", red$message, fixed = TRUE))
  # The per-arm numbers travel with the sentence.
  expect_true(grepl("HM27", red$message, fixed = TRUE))
  expect_true(grepl("HM450", red$message, fixed = TRUE))

  green <- fn_check_methyl_strata(
    helper_platform_matrix(7, plat_offset = 0, strata_sd = 1.0)$m,
    platform = helper_platform_matrix(7, plat_offset = 0, strata_sd = 1.0)$platform
  )
  expect_true(green$pass)
  expect_true(grepl("recovered", green$message, fixed = TRUE))
  expect_false(grepl("NOT recovered", green$message, fixed = TRUE))
})

test_that("a veto-only m1-m4 red does NOT blame the silhouette floor", {
  # THE CASE THE MESSAGE USED TO MISATTRIBUTE. When the cluster-vs-platform ARI
  # sits UNDER its ceiling but the merged silhouette still beats both arms,
  # `merged_exceeds_within` is the sole failing term. Before this branch existed
  # the function fell through to "m1-m4 strata NOT recovered (silhouette %.4f,
  # floor %.2f)" — naming a floor the partition had CLEARED (0.1197 > 0.10) as
  # the cause, which is the reading that invites someone to lower the floor.
  #
  # Unreachable on the frozen snapshot (ARI 0.583 fires the platform branch), so
  # the message function is called DIRECTLY rather than through
  # fn_check_methyl_strata. That is the point: it becomes reachable as soon as
  # platform handling changes downstream.
  within <- data.frame(
    platform   = METHYL_PLATFORMS,
    n          = c(214L, 310L),
    n_cpg      = c(4658L, 4905L),
    silhouette = c(0.0858, 0.0489),
    stringsAsFactors = FALSE
  )
  mean_sil <- 0.1197
  clean_ari <- 0.10

  # Arrange the premise explicitly, so this test fails loudly if the constants
  # move underneath it rather than silently testing a different case.
  expect_lt(clean_ari, SANITY_MAX_PLATFORM_ARI)
  expect_gt(mean_sil, SANITY_MIN_SILHOUETTE)

  msg <- fn_methyl_strata_message(
    pass = FALSE, mean_sil = mean_sil, platform_ari = clean_ari,
    max_platform_ari = SANITY_MAX_PLATFORM_ARI, within_platform = within,
    merged_exceeds_within = TRUE
  )

  # It must NOT read as a floor breach.
  expect_false(grepl(sprintf("silhouette %.4f, floor", mean_sil), msg,
                     fixed = TRUE))
  # It must name the actual finding, and say the floor is not the reason.
  expect_true(grepl("separates the ASSAYS", msg, fixed = TRUE))
  expect_true(grepl("is NOT the reason", msg, fixed = TRUE))
  expect_true(grepl("stays RED", msg, fixed = TRUE))
  # The per-arm evidence still travels with it.
  expect_true(grepl("HM27", msg, fixed = TRUE))
  expect_true(grepl("HM450", msg, fixed = TRUE))

  # ... and the floor message is still emitted when the floor IS the breach.
  floor_msg <- fn_methyl_strata_message(
    pass = FALSE, mean_sil = 0.004, platform_ari = clean_ari,
    max_platform_ari = SANITY_MAX_PLATFORM_ARI, within_platform = within,
    merged_exceeds_within = FALSE
  )
  expect_true(grepl("floor", floor_msg, fixed = TRUE))
  expect_true(grepl("NOT recovered", floor_msg, fixed = TRUE))
})

# --- Subtype platform-cleanliness guard -------------------------------------

test_that("fn_check_subtype_platform passes an assignment independent of the assay", {
  # MEASURED shape on the real snapshot: ARI 0.0058 over 524 cases, S1 11/9,
  # S2 120/186, S3 33/43, S4 50/72 across HM27/HM450.
  set.seed(101)
  n <- 524L
  ids <- paste0("P", seq_len(n))
  platform <- stats::setNames(
    factor(c(rep("HM27", 214L), rep("HM450", 310L)), levels = METHYL_PLATFORMS),
    ids
  )
  # Subtypes assigned WITHOUT reference to platform.
  subtypes <- stats::setNames(
    factor(sample(paste0("S", 1:4), n, replace = TRUE)), ids
  )

  res <- fn_check_subtype_platform(subtypes, platform)

  expect_true(res$pass)
  expect_lt(abs(res$ari), SUBTYPE_MAX_PLATFORM_ARI)
  expect_identical(res$n, n)
  expect_identical(res$n_subtypes, 4L)
  expect_identical(sum(res$cross_tab), n)
  expect_type(res$label, "character")
})

test_that("fn_check_subtype_platform FAILS when the subtypes ARE the platform", {
  # THE test that makes this guard worth having. Individual MOFA factors reach
  # AUC 0.888 against platform on this snapshot, so a subtype assignment that
  # simply recovered the assay is a live possibility, not a hypothetical.
  ids <- paste0("P", seq_len(524L))
  platform <- stats::setNames(
    factor(c(rep("HM27", 214L), rep("HM450", 310L)), levels = METHYL_PLATFORMS),
    ids
  )
  # S1/S2 inside HM27, S3/S4 inside HM450: the subtype split IS the assay split.
  subtypes <- stats::setNames(
    factor(c(rep(c("S1", "S2"), length.out = 214L),
             rep(c("S3", "S4"), length.out = 310L))),
    ids
  )

  res <- fn_check_subtype_platform(subtypes, platform)

  expect_false(res$pass)
  expect_gt(res$ari, SUBTYPE_MAX_PLATFORM_ARI)
  expect_lt(res$p_value, SANITY_MAX_P)

  # A perfect 1:1 correspondence is the extreme of the same failure.
  perfect <- stats::setNames(factor(as.character(platform)), ids)
  expect_gt(fn_check_subtype_platform(perfect, platform)$ari, 0.9)
  expect_false(fn_check_subtype_platform(perfect, platform)$pass)
})

test_that("fn_check_subtype_platform refuses comparisons it cannot honestly score", {
  ids <- paste0("P", seq_len(100L))
  platform <- stats::setNames(
    factor(c(rep("HM27", 40L), rep("HM450", 60L)), levels = METHYL_PLATFORMS),
    ids
  )
  subtypes <- stats::setNames(factor(rep(paste0("S", 1:4), 25L)), ids)

  # Unnamed vectors would be compared in whatever order they arrive in.
  expect_error(fn_check_subtype_platform(unname(subtypes), platform), "NAMED")
  expect_error(fn_check_subtype_platform(subtypes, unname(platform)), "NAMED")

  # A partial overlap would score the guard on a silently reduced cohort, and a
  # near-zero ARI from a handful of cases is free.
  expect_error(fn_check_subtype_platform(subtypes[1:50], platform),
               "different samples")

  # A constant on either side gives ARI 0 for nothing.
  const_sub <- stats::setNames(factor(rep("S1", 100L)), ids)
  expect_error(fn_check_subtype_platform(const_sub, platform), "two levels")
  const_plat <- stats::setNames(factor(rep("HM27", 100L)), ids)
  expect_error(fn_check_subtype_platform(subtypes, const_plat), "two levels")

  # Order must not matter once both sides are named: the function aligns by id.
  shuffled <- platform[sample(ids)]
  expect_equal(fn_check_subtype_platform(subtypes, shuffled)$ari,
               fn_check_subtype_platform(subtypes, platform)$ari)
})

# --- Credibility anchor: real pipeline results vs published ccRCC literature -
# Reads the frozen sanity_results target. Executes wherever the _targets store
# is present (locally after tar_make; in CI after the release-asset restore).
#
# STORE RESOLUTION — a deliberate departure from the plan's one-liner
# `tryCatch(targets::tar_read(sanity_results), error = function(e) NULL)`:
#
#  1. testthat runs with the working directory set to tests/testthat/, so a
#     bare tar_read() looks for ./_targets, never finds it, and SKIPS — INCLUDING
#     in the container, where the release-asset store has been restored and the
#     anchor is supposed to run for real. An anchor that stays green by silently
#     skipping is exactly the false confidence this suite exists to prevent, so
#     the store is resolved explicitly against the repo root.
#  2. A blanket tryCatch also turns an ERRORED target into a skip. If the store
#     records sanity_results as failed, that is a finding and must surface as a
#     FAILURE, never as a skip.
#  3. A store that HAS pipeline metadata but no sanity_results ROW is likewise a
#     failure, not a skip. `sanity_results` is the LAST target in _targets.R, so
#     any upstream failure (a MOFA OOM, an ExperimentHub timeout, a
#     `continue-on-error: true` tar_make step) leaves a fully-populated store
#     with the row missing — and so does restoring a release-asset store created
#     before Module 3 existed. Returning NULL there made every anchor skip GREEN
#     against real data, which is the exact false confidence this suite exists
#     to prevent. VERIFIED: a store recording unrelated targets and no
#     sanity_results row produced SKIP 9 / FAIL 0 before this guard.
#
# Skipping is therefore reserved for the one honest case: no pipeline metadata
# at all (a fresh clone, or a local machine with no real-data store — which is
# the expected local outcome, since Modules 1-2 only run in the container).

read_sanity_results <- function() {
  store <- testthat::test_path("..", "..", "_targets")
  if (!file.exists(file.path(store, "meta", "meta"))) {
    return(NULL)
  }
  meta <- targets::tar_meta(store = store)
  if (!"sanity_results" %in% meta$name) {
    stop("the _targets store at ", store, " has pipeline metadata but records ",
         "no `sanity_results` target: the restore predates Module 3, or an ",
         "upstream target failed before it. The credibility anchors must not ",
         "skip against a populated store.")
  }
  err <- meta$error[meta$name == "sanity_results"]
  if (!is.na(err[[1]])) {
    stop("sanity_results is recorded in the _targets store but ERRORED: ",
         err[[1]])
  }
  targets::tar_read(sanity_results, store = store)
}

test_that("ANCHOR: ccRCC driver mutation frequencies match published ranges", {
  sr <- read_sanity_results()
  skip_if(is.null(sr), "sanity_results not in _targets store (run tar_make)")

  mf <- sr$mutation_freq
  expect_true(mf$pass)
  vhl <- mf$per_gene$observed[mf$per_gene$gene == "VHL"]
  expect_gte(vhl, PUBLISHED_MUT_FREQ_RANGES$VHL["low"])
  expect_lte(vhl, PUBLISHED_MUT_FREQ_RANGES$VHL["high"])
})

test_that("ANCHOR: every published-range driver was scored, and each is in range", {
  # `all()` of an empty vector is TRUE, so a panel that silently lost genes
  # would report pass = TRUE having tested nothing. Name the genes that must be
  # present, then re-derive the verdict from the OBSERVED frequencies rather
  # than trusting the flag the function set.
  sr <- read_sanity_results()
  skip_if(is.null(sr), "sanity_results not in _targets store (run tar_make)")

  pg <- sr$mutation_freq$per_gene
  expect_setequal(pg$gene, names(PUBLISHED_MUT_FREQ_RANGES))
  expect_true(all(is.finite(pg$observed)))
  for (g in names(PUBLISHED_MUT_FREQ_RANGES)) {
    rng <- PUBLISHED_MUT_FREQ_RANGES[[g]]
    obs <- pg$observed[pg$gene == g]
    expect_gte(obs, rng[["low"]])
    expect_lte(obs, rng[["high"]])
  }

  # ... and on how many tumours. Published-range membership is a WEAK constraint
  # for a small n, so without the denominator a regression shrinking mut_annot
  # from 417 to a few dozen samples passes every anchor in the suite. The
  # mutation subset is MEASURED at 417 (2016 legacy MAF), and nothing in the DAG
  # subsets it further before this check.
  expect_gte(sr$mutation_freq$n, 0.95 * COHORT_SIZES$mutation_subset)
  expect_lte(sr$mutation_freq$n, COHORT_SIZES$mutation_subset)
})

test_that("ANCHOR: BAP1-mutant tumours have worse OS (HR > 1)", {
  # RE-SPECIFIED ONCE, ON ARITHMETIC — read this before touching it.
  #
  # This anchor used to demand `expect_lt(bs$p_value, 0.05)` (and, below,
  # `expect_gt(bs$ci_low, 1)`). MEASURED on the frozen snapshot, run
  # 30840373033 (printed in docs/results/phase3-anchors-run-30840373033.txt):
  # HR 1.584, CI 0.967-2.595, p = 0.0677, n = 417, 8.63% mutant (36 cases).
  # DERIVED BUT NOT RECORDED ANYWHERE: 138 OS events in the fitted subset, ~12
  # of them in the mutant arm — no transcript prints an event count, and the
  # recorded CI actually implies ~18 in the mutant arm. See the caveat block in
  # R/constants.R; the LEVEL 3 workflow step now prints these, so the next
  # container run settles it.
  # The DIRECTION matches the literature. The significance does not.
  #
  # It cannot. Schoenfeld's requirement at two-sided 0.05 and 80% power for an
  # HR of 1.584 at 8.63% exposure prevalence is ~470 EVENTS; this cohort has
  # 138 (derived) — about 3.4x short. The old requirement was therefore not a test of
  # this pipeline at all, it was a test of the size of TCGA KIRC, and no
  # correct implementation running on this snapshot could ever have satisfied
  # it. That is a MIS-SPECIFIED REQUIREMENT, which is the one and only
  # circumstance in which anything in this suite may be re-specified.
  #
  # WHAT REPLACED IT IS STRICTER, NOT LOOSER. The old assertion was `hr > 1`
  # plus a significance demand; the new one is `hr > 1` PLUS membership of the
  # published effect band PUBLISHED_BAP1_HR_RANGE (1.2-3.0), which the bare
  # `hr > 1` did not constrain at all — an HR of 1.0001 or of 40 both cleared
  # it. p, CI and n are now REPORTED and asserted WELL-FORMED rather than
  # asserted significant, and the under-power is asserted in its own right
  # below so the limitation is a tested claim rather than a footnote.
  #
  # Nothing else in this suite was re-thresholded. SANITY_MAX_PLATFORM_ARI,
  # SANITY_MIN_SILHOUETTE and every published range are untouched, and the
  # m1-m4 anchor is still RED.
  sr <- read_sanity_results()
  skip_if(is.null(sr), "sanity_results not in _targets store (run tar_make)")

  bs <- sr$bap1_survival
  # HARD requirement 1: the published DIRECTION.
  expect_gt(bs$hr, 1)
  # HARD requirement 2: a published-plausible MAGNITUDE. This is the term that
  # catches the failure modes `p < 0.05` was standing in for — a merge blow-up,
  # a fit on a handful of cases, or a coefficient read on the wrong scale all
  # produce an HR far outside the band.
  expect_gte(bs$hr, PUBLISHED_BAP1_HR_RANGE[["low"]])
  expect_lte(bs$hr, PUBLISHED_BAP1_HR_RANGE[["high"]])
  expect_true(bs$pass)

  # REPORTED, not required: the inference must be well-formed, and that is all
  # this cohort licenses anyone to ask of it.
  expect_true(is.finite(bs$p_value))
  expect_gte(bs$p_value, 0)
  expect_lte(bs$p_value, 1)
  expect_true(is.finite(bs$ci_low) && is.finite(bs$ci_high))
  expect_lt(bs$ci_low, bs$hr)
  expect_gt(bs$ci_high, bs$hr)
  expect_gt(bs$ci_low, 0)
  expect_true(is.finite(bs$n) && bs$n > 0)
})

test_that("ANCHOR: the BAP1 control is UNDERPOWERED, and that is asserted", {
  # The limitation, tested rather than written down. If this ever fails it means
  # EITHER the event count rose OR the observed effect size grew until this
  # cohort could detect it — the comment block in the anchor above (and in
  # R/constants.R) has then become false and must be rewritten, which is the
  # point of asserting it. It does NOT mean "the cohort grew": the observed
  # hazard ratio enters the arithmetic too, and it can move on a byte-identical
  # cohort.
  #
  # This is a DESIGN-ADEQUACY statement: how many deaths an effect this size
  # needs before alpha 0.05 is reachable. It is NOT post-hoc power, and it is
  # NOT an argument that BAP1 has no effect. Nothing here interprets p = 0.0677.
  sr <- read_sanity_results()
  skip_if(is.null(sr), "sanity_results not in _targets store (run tar_make)")

  bs <- sr$bap1_survival

  # The pipeline's own arithmetic, recomputed here from the reported inputs so
  # the assertion does not simply trust a field the function set.
  expected_required <- fn_schoenfeld_events(bs$hr, bs$mutant_frac)
  expect_equal(bs$events_required, expected_required)

  # ~470 events needed against the events this fit actually saw.
  expect_lt(bs$n_events, bs$events_required)
  expect_true(bs$underpowered)

  # Not marginally short — stated as a property of the DESIGN, not of today's
  # point estimate.
  #
  # `expect_gt(bs$events_required, 2 * bs$n_events)` USED TO BE HERE and has
  # been REMOVED as MIS-SPECIFIED. `events_required` is computed from the
  # OBSERVED hr, so with n_events and mutant_frac held fixed the assertion flips
  # at HR = 1.8229 (COMPUTED here by uniroot on fn_schoenfeld_events) — only 15%
  # above the measured 1.584 and well inside PUBLISHED_BAP1_HR_RANGE (1.2-3.0),
  # which the anchor above accepts. For any HR in [1.823, 2.338) the study is
  # still genuinely under-powered — `bs$underpowered` only flips at HR = 2.3377
  # — yet the old line went RED purely because the shortfall was a factor of 1.9
  # rather than 2. That is a false red on a result this same suite certifies as
  # correct, and pinning an incidental magnitude of an estimate is exactly what
  # an anchor must not do.
  #
  # What replaces it says the same thing without depending on the estimate: the
  # SMALLEST hazard ratio this design could detect at SURVIVAL_TARGET_POWER is
  # larger than the WEAKEST effect the literature reports. That is a function of
  # the exposure prevalence and the event count only — the two things that
  # actually determine whether the cohort can speak — so it stays true wherever
  # the point estimate lands inside the band it is allowed to occupy. MEASURED:
  # minimum detectable HR 2.338 against a published floor of 1.2.
  mdhr <- stats::uniroot(
    function(h) fn_schoenfeld_events(h, bs$mutant_frac) - bs$n_events,
    c(1.001, 50)
  )$root
  expect_gt(mdhr, PUBLISHED_BAP1_HR_RANGE[["low"]])

  # The composition the arithmetic rests on. MEASURED: 8.63% mutant (36 of
  # 417). DERIVED but NOT YET RECORDED in any committed transcript: 138 OS
  # events, ~12 in the mutant arm (see R/constants.R). A regression that changed
  # the exposure prevalence would change the event requirement without touching
  # anything else in this file.
  expect_gt(bs$mutant_frac, 0)
  expect_lt(bs$mutant_frac, 1)
  expect_equal(bs$mutant_frac, bs$n_mutant / bs$n)
  expect_lte(bs$n_events_mutant, bs$n_events)
  expect_lte(bs$n_mutant, bs$n)
  # The mutant arm carries the information, and it is thin. Named explicitly so
  # nobody reads the 138-event total as if it were the effective sample size.
  expect_lt(bs$n_events_mutant, bs$n_events / 2)

  # ... AND IT MUST NOT BE EMPTY. `expect_lt(n_events_mutant, n_events / 2)` is
  # an UPPER bound only, and `mutant_frac` was bounded solely by `> 0`. Together
  # with `hr > 1` and the published HR band, that admitted a fit on a COLLAPSED
  # exposed arm: VERIFIED by replaying these anchor bodies, a stub with
  # hr 1.6, n_mutant 5, n_events_mutant 2, CI 0.35-7.3, p 0.51 — an
  # uninformative fit of exactly the class `expect_gt(ci_low, 1)` used to catch
  # before it was removed — passed all three BAP1 anchors. The two floors below
  # restore that guarantee WITHOUT restoring a significance demand: neither is a
  # statement about p, and neither can turn a red check green.
  #
  # (i) The exposed arm must be the BAP1 arm the suite already anchors. The
  # mutation_freq anchor independently certifies BAP1 at 8.63% against the
  # published 6-18% band; the survival fit's own exposure prevalence must agree
  # with it, or the two halves of the suite are describing different cohorts.
  # This is the same kind of published-literature term as PUBLISHED_BAP1_HR_RANGE.
  expect_gte(bs$mutant_frac, PUBLISHED_MUT_FREQ_RANGES$BAP1[["low"]])
  expect_lte(bs$mutant_frac, PUBLISHED_MUT_FREQ_RANGES$BAP1[["high"]])
  # (ii) The informative arm must clear the repo's own events-per-variable rule.
  # A one-covariate Cox needs EPV_CAP events in the arm carrying the contrast;
  # below that the estimate is noise regardless of where the point estimate
  # lands. EPV_CAP, not a fresh literal, so this floor moves with the model
  # budget the rest of the pipeline is designed against.
  expect_gte(bs$n_events_mutant, EPV_CAP)
})

test_that("ANCHOR: the BAP1 survival anchor was fitted on the mutation subset", {
  # An inflated HR from a handful of cases, or a merge that fanned out into a
  # cartesian product, would both satisfy "HR > 1" while meaning nothing. The
  # join cannot legitimately exceed the mutation-annotated subset (n = 417).
  sr <- read_sanity_results()
  skip_if(is.null(sr), "sanity_results not in _targets store (run tar_make)")

  bs <- sr$bap1_survival
  # The floor must be the MEASURED cohort, not a Module 2 guard. It used to be
  # MIN_MUT_ANNOT_SAMPLES = 50L — a constant with a different meaning ("the
  # mutation subset must overlap the factors"), 8x below the truth. VERIFIED by
  # replaying these anchors against a stubbed sanity_results: an HR of 2.1 fitted
  # on n = 60 (an ID-harmonisation or merge regression dropping 85% of cases)
  # passed ALL anchors, FAIL 0 — exactly the "inflated HR from a handful of
  # cases" this test says it exists to exclude.
  #
  # MEASURED: the mutation subset is 417; `clinical` covers all 536 colData
  # cases, so the inner join loses only the handful with missing or non-positive
  # os_time (2 of 524 in the main cohort).
  expect_gte(bs$n, 0.95 * COHORT_SIZES$mutation_subset)
  expect_lte(bs$n, COHORT_SIZES$mutation_subset)
  # The event count too. `n` alone does not bound what the fit saw: a cohort of
  # 417 with 6 deaths carries no more information than a cohort of 20, and the
  # Cox likelihood is driven by events, not cases. DERIVED, NOT RECORDED: 138
  # OS events in the mutation subset (see R/constants.R). The floor below is
  # MIN_OS_EVENTS, a constant, precisely so this anchor does not depend on that
  # unrecorded figure.
  expect_gte(bs$n_events, MIN_OS_EVENTS)
  expect_lte(bs$n_events, bs$n)

  # `expect_gt(bs$ci_low, 1)` USED TO BE HERE and has been REMOVED, for exactly
  # the reason set out in the preceding anchor: with 138 events (derived, not
  # recorded) against the ~470
  # Schoenfeld requires, a 95% interval excluding 1 is arithmetically out of
  # this cohort's reach (MEASURED CI 0.967-2.595). Requiring it made the anchor
  # a test of TCGA KIRC's size rather than of this pipeline. It is replaced by
  # WELL-FORMEDNESS — the interval must be finite, ordered, positive, and must
  # bracket the point estimate — which is what still catches a broken fit. The
  # thing the old assertion was really guarding (an implausible HR) is now
  # caught by the PUBLISHED_BAP1_HR_RANGE band, which is strictly tighter than
  # the `hr > 1` it sat beside. Do not reinstate it; see the under-power anchor.
  expect_true(is.finite(bs$ci_low) && is.finite(bs$ci_high))
  expect_gt(bs$ci_low, 0)
  expect_lt(bs$ci_low, bs$ci_high)
  expect_lt(bs$ci_low, bs$hr)
  expect_gt(bs$ci_high, bs$hr)
  expect_true(bs$pass)

  # WELL-FORMEDNESS IS NOT INFORMATIVENESS, and the removal above left only the
  # former. An interval can be finite, ordered, positive and bracket the point
  # estimate while being 21x wide — which is what a fit on five mutant cases
  # returns, and which the old `ci_low > 1` implicitly excluded. BAP1_MAX_CI_RATIO
  # restores that exclusion as a PRECISION bound rather than as significance:
  # ci_high / ci_low = exp(2 * 1.96 * SE) does not depend on where the interval
  # sits relative to 1, so no value of it can make a failing check pass.
  # MEASURED: 2.595 / 0.967 = 2.68 against a ceiling of 5.
  expect_lt(bs$ci_high / bs$ci_low, BAP1_MAX_CI_RATIO)
})

test_that("ANCHOR: methylation recovers four strata (m1-m4) -- RED, a real negative result", {
  # THIS ANCHOR IS EXPECTED TO FAIL ON THE FROZEN SNAPSHOT, AND MUST KEEP
  # FAILING. It is not broken and it is not waiting to be fixed.
  #
  # WHAT WAS MEASURED (run 30840373033): silhouette 0.1197, Kruskal p 1.3e-82,
  # 4 non-empty groups over 4568 complete CpGs and 524 cases — and platform_ari
  # 0.583 against the 0.25 ceiling, platform_p 3.0e-113. The partition tracks
  # the HM27/HM450 assay, not the published m1-m4 biology.
  #
  # WHAT THE FOLLOW-UP DIAGNOSTIC ADDED (run 30911448546): re-running the same
  # 4-means INSIDE each platform gives HM27 0.0858 (n = 214, 4658 CpGs) and
  # HM450 0.0489 (n = 310, 4905 CpGs). Both are BELOW the merged 0.1197. The
  # merged "structure" is substantially the assay split itself; there is no
  # m1-m4 stratification hiding underneath it on this snapshot.
  #
  # WHAT WAS DECIDED (and is NOT to be re-litigated here): keep all 524 cases,
  # do NOT restrict to one platform, do NOT ComBat — too few cases overlap the
  # two assays and the probe sets differ, so a correction could be neither
  # validated nor trusted not to erase real signal. Platform is adjusted for as
  # a covariate downstream instead. (The overlap count has been stated as 3 but
  # is not yet in any committed transcript; the methyl_platform_overlap target
  # computes it and the CI transcript will record it.)
  #
  # DO NOT turn this green. Raising SANITY_MAX_PLATFORM_ARI or lowering
  # SANITY_MIN_SILHOUETTE would convert a documented negative result into a
  # fabricated positive one, which is the single failure mode this whole phase
  # exists to prevent. The correct response to a red light here is to report it.
  sr <- read_sanity_results()
  skip_if(is.null(sr), "sanity_results not in _targets store (run tar_make)")

  ms <- sr$methyl_strata
  expect_identical(ms$n_strata, 4L)
  expect_true(ms$pass,
              info = paste("EXPECTED RED on the frozen snapshot --", ms$message))
})

test_that("ANCHOR: the m1-m4 red state carries its within-platform evidence", {
  # A bare FALSE invites a future reader to "fix" it by moving a threshold. The
  # verdict must therefore travel with the measurement that explains it, so this
  # asserts the DIAGNOSTIC EVIDENCE is present and well-formed. It deliberately
  # does NOT assert which way the comparison came out — that is the finding, and
  # a finding that a test forces to stay the same is not a finding.
  sr <- read_sanity_results()
  skip_if(is.null(sr), "sanity_results not in _targets store (run tar_make)")

  ms <- sr$methyl_strata
  expect_type(ms$message, "character")
  expect_true(nzchar(ms$message))

  # One row per assay, each naming the cohort it was computed on.
  expect_s3_class(ms$within_platform, "data.frame")
  expect_setequal(ms$within_platform$platform, METHYL_PLATFORMS)
  expect_true(all(is.finite(ms$within_platform$silhouette)))
  expect_true(all(ms$within_platform$n > METHYL_N_STRATA))
  # Both arms together are the whole cohort — no case is scored twice or lost.
  expect_identical(sum(ms$within_platform$n), length(ms$cluster))

  # The signature of an assay-driven partition, recorded as a boolean so it can
  # be read off the object without recomputing anything.
  expect_type(ms$merged_exceeds_within, "logical")
  expect_false(is.na(ms$merged_exceeds_within))
  expect_identical(
    ms$merged_exceeds_within,
    ms$silhouette > max(ms$within_platform$silhouette)
  )
})

test_that("ANCHOR: the m1-m4 strata are separated, not just k-means returning k", {
  # n_strata == k is tautological for Hartigan-Wong k-means, so the stratum
  # anchor is only evidence through the silhouette; assert the discriminating
  # quantity directly instead of the aggregate flag.
  sr <- read_sanity_results()
  skip_if(is.null(sr), "sanity_results not in _targets store (run tar_make)")

  ms <- sr$methyl_strata
  # NOTE, so this passing line is not misread as support for m1-m4: the merged
  # silhouette does clear the floor (0.1197 vs 0.10), but the within-platform
  # diagnostic shows it clears it BECAUSE of the assay split — HM27 alone gives
  # 0.0858 and HM450 alone 0.0489, both below the merged value. Silhouette
  # measures separation, not what is being separated. The platform_ari anchor
  # below is the one that names what this partition actually found.
  expect_gt(ms$silhouette, SANITY_MIN_SILHOUETTE)
  expect_lt(ms$kw_p_value, SANITY_MAX_P)
  # Assert the MEASURED completeness, not the function's own internal floor.
  # `expect_gte(ms$n_cpg_used, SANITY_MIN_COMPLETE_FRAC * N_TOP_CPGS)` was
  # tautological: fn_complete_cpgs already stop()s below that fraction, and
  # nrow(methyl_mat) is exactly N_TOP_CPGS, so if the check ran at all the
  # assertion could not fail. VERIFIED: a stub with n_cpg_used = 2500 — half the
  # CpGs lost versus the measured 4568 — passed every anchor, FAIL 0.
  #
  # methyl_mat is 91.4% complete on the frozen snapshot (432 of 5000 CpGs carry
  # a non-finite value), so a large drop is a data regression, not tolerable
  # loss.
  expect_gte(ms$n_cpg_used, 0.85 * N_TOP_CPGS)
  expect_lte(ms$n_cpg_used, N_TOP_CPGS)
  # Every case in the main cohort gets a stratum; nothing is dropped, and the
  # check ran on the real cohort rather than a fixture.
  expect_gte(length(ms$cluster), COHORT_MIN)
  expect_lte(length(ms$cluster), COHORT_MAX)
})

test_that("ANCHOR: the m1-m4 strata are biology, not the HM27/HM450 batch -- RED", {
  # methyl_mat is cbind(HM27, HM450) with NO batch correction, and platform is
  # the strongest single axis in merged 27k/450k M-values. MEASURED on
  # constructed data carrying ONLY a per-platform offset and no biological
  # strata at all, this check reported pass = TRUE from 1.5 SD upward. So the
  # platform term must be present AND satisfied — an NA here means the DAG
  # stopped supplying methyl_platform, which is a silent loss of the guard.
  #
  # ON THE FROZEN SNAPSHOT THIS FAILS, AND MUST: platform_ari 0.583 against the
  # 0.25 ceiling (run 30840373033). It is the term that puts the m1-m4 verdict
  # in the red, and it is doing precisely the job it was written to do. The
  # ceiling was calibrated in R/constants.R against a cleanly bimodal pair of
  # regimes (-0.003 when k-means recovered true strata, 0.532 when it locked
  # onto the assay); 0.583 sits in the second regime. Do not raise it.
  sr <- read_sanity_results()
  skip_if(is.null(sr), "sanity_results not in _targets store (run tar_make)")

  ms <- sr$methyl_strata
  expect_true(is.finite(ms$platform_ari))
  expect_lt(
    ms$platform_ari, SANITY_MAX_PLATFORM_ARI,
    label = paste("cluster-vs-platform ARI (EXPECTED RED: the merged",
                  "methylation partition tracks assay platform, so the m1-m4",
                  "strata are not recovered on this snapshot)")
  )
})

test_that("ANCHOR: ccA/ccB expression signatures separate", {
  sr <- read_sanity_results()
  skip_if(is.null(sr), "sanity_results not in _targets store (run tar_make)")

  expect_true(sr$ccab_signature$pass)
})

test_that("ANCHOR: proliferation and angiogenesis are OPPOSED (Brannon 2010)", {
  # The non-circular half of the ccA/ccB check: the two programmes must run in
  # opposite directions across tumours. The silhouette and Wilcoxon terms are
  # computed from the same k-means split they then evaluate, so on their own
  # they greenlight structureless data (measured: 3 of 5 noise matrices).
  sr <- read_sanity_results()
  skip_if(is.null(sr), "sanity_results not in _targets store (run tar_make)")

  cs <- sr$ccab_signature
  expect_lt(cs$anticorr_rho, 0)
  expect_lt(cs$anticorr_p_value, SANITY_MAX_P)
  # The 2-D constant: 0.10 is the methylation floor and pure 2-D noise clears it
  # every time, so asserting it here would be asserting nothing.
  expect_gt(cs$silhouette, SANITY_MIN_SILHOUETTE_2D)
  expect_gte(length(cs$axis_score), COHORT_MIN)
  expect_lte(length(cs$axis_score), COHORT_MAX)
})

test_that("ANCHOR: the ccA/ccB axis was scored from the FULL published panels", {
  # A truncated panel returns an object indistinguishable from the full run, so
  # this must be asserted rather than assumed. The DAG feeds this check
  # `rna_full` (before fn_top_variable), precisely so the published panels are
  # never subject to a data-driven feature filter — if a marker is missing here,
  # it is genuinely absent from the RSEM gene set, and that is a FINDING about
  # the data, not a reason to lower the floor.
  sr <- read_sanity_results()
  skip_if(is.null(sr), "sanity_results not in _targets store (run tar_make)")

  cs <- sr$ccab_signature
  expect_setequal(cs$markers_used$ccb, CCB_PROLIFERATION_MARKERS)
  expect_setequal(cs$markers_used$cca, CCA_ANGIOGENESIS_MARKERS)
  expect_identical(cs$n_ccb_used, length(CCB_PROLIFERATION_MARKERS))
  expect_identical(cs$n_cca_used, length(CCA_ANGIOGENESIS_MARKERS))
  expect_gte(cs$n_ccb_used, SANITY_MIN_MARKERS_PER_PANEL)
  expect_gte(cs$n_cca_used, SANITY_MIN_MARKERS_PER_PANEL)
})

test_that("ANCHOR: the MOFA subtypes stay INDEPENDENT of the assay platform", {
  # The one clean result in the platform diagnostic, pinned so it cannot break
  # silently. MEASURED (run 30911448546): ARI(subtypes_mofa, platform) = 0.0058
  # over 524 cases (HM27 214 / HM450 310), cross-tab S1 11/9, S2 120/186,
  # S3 33/43, S4 50/72.
  #
  # This is NOT a foregone conclusion. Individual MOFA factors are heavily
  # platform-loaded on this snapshot — Factor2 separates the assays at AUC
  # 0.888 (q 2.9e-50 corrected; the transcript's 1.4e-50 was depressed by an
  # impossible Factor6 p = 0 — see the annotation on that transcript), Factor5 at
  # 0.818, Factor6 at 0.735 — so the material for
  # a platform-driven subtype assignment is present and the 4-means over the
  # 15-factor space simply does not use it. Contrast the methylation k-means,
  # which does (ARI 0.583) and whose anchor is red.
  #
  # SUBTYPE_MAX_PLATFORM_ARI can only turn a green verdict red; never raise it.
  sr <- read_sanity_results()
  skip_if(is.null(sr), "sanity_results not in _targets store (run tar_make)")

  sp <- sr$subtype_platform
  expect_true(is.finite(sp$ari))
  expect_lt(sp$ari, SUBTYPE_MAX_PLATFORM_ARI)
  expect_true(sp$pass)

  # Scored on the WHOLE cohort and on the real subtype set — an ARI near zero
  # computed from a handful of cases, or from a degenerate one-level labelling,
  # would be free and would mean nothing.
  expect_gte(sp$n, COHORT_MIN)
  expect_lte(sp$n, COHORT_MAX)
  expect_identical(sp$n_subtypes, K_SUBTYPES)
  expect_identical(dim(sp$cross_tab),
                   c(K_SUBTYPES, length(METHYL_PLATFORMS)))
  # Every case falls in exactly one (subtype, platform) cell.
  expect_identical(sum(sp$cross_tab), sp$n)
  # Both assays are represented in the comparison at all.
  expect_true(all(colSums(sp$cross_tab) > 0))
})

test_that("ANCHOR: sanity_results carries every check with a real verdict", {
  # A missing element would make `sr$<name>$pass` NULL, and expect_true(NULL)
  # errors rather than passing — but the failure would name the wrong thing.
  # Assert the shape of the credibility anchor itself.
  sr <- read_sanity_results()
  skip_if(is.null(sr), "sanity_results not in _targets store (run tar_make)")

  expect_setequal(names(sr), c("mutation_freq", "bap1_survival",
                               "methyl_strata", "ccab_signature",
                               "subtype_platform"))
  for (nm in names(sr)) {
    expect_type(sr[[nm]]$label, "character")
    expect_type(sr[[nm]]$pass, "logical")
    expect_length(sr[[nm]]$pass, 1L)
    expect_false(is.na(sr[[nm]]$pass))
  }
})
