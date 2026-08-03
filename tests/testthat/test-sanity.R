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
  sr <- read_sanity_results()
  skip_if(is.null(sr), "sanity_results not in _targets store (run tar_make)")

  bs <- sr$bap1_survival
  expect_gt(bs$hr, 1)
  expect_lt(bs$p_value, 0.05)
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
  # The interval must exclude the null in the HARMFUL direction, not merely
  # differ from it.
  expect_true(is.finite(bs$ci_low) && is.finite(bs$ci_high))
  expect_gt(bs$ci_low, 1)
  expect_true(bs$pass)
})

test_that("ANCHOR: methylation recovers four strata (m1-m4)", {
  sr <- read_sanity_results()
  skip_if(is.null(sr), "sanity_results not in _targets store (run tar_make)")

  ms <- sr$methyl_strata
  expect_identical(ms$n_strata, 4L)
  expect_true(ms$pass)
})

test_that("ANCHOR: the m1-m4 strata are separated, not just k-means returning k", {
  # n_strata == k is tautological for Hartigan-Wong k-means, so the stratum
  # anchor is only evidence through the silhouette; assert the discriminating
  # quantity directly instead of the aggregate flag.
  sr <- read_sanity_results()
  skip_if(is.null(sr), "sanity_results not in _targets store (run tar_make)")

  ms <- sr$methyl_strata
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

test_that("ANCHOR: the m1-m4 strata are biology, not the HM27/HM450 batch", {
  # methyl_mat is cbind(HM27, HM450) with NO batch correction, and platform is
  # the strongest single axis in merged 27k/450k M-values. MEASURED on
  # constructed data carrying ONLY a per-platform offset and no biological
  # strata at all, this check reported pass = TRUE from 1.5 SD upward. So the
  # platform term must be present AND satisfied — an NA here means the DAG
  # stopped supplying methyl_platform, which is a silent loss of the guard.
  sr <- read_sanity_results()
  skip_if(is.null(sr), "sanity_results not in _targets store (run tar_make)")

  ms <- sr$methyl_strata
  expect_true(is.finite(ms$platform_ari))
  expect_lt(ms$platform_ari, SANITY_MAX_PLATFORM_ARI)
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

test_that("ANCHOR: sanity_results carries all four checks with real verdicts", {
  # A missing element would make `sr$<name>$pass` NULL, and expect_true(NULL)
  # errors rather than passing — but the failure would name the wrong thing.
  # Assert the shape of the credibility anchor itself.
  sr <- read_sanity_results()
  skip_if(is.null(sr), "sanity_results not in _targets store (run tar_make)")

  expect_setequal(names(sr), c("mutation_freq", "bap1_survival",
                               "methyl_strata", "ccab_signature"))
  for (nm in names(sr)) {
    expect_type(sr[[nm]]$label, "character")
    expect_type(sr[[nm]]$pass, "logical")
    expect_length(sr[[nm]]$pass, 1L)
    expect_false(is.na(sr[[nm]]$pass))
  }
})
