# The platform covariate: fn_attach_platform (R/functions_clinical.R) and the
# `clinical` target that carries it.
#
# The HM27/HM450 split is NOT corrected anywhere in this DAG and — per run
# 30911448546 — is a large, measurable axis of the merged methylation data
# (cohort 214 HM27 / 310 HM450; MOFA Factor2 separates the platforms at
# AUC 0.888; the merged 4-means silhouette 0.1197 exceeds either arm alone,
# 0.0858 / 0.0489). Phase 4 therefore ADJUSTS for platform instead of
# correcting it, and that requires platform to reach the survival model
# attached to the same frame as the outcome. These tests hold that plumbing to
# its contract; none of them touches a published threshold, and none of them
# can turn the still-red m1-m4 anchor green.

# --- Fixtures --------------------------------------------------------------

make_clinical_frame <- function(ids) {
  data.frame(
    sample_id = ids,
    os_time   = seq_along(ids) * 100,
    os_event  = rep(c(0L, 1L), length.out = length(ids)),
    stringsAsFactors = FALSE
  )
}

make_platform_vector <- function(ids, values) {
  p <- factor(values, levels = METHYL_PLATFORMS)
  names(p) <- ids
  p
}

# --- Unit tests: fn_attach_platform ----------------------------------------

test_that("fn_attach_platform appends platform and leaves the survival columns alone", {
  # Arrange
  ids <- paste0("TCGA-AA-000", 1:4)
  clinical <- make_clinical_frame(ids)
  platform <- make_platform_vector(ids, c("HM27", "HM450", "HM450", "HM27"))

  # Act
  out <- fn_attach_platform(clinical, platform)

  # Assert — Module 4 depends on these three names, in this order, unchanged.
  expect_identical(names(out), c("sample_id", "os_time", "os_event", "platform"))
  expect_identical(out$sample_id, clinical$sample_id)
  expect_identical(out$os_time, clinical$os_time)
  expect_identical(out$os_event, clinical$os_event)
  expect_identical(nrow(out), nrow(clinical))

  expect_true(is.factor(out$platform))
  expect_identical(levels(out$platform), METHYL_PLATFORMS)
  expect_false(anyNA(out$platform))
})

test_that("fn_attach_platform joins by ID, not by row position", {
  # Arrange — the platform vector is in a DIFFERENT order from the frame, which
  # is the realistic case (clinical follows colData, methyl_platform follows
  # common_ids). A positional cbind would mislabel every case here.
  ids <- paste0("TCGA-AA-000", 1:4)
  clinical <- make_clinical_frame(ids)
  shuffled <- ids[c(3L, 1L, 4L, 2L)]
  platform <- make_platform_vector(shuffled, c("HM450", "HM27", "HM27", "HM450"))

  # Act
  out <- fn_attach_platform(clinical, platform)

  # Assert
  expect_identical(as.character(out$platform),
                   as.character(platform[out$sample_id]))
  expect_identical(as.character(out$platform),
                   c("HM27", "HM450", "HM450", "HM27"))
})

test_that("cases outside the platform's cohort are NA, and cohort cases never are", {
  # Arrange — `clinical` covers all of colData (536 on the snapshot) while
  # methyl_platform covers the 524-case main cohort, so this asymmetry is the
  # real one. Inventing HM450 for an uncovered case would be a fabrication.
  ids <- paste0("TCGA-AA-000", 1:5)
  clinical <- make_clinical_frame(ids)
  cohort <- ids[1:3]
  platform <- make_platform_vector(cohort, c("HM27", "HM450", "HM27"))

  # Act
  out <- fn_attach_platform(clinical, platform)

  # Assert
  expect_false(anyNA(out$platform[out$sample_id %in% cohort]))
  expect_true(all(is.na(out$platform[!out$sample_id %in% cohort])))
  expect_identical(sum(is.na(out$platform)), 2L)
})

test_that("fn_attach_platform refuses a cohort case missing from clinical", {
  # A platform entry with no clinical row means the join silently loses a case
  # the model was supposed to adjust; that must stop, not shrink the cohort.
  ids <- paste0("TCGA-AA-000", 1:3)
  clinical <- make_clinical_frame(ids)
  platform <- make_platform_vector(c(ids, "TCGA-AA-0009"),
                                   c("HM27", "HM450", "HM27", "HM450"))

  expect_error(fn_attach_platform(clinical, platform), "absent from clinical")
})

test_that("fn_attach_platform refuses malformed platform vectors", {
  ids <- paste0("TCGA-AA-000", 1:3)
  clinical <- make_clinical_frame(ids)
  good <- make_platform_vector(ids, c("HM27", "HM450", "HM27"))

  # Not a factor / wrong levels: a character vector or a re-levelled factor
  # would silently become a differently-coded covariate downstream.
  expect_error(fn_attach_platform(clinical, as.character(good)),
               "must be a factor with levels")
  expect_error(
    fn_attach_platform(clinical, factor(as.character(good),
                                        levels = c("HM450", "HM27", "EPIC"))),
    "must be a factor with levels"
  )

  # NA inside the cohort it claims to cover.
  na_platform <- good
  na_platform[2L] <- NA
  expect_error(fn_attach_platform(clinical, na_platform), "must be known")

  # Unnamed / duplicated names: the join key must be unique.
  unnamed <- good
  names(unnamed) <- NULL
  expect_error(fn_attach_platform(clinical, unnamed), "unique, non-missing")
  dup <- make_platform_vector(c(ids[1], ids[1], ids[2]),
                              c("HM27", "HM27", "HM450"))
  expect_error(fn_attach_platform(clinical, dup), "unique, non-missing")
})

test_that("fn_attach_platform refuses a malformed clinical frame", {
  ids <- paste0("TCGA-AA-000", 1:3)
  clinical <- make_clinical_frame(ids)
  platform <- make_platform_vector(ids, c("HM27", "HM450", "HM27"))

  expect_error(fn_attach_platform(clinical[, c("sample_id", "os_time")], platform),
               "lacks required survival columns")
  expect_error(fn_attach_platform(as.list(clinical), platform),
               "needs a clinical data.frame")

  # Duplicate case IDs would fan the join out; and re-attaching would let two
  # different platform codings coexist under one name.
  dup_clinical <- make_clinical_frame(c(ids[1], ids))
  expect_error(fn_attach_platform(dup_clinical, platform), "not unique")
  expect_error(fn_attach_platform(fn_attach_platform(clinical, platform), platform),
               "already carries")
})

test_that("fn_attach_platform mutates neither argument", {
  ids <- paste0("TCGA-AA-000", 1:3)
  clinical <- make_clinical_frame(ids)
  platform <- make_platform_vector(ids, c("HM27", "HM450", "HM27"))
  clinical_before <- clinical
  platform_before <- platform

  invisible(fn_attach_platform(clinical, platform))

  expect_identical(clinical, clinical_before)
  expect_identical(platform, platform_before)
})

# --- Unit tests: fn_factor_platform_association -----------------------------

make_factor_matrix <- function(ids, platform, loaded = character(0),
                               shift = 3, n_factors = 6L, seed = 11L) {
  set.seed(seed)
  m <- matrix(stats::rnorm(length(ids) * n_factors), nrow = length(ids),
              dimnames = list(ids, paste0("Factor", seq_len(n_factors))))
  # A "loaded" factor is shifted by platform: that is what a confounded axis is.
  for (f in loaded) {
    m[, f] <- m[, f] + ifelse(platform[ids] == "HM27", 0, shift)
  }
  m
}

test_that("fn_factor_platform_association separates loaded from clean factors", {
  # Arrange — Factor2 is a pure assay effect, everything else is iid noise.
  ids <- paste0("P", seq_len(300L))
  platform <- make_platform_vector(ids, rep(c("HM27", "HM450"), c(120L, 180L)))
  factors <- make_factor_matrix(ids, platform, loaded = "Factor2")

  # Act
  res <- fn_factor_platform_association(factors, platform)

  # Assert — shape first: one row per factor, in matrix column order.
  expect_identical(res$factor, colnames(factors))
  expect_true(all(res$n == length(ids)))
  # AUC is a separation measure, so it never drops below 0.5 regardless of
  # which platform is the first level.
  expect_true(all(res$auc >= 0.5 & res$auc <= 1))
  expect_true(all(res$q_value >= res$p_value - 1e-12))

  loaded <- res[res$factor == "Factor2", ]
  clean <- res[res$factor != "Factor2", ]
  expect_gt(loaded$auc, 0.9)
  expect_lt(loaded$q_value, SANITY_MAX_P)
  expect_true(all(clean$q_value > SANITY_MAX_P))

  # Nothing here should be degenerate: the normal approximation has a strictly
  # positive p floor at these group sizes, so an exact zero means its variance
  # collapsed. Run 30911448546 recorded one (Factor6, p = q = 0 at AUC 0.735),
  # which is why the flag exists.
  expect_false(any(res$degenerate))
  expect_true(all(res$p_value > 0))
})

test_that("fn_factor_platform_association is invariant to platform level order", {
  # The AUC must describe the separation, not which assay happens to be first.
  ids <- paste0("P", seq_len(200L))
  platform <- make_platform_vector(ids, rep(c("HM27", "HM450"), c(80L, 120L)))
  factors <- make_factor_matrix(ids, platform, loaded = "Factor3")
  flipped <- factor(as.character(platform), levels = rev(METHYL_PLATFORMS))
  names(flipped) <- names(platform)

  expect_equal(fn_factor_platform_association(factors, platform)$auc,
               fn_factor_platform_association(factors, flipped)$auc)
})

test_that("fn_factor_platform_association refuses inputs it cannot score", {
  ids <- paste0("P", seq_len(60L))
  platform <- make_platform_vector(ids, rep(c("HM27", "HM450"), c(25L, 35L)))
  factors <- make_factor_matrix(ids, platform)

  unnamed <- factors
  rownames(unnamed) <- NULL
  expect_error(fn_factor_platform_association(unnamed, platform),
               "unique rownames")
  expect_error(fn_factor_platform_association(as.data.frame(factors), platform),
               "numeric samples x factors matrix")
  expect_error(
    fn_factor_platform_association(factors, stats::setNames(platform, NULL)),
    "must be NAMED"
  )
  expect_error(fn_factor_platform_association(factors, platform[1:10]),
               "carry no platform label")
  one_level <- make_platform_vector(ids, rep("HM27", length(ids)))
  expect_error(fn_factor_platform_association(factors, one_level),
               "exactly two assay levels")
  # Non-finite scores would silently shrink the comparison behind the AUC.
  bad <- factors
  bad[3L, 2L] <- NA_real_
  expect_error(fn_factor_platform_association(bad, platform), "non-finite")
})

# --- Anchors on the real `clinical` target ---------------------------------
# Same discipline as read_sanity_results() in test-sanity.R. The reader itself,
# `read_pipeline_target()`, now lives in helper-fixtures.R so the Module-4
# anchors in test-survival.R can use the same one honest skip rule: a store
# with no pipeline metadata at all skips; a populated store missing the target,
# or recording it as errored, FAILS.

test_that("ANCHOR: the clinical target carries a platform factor for the whole cohort", {
  clinical <- read_pipeline_target("clinical")
  skip_if(is.null(clinical), "clinical not in _targets store (run tar_make)")
  platform <- read_pipeline_target("methyl_platform")

  # The Module 4 contract first: the survival columns are unchanged.
  expect_true(all(c("sample_id", "os_time", "os_event") %in% names(clinical)))
  expect_true("platform" %in% names(clinical))
  expect_true(is.factor(clinical$platform))
  expect_identical(levels(clinical$platform), METHYL_PLATFORMS)
  expect_identical(levels(clinical$platform), c("HM27", "HM450"))

  # Every case in the main cohort has a known platform. NA here would mean the
  # survival model silently drops cases (coxph deletes incomplete rows), so
  # this is the property Phase 4 actually depends on.
  #
  # The cohort is `common_ids` — Module 1's own definition — NOT
  # names(methyl_platform). Taking the cohort from the platform vector would
  # make coverage self-satisfying: a methyl_platform that lost 24 cases would
  # define them out of the cohort and the assertion could not fail.
  cohort_ids <- read_pipeline_target("common_ids")
  expect_setequal(names(platform), cohort_ids)
  in_cohort <- clinical$platform[match(cohort_ids, clinical$sample_id)]
  expect_identical(length(in_cohort), length(cohort_ids))
  expect_false(anyNA(in_cohort))
  expect_identical(as.character(in_cohort),
                   as.character(platform[cohort_ids]))

  # ... on the real cohort, not a fixture, and with both arms genuinely present:
  # a covariate that is constant carries no adjustment at all.
  expect_gte(length(cohort_ids), COHORT_MIN)
  expect_lte(length(cohort_ids), COHORT_MAX)
  n_hm27  <- sum(in_cohort == "HM27", na.rm = TRUE)
  n_hm450 <- sum(in_cohort == "HM450", na.rm = TRUE)
  expect_gt(n_hm27, 0L)
  expect_gt(n_hm450, 0L)
  expect_identical(n_hm27 + n_hm450, length(cohort_ids))
})

test_that("ANCHOR: the clinical platform column reproduces the measured HM27/HM450 split", {
  # MEASURED, NOT ASSUMED: GitHub Actions run 30911448546 on the frozen
  # curatedTCGAData KIRC snapshot — 214 HM27 / 310 HM450 across the 524-case
  # cohort. Agreement with methyl_platform (asserted above) cannot catch a
  # regression INSIDE methyl_platform's own rule, e.g. an inverted ifelse that
  # would relabel every case while keeping the two targets consistent.
  clinical <- read_pipeline_target("clinical")
  skip_if(is.null(clinical), "clinical not in _targets store (run tar_make)")

  expect_identical(sum(clinical$platform == "HM27", na.rm = TRUE), 214L)
  expect_identical(sum(clinical$platform == "HM450", na.rm = TRUE), 310L)
})

test_that("ANCHOR: every factor wired as a survival predictor is platform-clean", {
  # THE SELECTION CRITERION, MADE FALSIFIABLE. PLATFORM_CLEAN_MOFA_FACTORS is
  # the outcome-blind basis on which Phase 4 chooses its predictors, and it was
  # justified only by hard-coded numbers in comments citing run 30911448546 —
  # computed by a manually-triggered workflow that touched no target. A MOFA
  # seed / version / upstream-matrix change could make Factor1 or Factor4
  # platform-loaded and nothing would say so; re-adding "Factor2" to the vector
  # would have broken no check.
  #
  # This can only turn a green verdict RED. It never selects a factor and it
  # never looks at survival — doing either would be selection bias.
  fp <- read_pipeline_target("factor_platform")
  skip_if(is.null(fp), "factor_platform not in _targets store (run tar_make)")

  expect_s3_class(fp, "data.frame")
  expect_true(all(c("factor", "auc", "p_value", "q_value") %in% names(fp)))
  expect_true(all(is.finite(fp$auc)))
  # Scored on the real cohort, both arms present.
  expect_gte(unique(fp$n), COHORT_MIN)
  expect_lte(unique(fp$n), COHORT_MAX)

  # Every wired factor must exist and must be clean.
  expect_length(PLATFORM_CLEAN_MOFA_FACTORS, N_SURVIVAL_MOFA_FACTORS)
  expect_true(all(PLATFORM_CLEAN_MOFA_FACTORS %in% fp$factor))
  for (f in PLATFORM_CLEAN_MOFA_FACTORS) {
    q <- fp$q_value[fp$factor == f]
    expect_gt(q, SANITY_MAX_P)
    # A p of exactly 0 means the normal approximation degenerated, so both that
    # row's p AND its AUC are unusable — run 30911448546 recorded exactly that
    # for Factor6. A wired predictor must not be resting on such a row.
    expect_false(fp$degenerate[fp$factor == f])
  }

  # ... and the diagnostic must still have TEETH: if no factor at all came out
  # platform-loaded, the comparison has stopped discriminating and a clean
  # verdict for Factor1/Factor4 would be free. MEASURED (run 30911448546): 10 of
  # 15 factors are significant after BH, Factor2 at AUC 0.888.
  expect_gt(sum(fp$q_value <= SANITY_MAX_P), 0L)
})

test_that("ANCHOR: the wired factors are what the SELECTION RULE mechanically yields", {
  # THE OTHER HALF OF THE RULE. The anchor above executes clause (1) —
  # eligibility, q > SANITY_MAX_P and not degenerate. Clauses (2) RANK the
  # eligible factors by variance explained SUMMED across RNA + Methylation +
  # CNV and (3) TAKE THE TOP N_SURVIVAL_MOFA_FACTORS, and until now neither was
  # asserted anywhere: `grep -rn 'mofa_varexp' tests/` returned nothing. That
  # meant swapping Factor4 for Factor7, Factor13 or Factor15 passed every check
  # in the repo — INCLUDING swapping it because it gave a better held-out
  # C-index, which is exactly the selection bias the _targets.R comment
  # forbids. The rule claims to be DETERMINISTIC; this applies it.
  #
  # It can only turn a green verdict RED. It never chooses a factor itself, and
  # it never reads an outcome.
  fp <- read_pipeline_target("factor_platform")
  skip_if(is.null(fp), "factor_platform not in _targets store (run tar_make)")
  ve <- read_pipeline_target("mofa_varexp")
  skip_if(is.null(ve), "mofa_varexp not in _targets store (run tar_make)")

  ve <- as.matrix(ve)                       # rows = factors, cols = views
  expect_true(is.numeric(ve))
  expect_true(all(is.finite(ve)))
  expect_false(is.null(rownames(ve)))
  expect_true(all(fp$factor %in% rownames(ve)))

  # (1) eligibility -> (2) rank on summed variance explained -> (3) take top N.
  eligible <- fp$factor[fp$q_value > SANITY_MAX_P & !fp$degenerate]
  expect_gte(length(eligible), N_SURVIVAL_MOFA_FACTORS)
  total_ve <- rowSums(ve)[eligible]
  expect_false(anyNA(total_ve))
  ranked <- names(sort(total_ve, decreasing = TRUE))

  expect_identical(PLATFORM_CLEAN_MOFA_FACTORS,
                   head(ranked, N_SURVIVAL_MOFA_FACTORS))
})

test_that("ANCHOR: the two-platform overlap is MEASURED, not asserted", {
  # "Only 3 cases overlap the two platforms" is the SOLE stated justification
  # for refusing a batch correction, and it appeared in six places while being
  # computed in none. `methyl_platform_overlap` makes it executable.
  #
  # This anchor deliberately does NOT pin a magnitude — the number has never
  # been recorded, so pinning it would fabricate the very evidence it is meant
  # to supply. It asserts the two things that ARE derivable: the count is
  # well-formed, and it is consistent with methyl_platform's own rule, which
  # labels every dual-assayed case HM27. Overlap exceeding the HM27 arm would
  # mean that rule is not doing what its comment says.
  overlap <- read_pipeline_target("methyl_platform_overlap")
  skip_if(is.null(overlap),
          "methyl_platform_overlap not in _targets store (run tar_make)")
  cohort_ids <- read_pipeline_target("common_ids")
  platform <- read_pipeline_target("methyl_platform")

  expect_length(overlap, 1L)
  expect_true(is.numeric(overlap) && is.finite(overlap))
  expect_identical(overlap, as.integer(overlap))
  expect_gte(overlap, 0L)
  expect_lte(overlap, length(cohort_ids))
  # Dual-assayed cases are carried by their HM27 column, so they are a subset of
  # the HM27 arm.
  expect_lte(overlap, sum(platform == "HM27"))
})

test_that("ANCHOR: clinical covers colData, so platform is NA only outside the cohort", {
  # The frame deliberately spans every colData case (536) while methyl_platform
  # spans the 524-case main cohort; the difference, and ONLY the difference, is
  # allowed to be NA.
  clinical <- read_pipeline_target("clinical")
  skip_if(is.null(clinical), "clinical not in _targets store (run tar_make)")
  cohort_ids <- read_pipeline_target("common_ids")

  na_ids <- clinical$sample_id[is.na(clinical$platform)]
  expect_length(intersect(na_ids, cohort_ids), 0L)
  expect_identical(sum(!is.na(clinical$platform)), length(cohort_ids))
  # ... and the cohort really is a strict subset of the clinical frame.
  expect_length(setdiff(cohort_ids, clinical$sample_id), 0L)
})
