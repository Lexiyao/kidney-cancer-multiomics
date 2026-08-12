# Task 5.5 presentation helpers for dashboard/survival.qmd.
#
# The survival page is where this project is most likely to overclaim, so the
# two things that qualify the headline number are computed here rather than
# hand-typed into prose:
#
#   * OPTIMISM IS REPORTED, NEVER SUBTRACTED AWAY. `survival_metrics$optimism`
#     is apparent(train) - validated(held-out) for the Cox arm. The gap is the
#     finding; a page that showed only the held-out figure would be hiding how
#     much of the fit is in-sample flattery, and one that showed only the
#     apparent figure would be quoting the flattery itself.
#   * THE HELD-OUT FIGURE IS NOT FULLY OUT-OF-SAMPLE. `mofa_factors` was fitted
#     on all 524 cases, including the rows later called TEST, so the latent axes
#     were defined with held-out rows contributing. MOFA never sees an outcome,
#     so this is not a label leak -- it is unsupervised-transductive optimism,
#     and the C-index therefore bounds the SUPERVISED component only.
#
# Neither helper may invent a number: a missing metric comes back NA and the
# page prints a stated gap.

sm_fixture <- function() {
  list(
    cindex = list(cox = 0.62, penalised = 0.60, rsf = 0.58),
    optimism = list(cox = 0.07),
    calibration = data.frame(predicted = c(0.2, 0.8), observed = c(0.25, 0.75))
  )
}

test_that("fn_survival_cindex_table reports the held-out C-index per arm", {
  tab <- fn_survival_cindex_table(sm_fixture())

  expect_s3_class(tab, "data.frame")
  expect_equal(nrow(tab), 3L)
  expect_equal(tab$heldout, c(0.62, 0.60, 0.58))
})

test_that("fn_survival_cindex_table reports Cox optimism as its own column", {
  tab <- fn_survival_cindex_table(sm_fixture())

  expect_equal(tab$optimism[tab$arm == "Cox (primary)"], 0.07)
})

test_that("fn_survival_cindex_table derives the apparent figure, not the reverse", {
  # apparent = held-out + optimism. The held-out number is the one the pipeline
  # measured; the apparent one is shown only so the gap is visible.
  tab <- fn_survival_cindex_table(sm_fixture())

  expect_equal(tab$apparent[tab$arm == "Cox (primary)"], 0.69)
})

test_that("fn_survival_cindex_table leaves optimism NA where none was computed", {
  # Only the Cox arm has an optimism estimate. Copying the Cox gap onto the
  # penalised and RSF rows would be an invented number.
  tab <- fn_survival_cindex_table(sm_fixture())

  expect_true(all(is.na(tab$optimism[tab$arm != "Cox (primary)"])))
  expect_true(all(is.na(tab$apparent[tab$arm != "Cox (primary)"])))
})

test_that("fn_survival_cindex_table returns NA rather than guessing a missing arm", {
  sm <- sm_fixture()
  sm$cindex$rsf <- NULL

  tab <- fn_survival_cindex_table(sm)

  expect_true(is.na(tab$heldout[tab$arm == "Random survival forest"]))
})

test_that("fn_survival_cindex_table rejects a metrics object with no cindex", {
  expect_error(fn_survival_cindex_table(list(optimism = list(cox = 0.07))))
})

# --- The model-size panel -----------------------------------------------------
# Every figure here is already recorded by fn_fit_cox. Recomputing any of them
# on the page would create a second derivation of a number the pipeline owns.

cox_fixture <- function() {
  list(
    train = data.frame(time = 1:40, status = rep(0:1, 20)),
    test  = data.frame(time = 1:15, status = rep(0:1, length.out = 15)),
    predictors = c("Factor1", "Factor4", "age_years", "stage_num", "platform"),
    n_events = 124L, n_events_cohort = 171L,
    n_terms = 5L, max_predictors = 12L
  )
}

test_that("fn_cox_fit_summary reads the split and the EPV budget off the fit", {
  tab <- fn_cox_fit_summary(cox_fixture())

  expect_s3_class(tab, "data.frame")
  expect_named(tab, c("quantity", "value"))
  expect_true(any(grepl("^40$", tab$value)))   # training rows
  expect_true(any(grepl("^15$", tab$value)))   # held-out rows
  expect_true(any(grepl("124", tab$value)))    # training events (the EPV base)
  expect_true(any(grepl("5 of 12", tab$value)))
})

test_that("fn_cox_fit_summary names the predictors it was actually fitted on", {
  tab <- fn_cox_fit_summary(cox_fixture())

  expect_true(any(grepl("Factor1, Factor4", tab$value, fixed = TRUE)))
  # `subtype` is NOT a Cox term; the KM curves on the same page are descriptive
  # only, and a summary implying otherwise would misdescribe the model.
  expect_false(any(grepl("subtype", tab$value, fixed = TRUE)))
})

test_that("fn_cox_fit_summary refuses a fit missing the EPV record", {
  bad <- cox_fixture()
  bad$max_predictors <- NULL

  expect_error(fn_cox_fit_summary(bad))
})
