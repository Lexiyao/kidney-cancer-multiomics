test_that("phase-4 model constants are defined with sane values", {
  # Arrange / Act / Assert
  expect_identical(EPV_CAP, 10L)
  expect_true(HELDOUT_FRACTION > 0 && HELDOUT_FRACTION < 1)
  expect_identical(CV_FOLDS, 5L)
  expect_identical(MODEL_SEED, 20160128L)
  expect_identical(RSF_NTREE, 1000L)
  expect_identical(CALIBRATION_BINS, 5L)
  expect_true(SURVIVAL_HORIZON_DAYS == 5 * 365)
  expect_identical(BAP1_LABEL_COL, "BAP1")
})

test_that("fn_cindex returns 1 for perfectly ranked risk and 0 for reversed", {
  # Arrange
  time   <- c(1, 2, 3, 4)
  status <- c(1L, 1L, 1L, 0L)
  good   <- c(4, 3, 2, 1)   # highest risk dies earliest
  bad    <- c(1, 2, 3, 4)   # reversed
  # Act / Assert
  expect_equal(fn_cindex(time, status, good), 1)
  expect_equal(fn_cindex(time, status, bad), 0)
  expect_equal(fn_cindex(time, status, rep(0, 4)), 0.5)  # all tied
})

test_that("fn_cindex returns NA when no permissible pairs exist", {
  expect_true(is.na(fn_cindex(c(5, 6), c(0L, 0L), c(1, 2))))
})

test_that("fn_cindex matches survival::concordance on a random cohort", {
  # The from-scratch implementation is only credible if it reproduces the
  # reference. survival::concordance(Surv(t, s) ~ risk) treats a HIGHER
  # covariate as BETTER prognosis, so the reference for a mortality-style
  # risk score is taken on -risk.
  # Arrange
  skip_if_not_installed("survival")
  set.seed(MODEL_SEED)
  n <- 120
  risk <- stats::rnorm(n)
  time <- stats::rexp(n, rate = exp(risk) / 1000)
  status <- stats::rbinom(n, 1L, 0.6)
  # Act
  mine <- fn_cindex(time, as.integer(status), risk)
  ref <- survival::concordance(
    survival::Surv(time, status) ~ risk, reverse = TRUE
  )$concordance
  # Assert
  expect_equal(mine, ref, tolerance = 1e-8)
  expect_gt(mine, 0.5)  # risk was generated to be genuinely prognostic
})

test_that("fn_cindex scores tied event times and tied risks by the Harrell rule", {
  # Arrange: deliberate ties in BOTH time and risk.
  time   <- c(10, 10, 20, 20, 30, 30)
  status <- c(1L, 0L, 1L, 1L, 0L, 1L)
  risk   <- c(2, 2, 1, 3, 1, 1)
  # Act
  mine <- fn_cindex(time, status, risk)
  # Assert -- hand count. A pair is permissible when subject i has an event and
  # subject j OUTLIVES i: a strictly later time, OR the SAME recorded time with
  # j censored (j was known to be event-free when i died, so j outlived i --
  # the survival::concordance convention).
  #   i=1 (t=10, r=2) vs t=20,20,30,30 -> risks 1,3,1,1 : 3 concordant, 1 disc
  #                       vs t=10 CENSORED, r=2         : 1 tied risk
  #   i=3 (t=20, r=1) vs t=30,30       -> risks 1,1     : 0 concordant, 2 tied
  #   i=4 (t=20, r=3) vs t=30,30       -> risks 1,1     : 2 concordant
  #   i=6 (t=30, r=1) vs t=30 CENSORED, r=1             : 1 tied risk
  # permissible 10, concordant 5, tied-risk 4 -> (5 + 2) / 10 = 0.70
  expect_equal(mine, 0.70)
})

test_that("fn_cindex agrees with survival::concordance on same-time censoring", {
  # THE REGRESSION THIS EXISTS TO CATCH, and it is not hypothetical: an earlier
  # implementation required a STRICTLY later time and so discarded every pair
  # where subject i died at t and subject j was censored at t. That pinned a
  # documented "divergence" from the reference as if it were intentional.
  # MEASURED (R 4.6.0 / survival 3.8.6): on continuous times the two rules
  # agree to 1.1e-16, but on INTEGER-day times (n=250, 33% events, 25 cohorts)
  # the old rule drifted from the reference by up to 4.7e-04 and on a 50-day
  # grid by up to 1.1e-02. TCGA days_to_death / days_to_last_followup ARE
  # integer days, so the real 522-case fit sits squarely in that regime and the
  # reported C-index would not be reproducible with the standard reference.
  # Arrange
  skip_if_not_installed("survival")
  time   <- c(10, 10, 20, 20, 30, 30)
  status <- c(1L, 0L, 1L, 1L, 0L, 1L)
  risk   <- c(2, 2, 1, 3, 1, 1)
  # Act
  cc <- survival::concordance(
    survival::Surv(time, status) ~ risk, reverse = TRUE
  )
  # Assert
  expect_equal(unname(cc$count[["concordant"]]), 5)
  expect_equal(unname(cc$count[["discordant"]]), 1)
  expect_equal(cc$concordance, 0.70, tolerance = 1e-8)
  expect_equal(fn_cindex(time, status, risk), cc$concordance, tolerance = 1e-8)
})

test_that("fn_cindex matches survival::concordance on INTEGER-day cohorts", {
  # The random-cohort test above uses rexp(), which produces no ties at all, so
  # it cannot exercise the tie rule. TCGA follow-up is recorded in whole days;
  # these cohorts are rounded to reproduce that, and one is put on a coarse
  # 50-day grid so tied times are unmissable.
  skip_if_not_installed("survival")
  for (seed in 1:6) {
    set.seed(seed)
    n <- 200
    risk <- stats::rnorm(n)
    raw <- stats::rexp(n, rate = exp(risk) / 1000)
    for (times in list(round(raw), round(raw / 50) * 50)) {
      status <- as.integer(stats::rbinom(n, 1L, 0.33))
      expect_true(anyDuplicated(times) > 0)  # the ties are really there
      ref <- survival::concordance(
        survival::Surv(times, status) ~ risk, reverse = TRUE
      )$concordance
      expect_equal(fn_cindex(times, status, risk), ref, tolerance = 1e-10)
    }
  }
})

test_that("fn_calibration returns one row per non-empty bin with valid probabilities", {
  # Arrange
  set.seed(1)
  n <- 200
  pred_surv <- runif(n, 0.1, 0.9)
  time <- rexp(n, rate = 1 / (pred_surv * 2000))  # better predicted -> longer time
  status <- rbinom(n, 1, 0.5)
  # Act
  cal <- fn_calibration(time, as.integer(status), pred_surv,
                        horizon = SURVIVAL_HORIZON_DAYS, n_bins = 5L)
  # Assert
  expect_s3_class(cal, "data.frame")
  expect_equal(nrow(cal), 5L)
  expect_named(cal, c("bin", "n", "predicted", "observed"))
  expect_true(all(cal$predicted >= 0 & cal$predicted <= 1))
  expect_true(all(cal$observed >= 0 & cal$observed <= 1, na.rm = TRUE))
  expect_equal(sum(cal$n), n)
})

test_that("fn_calibration reproduces a hand-computed Kaplan-Meier at the horizon", {
  # A shape test cannot catch a wrong `observed`. This one pins the number.
  # Arrange: one bin (n_bins = 1L) so the whole cohort is a single KM curve.
  # Deaths at 100 and 300, censoring at 200, horizon 400.
  #   KM: S(100) = 1 - 1/4 = 0.75
  #       t=200 censored -> no drop, 2 at risk after
  #       S(300) = 0.75 * (1 - 1/2) = 0.375
  #   S(400) = 0.375 (carried forward past the last event by extend = TRUE)
  # The four predicted probabilities differ (mean exactly 0.4) but all land in
  # the single bin. They must NOT be constant: `fn_calibration` cuts on
  # quantiles, and a perfectly constant prediction vector yields non-unique
  # breaks and errors out. That is a real edge case of the implementation --
  # a degenerate model predicting one probability for everyone stops loudly
  # rather than returning a one-row table.
  time      <- c(100, 200, 300, 500)
  status    <- c(1L, 0L, 1L, 0L)
  pred_surv <- c(0.30, 0.35, 0.45, 0.50)
  # Act
  cal <- fn_calibration(time, status, pred_surv, horizon = 400, n_bins = 1L)
  # Assert
  expect_equal(nrow(cal), 1L)
  expect_equal(cal$n, 4L)
  expect_equal(cal$predicted, 0.4)
  expect_equal(cal$observed, 0.375, tolerance = 1e-8)
})

test_that("fn_calibration bins are monotone in predicted survival and match survfit", {
  # Arrange: predicted survival ordered by construction, so bin 1 must carry
  # the lowest mean predicted probability and bin 5 the highest.
  skip_if_not_installed("survival")
  set.seed(MODEL_SEED)
  n <- 250
  pred_surv <- stats::runif(n, 0.05, 0.95)
  time <- stats::rexp(n, rate = 1 / (pred_surv * 2000))
  status <- stats::rbinom(n, 1L, 0.7)
  # Act
  cal <- fn_calibration(time, as.integer(status), pred_surv,
                        horizon = SURVIVAL_HORIZON_DAYS, n_bins = CALIBRATION_BINS)
  # Assert: bin ordering, and `observed` recomputed independently for one bin
  expect_identical(cal$bin, seq_len(CALIBRATION_BINS))
  expect_false(is.unsorted(cal$predicted))
  brks <- stats::quantile(pred_surv, probs = seq(0, 1, length.out = CALIBRATION_BINS + 1))
  idx <- which(cut(pred_surv, breaks = brks, include.lowest = TRUE, labels = FALSE) == 3L)
  km <- survival::survfit(survival::Surv(time[idx], status[idx]) ~ 1)
  expect_equal(
    cal$observed[cal$bin == 3L],
    summary(km, times = SURVIVAL_HORIZON_DAYS, extend = TRUE)$surv,
    tolerance = 1e-10
  )
  expect_equal(sum(cal$n), n)
})

test_that("fn_calibration rejects mismatched input lengths", {
  expect_error(
    fn_calibration(c(1, 2, 3), c(1L, 0L, 1L), c(0.5, 0.5),
                   horizon = SURVIVAL_HORIZON_DAYS, n_bins = 2L)
  )
})

# Synthetic survival frame shared by the Module-4 fitting tests.
#
# DEVIATION FROM THE PLAN, RECORDED: the plan's Task 4.4 body for this helper
# emits only time/status/Factor1-3/age_years/stage_num, but the plan's own
# penalised-Cox test calls `make_surv_df(n = 400, signal = TRUE)` and asks for
# `Factor4` and a two-level `platform` factor. The helper is therefore written
# once, here, with the `signal` switch and the `platform` factor the later test
# requires. `platform` carries BOTH METHYL_PLATFORMS levels as a DECLARED level
# set so train and test designs code identically even when a split happens to
# contain one platform only.
make_surv_df <- function(n = 160, seed = 42, signal = FALSE) {
  set.seed(seed)
  factor1 <- rnorm(n)
  factor2 <- rnorm(n)
  factor3 <- rnorm(n)
  factor4 <- rnorm(n)
  age_years <- round(runif(n, 40, 80))
  stage_num <- sample(1:4, n, replace = TRUE)
  platform <- factor(sample(METHYL_PLATFORMS, n, replace = TRUE),
                     levels = METHYL_PLATFORMS)
  # `signal = TRUE` puts genuine prognostic structure on the log-hazard so that
  # LASSO has something to keep; on pure noise a shrunk-to-zero fit is a
  # legitimate outcome and cannot distinguish a correct design from a broken one.
  lp <- if (signal) {
    0.9 * factor1 + 0.7 * factor4 + 0.03 * (age_years - 60) +
      0.5 * (platform == "HM450")
  } else {
    rep(0, n)
  }
  data.frame(
    time   = rexp(n, rate = exp(lp) / 1500),
    status = rbinom(n, 1, 0.6),
    Factor1 = factor1, Factor2 = factor2, Factor3 = factor3, Factor4 = factor4,
    age_years = age_years,
    stage_num = stage_num,
    platform = platform
  )
}

test_that("fn_max_predictors applies the EPV cap by flooring events/EPV", {
  expect_identical(fn_max_predictors(100L), 10L)
  expect_identical(fn_max_predictors(95L), 9L)
  expect_identical(fn_max_predictors(9L), 0L)
})

test_that("fn_split_train_test is deterministic and partitions all rows", {
  df <- make_surv_df()
  a <- fn_split_train_test(df, 0.3, seed = 1L)
  b <- fn_split_train_test(df, 0.3, seed = 1L)
  expect_equal(a$test, b$test)                      # deterministic
  expect_equal(nrow(a$train) + nrow(a$test), nrow(df))
  expect_equal(nrow(a$test), floor(nrow(df) * 0.3))
})

test_that("fn_fit_cox fits on a held-out split and returns test risk scores", {
  df <- make_surv_df()
  fit <- fn_fit_cox(df, c("Factor1", "Factor2", "age_years"))
  expect_s3_class(fit$model, "coxph")
  expect_length(fit$risk_test, nrow(fit$test))
  expect_true(all(is.finite(fit$risk_test)))
})

test_that("fn_fit_cox throws when predictors exceed the EPV cap", {
  # MEASURED with this helper: n = 40 gives 25 cohort events (cap 2) and 16
  # events in the TRAINING partition the model is actually estimated on
  # (cap 1), so the guard throws under either counting rule.
  df <- make_surv_df(n = 40)
  expect_error(
    fn_fit_cox(df, c("Factor1", "Factor2", "Factor3", "age_years", "stage_num")),
    "EPV cap violated"
  )
})

test_that("the EPV cap counts the events in TRAIN, not the whole cohort", {
  # THE DEFECT THIS PINS: the cap used to be evaluated on the whole input frame
  # while coxph was fitted on the 70% training partition, so it licensed a
  # budget ~40% larger than the data supporting the fit. MEASURED with this
  # helper at n = 100: 65 cohort events (cap 6, which ADMITS these five
  # predictors) but 46 training events (cap 4, which does not).
  df <- make_surv_df(n = 100)
  preds <- c("Factor1", "Factor2", "Factor3", "age_years", "stage_num")
  train <- fn_split_train_test(df, HELDOUT_FRACTION)$train
  expect_gte(fn_max_predictors(sum(df$status == 1L)), length(preds))   # cohort admits
  expect_lt(fn_max_predictors(sum(train$status == 1L)), length(preds))  # train does not
  expect_error(fn_fit_cox(df, preds), "EPV cap violated")

  # ... and the counts the fit REPORTS describe the rows it was fitted on.
  ok <- fn_fit_cox(df, c("Factor1", "age_years"))
  expect_identical(ok$n_events, sum(train$status == 1L))
  expect_identical(ok$max_predictors, fn_max_predictors(sum(train$status == 1L)))
  expect_identical(ok$n_events_cohort, sum(df$status == 1L))
})

test_that("the EPV cap spends DESIGN DF, not model terms, on a multi-level factor", {
  # A k-level factor costs k-1 coefficients but is one term. Counting terms let
  # a 3-level covariate through a budget that could not afford it. The
  # penalised arm already codes the design this way (model.matrix(~ .)[, -1]),
  # so both arms are now counted identically.
  # MEASURED: n = 60 -> 25 training events -> cap 2. `Factor1 + site` is 2
  # TERMS (which the old rule allowed) but 1 + 2 = 3 DF.
  df <- make_surv_df(n = 60)
  set.seed(7)
  df$site <- factor(sample(c("A", "B", "C"), nrow(df), replace = TRUE),
                    levels = c("A", "B", "C"))
  train <- fn_split_train_test(df, HELDOUT_FRACTION)$train
  expect_identical(fn_max_predictors(sum(train$status == 1L)), 2L)
  expect_error(fn_fit_cox(df, c("Factor1", "site")), "EPV cap violated")
  # The message reports the DF spent, not the term count.
  expect_error(fn_fit_cox(df, c("Factor1", "site")), "3 predictors requested")
})

test_that("fn_fit_penalised_cox returns cv.glmnet fit and finite held-out risk", {
  df <- make_surv_df(n = 180)
  fit <- fn_fit_penalised_cox(df, c("Factor1", "Factor2", "Factor3", "age_years"))
  expect_s3_class(fit$model, "cv.glmnet")
  expect_true(is.finite(fit$lambda_min))
  expect_length(fit$risk_test, nrow(fit$test))
  expect_true(all(is.finite(fit$risk_test)))
})

test_that("fn_fit_penalised_cox dummy-codes the platform factor instead of NA-ing it", {
  # THE REGRESSION THIS EXISTS TO CATCH: `as.matrix` on a data.frame carrying a
  # factor yields a CHARACTER matrix; glmnet re-coerces to NA, fits an
  # all-zero model and returns a CONSTANT risk vector. Finiteness alone does
  # not catch it -- a vector of zeros is perfectly finite -- and fn_cindex then
  # reports ~0.5 without complaint.
  #
  # `make_surv_df` must therefore carry a two-level `platform` factor AND
  # genuine signal, otherwise LASSO can legitimately shrink everything to zero
  # at lambda.min and the sd assertion fails for an innocent reason.
  #
  # PLACEMENT NOTE: the plan prints this test inside Task 4.4's block, but it
  # exercises `fn_fit_penalised_cox`, which Task 4.4 does not create. Running it
  # there would make Task 4.4's own GREEN step impossible, so it lives with the
  # function it tests. Text otherwise unchanged.
  df <- make_surv_df(n = 400, signal = TRUE)
  fit <- fn_fit_penalised_cox(df, c("Factor1", "Factor4", "age_years", "platform"))
  expect_true("platformHM450" %in% fit$design_cols)
  expect_false("platform" %in% fit$design_cols)
  expect_gt(stats::sd(fit$risk_test), 0)
})

test_that("fn_fit_rsf fits a survival forest and returns finite mortality on test", {
  df <- make_surv_df(n = 160)
  fit <- fn_fit_rsf(df, c("Factor1", "Factor2", "age_years"), n_tree = 100L)
  expect_s3_class(fit$model, "rfsrc")
  expect_length(fit$risk_test, nrow(fit$test))
  expect_true(all(is.finite(fit$risk_test)))
})

test_that("fn_fit_rsf holds out the SAME rows as the Cox arm", {
  # The three arms are compared to each other on the held-out set, so they must
  # be scored on the SAME held-out rows. All three route the split through
  # fn_split_train_test with the same default seed; this pins that, because a
  # per-model split would make the C-index comparison meaningless without
  # anything failing.
  df <- make_surv_df(n = 200)
  preds <- c("Factor1", "Factor2", "age_years")
  rsf <- fn_fit_rsf(df, preds, n_tree = 50L)
  cox <- fn_fit_cox(df, preds)
  expect_identical(rownames(rsf$test), rownames(cox$test))
})

# --- Discrimination floor + null detector ----------------------------------
# WHY THESE TWO BLOCKS EXIST, MEASURED: replacing risk_train/risk_test in all
# three fitters with stats::rnorm() left the suite at 90 PASS / 0 FAIL —
# byte-identical to baseline. Replacing them with a CONSTANT failed only 2
# assertions, both `expect_gt(sd(...), 0)` on the penalised and RSF arms; the
# Cox arm was caught by NOTHING, and Cox alone feeds survival_metrics$cindex
# $cox, $optimism$cox and the whole 5-year calibration table. Every other Cox
# assertion is shape-only, and fn_cindex on noise returns ~0.5 without
# complaint — the silent-green failure mode the comment block in
# R/functions_survival.R says this repo exists to prevent. The Python side has
# had a null detector since Task 4.8; this is its R counterpart.
#
# NOTE ON THE THRESHOLDS: they are floors on a SYNTHETIC fixture whose
# log-hazard is known by construction, not a target for the real cohort. They
# can only fail a fit that has stopped ranking; nothing here is tuned, and no
# predictor is chosen on any outcome.

test_that("every arm discriminates on a fixture with known prognostic signal", {
  # make_surv_df(signal = TRUE) puts lp = 0.9*Factor1 + 0.7*Factor4 +
  # 0.03*(age-60) + 0.5*(platform == "HM450") on the hazard, so held-out
  # concordance is genuinely above chance. MEASURED at n = 400: cox 0.771
  # (train 0.746), penalised 0.770, rsf 0.768.
  df <- make_surv_df(n = 400, signal = TRUE)
  preds <- c("Factor1", "Factor4", "age_years", "platform")

  cox <- fn_fit_cox(df, preds)
  expect_gt(stats::sd(cox$risk_test), 0)
  expect_gt(fn_cindex(cox$test$time, cox$test$status, cox$risk_test), 0.6)
  # risk_train is pinned too: it is the apparent half of the reported optimism.
  expect_gt(fn_cindex(cox$train$time, cox$train$status, cox$risk_train), 0.6)

  pen <- fn_fit_penalised_cox(df, preds)
  expect_gt(stats::sd(pen$risk_test), 0)
  expect_gt(fn_cindex(pen$test$time, pen$test$status, pen$risk_test), 0.6)

  rsf <- fn_fit_rsf(df, preds, n_tree = 200L)
  expect_gt(stats::sd(rsf$risk_test), 0)
  expect_gt(fn_cindex(rsf$test$time, rsf$test$status, rsf$risk_test), 0.6)
})

test_that("no arm discriminates on a fixture with NO signal (leak detector)", {
  # The mirror image, and the reason the block above cannot be satisfied by
  # leaking the outcome into the features: on pure noise the held-out C-index
  # must sit at chance. MEASURED at n = 400, signal = FALSE: cox 0.530,
  # penalised 0.534, rsf 0.577 — all inside the band below. A fit that scored
  # 0.77 here would be reading the held-out rows.
  df <- make_surv_df(n = 400, signal = FALSE)
  preds <- c("Factor1", "Factor4", "age_years", "platform")

  cox <- fn_fit_cox(df, preds)
  c_cox <- fn_cindex(cox$test$time, cox$test$status, cox$risk_test)
  expect_gt(c_cox, 0.35)
  expect_lt(c_cox, 0.65)

  pen <- fn_fit_penalised_cox(df, preds)
  c_pen <- fn_cindex(pen$test$time, pen$test$status, pen$risk_test)
  expect_gt(c_pen, 0.35)
  expect_lt(c_pen, 0.65)

  rsf <- fn_fit_rsf(df, preds, n_tree = 200L)
  c_rsf <- fn_cindex(rsf$test$time, rsf$test$status, rsf$risk_test)
  expect_gt(c_rsf, 0.35)
  expect_lt(c_rsf, 0.65)
})

test_that("fn_fit_rsf carries the platform FACTOR through the formula, not as NA", {
  # The glmnet arm needs an explicit model.matrix (see fn_fit_penalised_cox);
  # the forest goes through a formula and must handle the two-level factor
  # itself. A degenerate forest would return a CONSTANT mortality vector, which
  # is finite and would slip past the finiteness assertion above while
  # fn_cindex quietly reported ~0.5.
  df <- make_surv_df(n = 300, signal = TRUE)
  fit <- fn_fit_rsf(df, c("Factor1", "Factor4", "age_years", "platform"),
                    n_tree = 100L)
  expect_true(all(is.finite(fit$risk_test)))
  expect_gt(stats::sd(fit$risk_test), 0)
})

# --- Anchors on the real Module-4 targets -----------------------------------
# WHAT THESE FIX: before them, `grep -rn 'survival_df|cox_fit|survival_metrics|
# bap1_auroc' tests/ .github/` returned ZERO hits outside prose. Every test
# above scores PURE FUNCTIONS on synthetic frames; the DAG targets themselves
# were built by an unrestricted `tar_make()` under `continue-on-error: true`
# and read by nothing, so a zero-event join, a degenerate design, an EPV throw
# or a missing scikit-learn could leave CI green with no artifact recording it.
# The repo has twice been bitten by a skip that never lifts; this was worse —
# there was not even a skip.
#
# They use the same read_pipeline_target() contract as the Module-3 anchors:
# NULL (no pipeline metadata at all) skips, a populated store missing or
# erroring the target FAILS. The `ANCHOR:` prefix is required — LEVEL 1 of
# verify-module2.yml hard-fails on any NON-anchor skip, and LEVEL 3 is where
# these must execute. They assert bands and well-formedness only; no threshold
# below is a target, and nothing here selects a predictor or tunes a metric.

test_that("ANCHOR: survival_df carries the measured OS event count", {
  df <- read_pipeline_target("survival_df")
  skip_if(is.null(df), "survival_df not in _targets store (run tar_make)")
  preds <- read_pipeline_target("survival_predictors")

  expect_s3_class(df, "data.frame")
  expect_true(all(c("time", "status") %in% names(df)))
  expect_true(all(preds %in% names(df)))
  # MEASURED (run 30708943504): 173 OS events among 522 usable cases on the
  # 524-case main cohort. survival_df additionally drops rows with a missing
  # predictor and rows with time <= 0, so it can only lose events, never gain
  # them; the lower bound leaves room for that filtering without letting a
  # collapsed join through.
  n_events <- sum(df$status == 1L)
  expect_gte(n_events, 140L)
  expect_lte(n_events, 173L)
  expect_gte(nrow(df), COHORT_MIN - 25L)
  expect_lte(nrow(df), COHORT_MAX)
  # The adjustment covariate must be complete and genuinely two-armed:
  # coxph deletes incomplete rows silently, and a constant covariate adjusts
  # for nothing.
  expect_false(anyNA(df$platform))
  expect_identical(nlevels(droplevels(df$platform)), 2L)
  expect_true(all(df$time > 0))
})

test_that("ANCHOR: cox_fit ranks the held-out rows instead of scoring a constant", {
  fit <- read_pipeline_target("cox_fit")
  skip_if(is.null(fit), "cox_fit not in _targets store (run tar_make)")

  expect_s3_class(fit$model, "coxph")
  expect_length(fit$risk_test, nrow(fit$test))
  expect_length(fit$risk_train, nrow(fit$train))
  expect_true(all(is.finite(fit$risk_test)))
  expect_true(all(is.finite(fit$risk_train)))
  # sd > 0, not just finiteness: an all-zero-coefficient fit returns a CONSTANT
  # risk vector, which is perfectly finite and makes fn_cindex report ~0.5
  # without complaint (the failure mode documented in R/functions_survival.R).
  expect_gt(stats::sd(fit$risk_test), 0)
  expect_gt(stats::sd(fit$risk_train), 0)
  # The EPV budget is measured on the rows the model was estimated on.
  expect_identical(fit$n_events, sum(fit$train$status == 1L))
  expect_lte(fit$n_terms, fit$max_predictors)
  expect_identical(fit$max_predictors, fn_max_predictors(fit$n_events))
})

test_that("ANCHOR: survival_metrics reports a held-out C-index and the optimism", {
  m <- read_pipeline_target("survival_metrics")
  skip_if(is.null(m), "survival_metrics not in _targets store (run tar_make)")

  expect_named(m, c("cindex", "optimism", "calibration"))
  for (arm in c("cox", "penalised", "rsf")) {
    expect_true(is.finite(m$cindex[[arm]]))
    expect_gte(m$cindex[[arm]], 0)
    expect_lte(m$cindex[[arm]], 1)
  }
  # Optimism is REPORTED, never subtracted away, so it only has to exist and be
  # finite. Its SIGN is a finding, not a requirement, and is not asserted.
  expect_true(is.finite(m$optimism$cox))
  expect_s3_class(m$calibration, "data.frame")
  expect_identical(nrow(m$calibration), as.integer(CALIBRATION_BINS))
  expect_named(m$calibration, c("bin", "n", "predicted", "observed"))
  expect_true(all(m$calibration$predicted >= 0 & m$calibration$predicted <= 1))
  expect_true(all(m$calibration$observed >= 0 & m$calibration$observed <= 1,
                  na.rm = TRUE))
})

test_that("ANCHOR: bap1_auroc is finite and scored on the mutation subset", {
  a <- read_pipeline_target("bap1_auroc")
  skip_if(is.null(a), "bap1_auroc not in _targets store (run tar_make)")

  expect_true(all(c("cv_auroc", "heldout_auroc", "n_samples",
                    "n_bap1_mutant") %in% names(a)))
  for (nm in c("cv_auroc", "heldout_auroc")) {
    expect_true(is.finite(a[[nm]]))
    expect_gte(a[[nm]], 0)
    expect_lte(a[[nm]], 1)
  }
  # The label is EXTERNAL and the subset is the n=417 mutation cohort; an
  # AUROC computed on a handful of samples, or with no mutant at all, is not a
  # result. Same guard fn_annotate_mutation applies to the same join.
  expect_gte(a$n_samples, MIN_MUT_ANNOT_SAMPLES)
  expect_gt(a$n_bap1_mutant, 0L)
  expect_lt(a$n_bap1_mutant, a$n_samples)
})
