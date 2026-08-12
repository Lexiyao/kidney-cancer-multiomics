# Module 4 survival: low-dimensional Cox / penalised Cox / RSF on
# factors/subtypes + a few clinical variables, EPV-capped, held-out split.

fn_max_predictors <- function(n_events, epv_cap = EPV_CAP) {
  # Events-per-variable cap: at EPV~10 the model is limited to
  # floor(events / EPV) predictors (spec §2, load-bearing). `n_events` is the
  # OS event count MEASURED on the rows the model is ESTIMATED on, which is the
  # TRAINING partition, not the whole cohort (see fn_enforce_epv_cap). The
  # cohort census — 173 events on the 524-case main cohort, 148 within 5 years
  # — is what licenses the design; the enforced budget is the training-partition
  # figure, roughly 0.7 of it. It is still computed rather than hard-coded, so
  # the budget stays honest if the cohort or the filtering changes; the design
  # spends only 5.
  as.integer(floor(n_events / epv_cap))
}

fn_split_train_test <- function(surv_df, split, seed = MODEL_SEED) {
  # Deterministic held-out split; `split` is the test-set fraction.
  set.seed(seed)
  n <- nrow(surv_df)
  test_idx <- sample.int(n, size = floor(n * split))
  list(train = surv_df[-test_idx, , drop = FALSE],
       test  = surv_df[test_idx, , drop = FALSE])
}

# TWO THINGS THIS COUNTS, AND WHY EACH IS THE FITTED QUANTITY.
#
# ROWS: `surv_df` here must be the TRAINING partition. The cap used to be
# evaluated on the whole input frame while coxph was fitted on the 70% train
# split, so it granted a budget ~40% larger than the data supporting the fit.
# MEASURED on a 522-row frame at the repo's own MODEL_SEED/HELDOUT_FRACTION:
# 168 cohort events -> cap 16, but 120 training events -> cap 12.
#
# COLUMNS: `length(predictors)` counts model TERMS; a k-level factor spends
# k-1 coefficients while counting as one. The design width is taken from the
# SAME coding fn_fit_penalised_cox uses (model.matrix(~ .)[, -1]) so both arms
# are charged identically. Neither error binds at the 5 wired predictors — 2
# factors + age + stage + a 2-level platform is 5 terms and 5 df — but both ran
# in the permissive direction, so the number printed was not the number the fit
# could afford.
fn_enforce_epv_cap <- function(surv_df, predictors) {
  n_events <- sum(surv_df$status == 1L)
  max_pred <- fn_max_predictors(n_events)
  design <- stats::model.matrix(~ ., data = surv_df[, predictors, drop = FALSE])
  n_terms <- ncol(design[, -1, drop = FALSE])
  if (n_terms > max_pred) {
    # Format string hoisted to a local only so the call satisfies the repo's
    # indentation linter (CI fails on ANY lint under R/); wording unchanged.
    fmt <- "EPV cap violated: %d predictors requested, only %d allowed at %d events (EPV=%d)"
    stop(sprintf(fmt, n_terms, max_pred, n_events, EPV_CAP), call. = FALSE)
  }
  list(n_events = n_events, max_predictors = max_pred, n_terms = n_terms)
}

fn_fit_cox <- function(surv_df, predictors, split = HELDOUT_FRACTION,
                       seed = MODEL_SEED) {
  stopifnot(all(c("time", "status") %in% names(surv_df)))
  # SPLIT FIRST, then cap: the guard must measure the rows coxph is estimated
  # on. `n_events_cohort` keeps the cohort-level census quotable for the record
  # without being the number enforced.
  parts <- fn_split_train_test(surv_df, split, seed)
  cap <- fn_enforce_epv_cap(parts$train, predictors)
  rhs <- paste(predictors, collapse = " + ")
  form <- stats::as.formula(paste("survival::Surv(time, status) ~", rhs))
  fit <- survival::coxph(form, data = parts$train, x = TRUE)
  list(model = fit,
       predictors = predictors,
       train = parts$train,
       test = parts$test,
       # stats::predict, NOT bare predict: tar_make() runs every target in one
       # process, so packages attached by earlier targets accumulate and MOFA2
       # (attached for the mofa_* targets) masks the generic with a signature
       # that has no `type`. predict.coxph is not exported, so qualifying the
       # METHOD is impossible -- qualify the GENERIC instead.
       risk_train = as.numeric(stats::predict(fit, type = "lp")),
       risk_test  = as.numeric(stats::predict(fit, newdata = parts$test, type = "lp")),
       n_events = cap$n_events,
       n_events_cohort = sum(surv_df$status == 1L),
       n_terms = cap$n_terms,
       max_predictors = cap$max_predictors)
}

# `as.matrix` IS NOT SAFE HERE ANY MORE. `survival_predictors` now contains the
# `platform` FACTOR (Task 4.7), and `as.matrix` on a data.frame holding a
# factor coerces the WHOLE design to a CHARACTER matrix; glmnet then
# re-coerces with as.numeric and the platform column becomes all-NA.
#
# MEASURED (R 4.6.0, glmnet 5.0, n = 200, 5 predictors incl. platform):
# `cv.glmnet` does NOT error. It emits only "NAs introduced by coercion" and
# returns a model with EVERY coefficient zero at lambda.min, and `risk_test`
# identically 0. Note what that does to this task's own acceptance test:
# `expect_true(all(is.finite(fit$risk_test)))` PASSES on that fit, and
# `fn_cindex` on a constant risk vector returns ~0.5 without complaint. That
# is the silent-green failure mode the rest of this repo is built to prevent.
# `fn_fit_cox` and `fn_fit_rsf` are unaffected — both go through a formula,
# which dummy-codes the factor correctly.
#
# CODING: `model.matrix(~ ., ...)[, -1]` (treatment contrasts, intercept
# dropped), NOT `model.matrix(~ . - 1, ...)`. The latter emits BOTH
# platformHM27 and platformHM450, costing 2 df where Task 4.7's comment and
# `fn_fit_cox`'s formula both spend 1, so the penalised arm would be modelling
# a differently-coded covariate from the Cox arm it is compared against.
# MEASURED: 6 columns vs 5, same held-out risk spread.
#
# TRAIN AND TEST MUST BE CODED IDENTICALLY. `platform` is a factor with both
# METHYL_PLATFORMS levels DECLARED, so model.matrix emits the same columns
# even for a split that happens to contain one level only (VERIFIED: an
# all-HM27 subset still yields 5 columns). Keep it a declared-level factor;
# a character column would silently produce a narrower test design.
#
# NO LEAKAGE ACROSS THE SPLIT: the design is built per-part, cv.glmnet sees the
# TRAIN matrix only, and lambda.min is chosen by CV inside TRAIN. The test part
# is only ever scored.
fn_fit_penalised_cox <- function(surv_df, predictors, split = HELDOUT_FRACTION,
                                 seed = MODEL_SEED, n_folds = CV_FOLDS) {
  stopifnot(all(c("time", "status") %in% names(surv_df)))
  parts <- fn_split_train_test(surv_df, split, seed)
  fn_design <- function(d) {
    stats::model.matrix(~ ., data = d[, predictors, drop = FALSE])[, -1, drop = FALSE]
  }
  x_train <- fn_design(parts$train)
  x_test  <- fn_design(parts$test)
  # HARD GUARD: an all-NA or character design must stop, never fit.
  stopifnot(is.numeric(x_train), is.numeric(x_test),
            !anyNA(x_train), !anyNA(x_test),
            identical(colnames(x_train), colnames(x_test)))
  y_train <- survival::Surv(parts$train$time, parts$train$status)
  set.seed(seed)
  # cox.ties pinned explicitly: glmnet 5.0 defaults to "breslow" but warns that
  # 5.1 will switch to "efron". Leaving it implicit means a routine dependency
  # bump silently changes every C-index this pipeline reports. "breslow" keeps
  # the numbers reproducible against the recorded runs; changing it is a
  # deliberate act that must be re-verified, not a side effect of an upgrade.
  cvfit <- glmnet::cv.glmnet(x_train, y_train, family = "cox",
                             alpha = 1, nfolds = n_folds,
                             cox.ties = "breslow")
  # stats::predict for the same masking reason as fn_fit_cox above.
  risk_test <- as.numeric(stats::predict(cvfit, newx = x_test, s = "lambda.min",
                                         type = "link"))
  list(model = cvfit,
       lambda_min = cvfit$lambda.min,
       predictors = predictors,
       design_cols = colnames(x_train),
       test = parts$test,
       risk_test = risk_test)
}

# The forest goes through a FORMULA, so unlike the glmnet arm it dummy-codes the
# `platform` factor itself and needs no model.matrix. Same split helper and the
# same default seed as fn_fit_cox / fn_fit_penalised_cox, so all three arms are
# scored on IDENTICAL held-out rows and their C-indices are comparable.
#
# NO LEAKAGE: the forest is grown on TRAIN only and the test part is scored by
# predict.rfsrc; no tuning touches the held-out rows. `seed = -abs(seed)`
# because rfsrc requires a NEGATIVE integer to seed its RNG reproducibly.
fn_fit_rsf <- function(surv_df, predictors, split = HELDOUT_FRACTION,
                       seed = MODEL_SEED, n_tree = RSF_NTREE) {
  stopifnot(all(c("time", "status") %in% names(surv_df)))
  parts <- fn_split_train_test(surv_df, split, seed)
  rhs <- paste(predictors, collapse = " + ")
  form <- stats::as.formula(paste("Surv(time, status) ~", rhs))
  fit <- randomForestSRC::rfsrc(form, data = parts$train,
                                ntree = n_tree, importance = "none",
                                seed = -abs(seed))
  pred <- randomForestSRC::predict.rfsrc(fit, newdata = parts$test)
  list(model = fit,
       predictors = predictors,
       test = parts$test,
       risk_test = as.numeric(pred$predicted))  # ensemble mortality (higher = worse)
}
