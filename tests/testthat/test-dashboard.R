# Module 5 -- dashboard consumption helpers. Offline; no store, no network.
#
# WHY THIS FILE EXISTS AT ALL. Task 5.2's `sanity_table` feeds
# `knitr::kable(san[, c("check", "passed", "detail")])` on the landing page, and
# the plan fills `detail` with `if (!is.null(el$detail)) el$detail else ""`. NO
# element of `sanity_results` carries a `detail` field, so every row would
# render blank -- and a blank detail turns the two results that most need words
# into a bare tick or cross:
#   * BAP1 survival PASSES on `hr > 1` while being UNDERPOWERED (HR 1.584,
#     CI 0.967-2.595, p 0.068). A green tick with no caveat reads as a
#     confirmation, which it is not.
#   * m1-m4 FAILS, and `fn_check_methyl_strata` already computed the sentence
#     explaining that the merge manufactures the structure. Dropping it invites
#     a reader -- or a future maintainer -- to treat the red light as noise.

test_that("fn_sanity_detail surfaces the m1-m4 message verbatim", {
  el <- list(label = "m1-m4", pass = FALSE,
             message = "RED: platform ARI 0.583 exceeds the 0.25 veto",
             platform_ari = 0.583)

  detail <- fn_sanity_detail(el)

  expect_type(detail, "character")
  expect_length(detail, 1L)
  expect_true(grepl("platform ARI 0.583", detail, fixed = TRUE))
})

test_that("fn_sanity_detail flags an underpowered check even though it passed", {
  el <- list(label = "BAP1 worse OS", pass = TRUE, underpowered = TRUE,
             hr = 1.584, p_value = 0.0677, n_events = 138, events_required = 470)

  detail <- fn_sanity_detail(el)

  expect_match(detail, "UNDERPOWERED")
  expect_match(detail, "hr = 1.58")
  expect_match(detail, "events_required = 470")
})

test_that("fn_sanity_detail reports the numbers a check measured", {
  el <- list(label = "subtype vs platform", pass = TRUE, ari = 0.0058, n = 524)

  detail <- fn_sanity_detail(el)

  expect_match(detail, "ari = 0.0058")
  expect_match(detail, "n = 524")
})

test_that("fn_sanity_detail never emits an empty string for a real element", {
  el <- list(label = "some check", pass = TRUE, silhouette = 0.1197)

  expect_true(nzchar(fn_sanity_detail(el)))
})

test_that("fn_sanity_detail says so when a check recorded no scalar evidence", {
  el <- list(label = "evidence-free", pass = TRUE)

  expect_match(fn_sanity_detail(el), "no scalar evidence")
})

test_that("fn_sanity_detail prefers an explicit detail field when present", {
  el <- list(label = "x", pass = TRUE, detail = "hand-written detail", n = 1)

  expect_equal(fn_sanity_detail(el), "hand-written detail")
})

test_that("fn_sanity_detail rejects a non-list element", {
  expect_error(fn_sanity_detail("not a list"))
})

# --- fn_dashboard_read: the store guard the .qmd pages render through ---------
#
# There is no `_targets` store on a dev machine, and in CI there is none until
# the release asset has been restored. The plan's pages read the store with a
# bare `targets::tar_read_raw(name, store = "../_targets")`, which throws — so a
# missing store takes the WHOLE render down and Pages has nothing to deploy.
# This guard mirrors read_pipeline_target() in helper-fixtures.R exactly:
# absent store -> NULL so the page degrades to a stated gap; PRESENT store that
# cannot supply the target -> an error, because that is a broken restore
# masquerading as a fresh clone and must not be papered over.

build_tiny_store <- function() {
  dir <- tempfile("tinystore")
  dir.create(dir)
  old <- setwd(dir)
  on.exit(setwd(old), add = TRUE)
  writeLines(c(
    "library(targets)",
    "list(",
    "  tar_target(alpha, 41L + 1L),",
    "  tar_target(broken, stop('deliberate'), error = 'null')",
    ")"
  ), "_targets.R")
  suppressWarnings(suppressMessages(
    targets::tar_make(callr_function = NULL, reporter = "silent")
  ))
  file.path(dir, "_targets")
}

test_that("fn_dashboard_read returns NULL when the store was never restored", {
  empty <- tempfile("nostore")
  dir.create(empty)

  expect_null(fn_dashboard_read("alpha", store = empty))
})

test_that("fn_dashboard_read reads a target out of a real store", {
  store <- build_tiny_store()

  expect_identical(fn_dashboard_read("alpha", store = store), 42L)
})

test_that("fn_dashboard_read errors on a populated store missing the target", {
  store <- build_tiny_store()

  # A restore that predates a target is a BROKEN restore, not an absent one.
  # Returning NULL here would let a page silently render a gap where a real
  # number exists upstream.
  expect_error(fn_dashboard_read("survival_metrics", store = store),
               "records no")
})

test_that("fn_dashboard_read errors when the target is recorded as errored", {
  store <- build_tiny_store()

  expect_error(fn_dashboard_read("broken", store = store), "ERRORED")
})

# --- fn_factor_platform_labels: no factor is shown without its platform AUC ---
#
# Factor2 explains more methylation variance than anything except Factor1, and
# its platform AUC is 0.888 -- it is largely the HM27/HM450 assay split. A page
# that plots factors without that number invites a reader to take a confounded
# axis at face value, which is the single most likely way this project
# overclaims. The verdict uses the SAME q <= SANITY_MAX_P rule Phase 4 used to
# pick its predictors, so the page and the model cannot disagree.

fp_fixture <- function() {
  data.frame(
    factor = c("Factor1", "Factor2", "Factor6"),
    n = 524L,
    auc = c(0.500, 0.888, 0.735),
    p_value = c(0.993, 1.91e-51, 0),
    q_value = c(0.993, 1.43e-50, 0),
    degenerate = c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
}

test_that("fn_factor_platform_labels calls a significant factor confounded", {
  out <- fn_factor_platform_labels(c("Factor1", "Factor2"), fp_fixture())

  expect_equal(out$verdict[out$factor == "Factor2"], FACTOR_VERDICT_CONFOUNDED)
  expect_equal(out$verdict[out$factor == "Factor1"], FACTOR_VERDICT_CLEAN)
})

test_that("fn_factor_platform_labels prints the AUC in every label", {
  out <- fn_factor_platform_labels(c("Factor1", "Factor2"), fp_fixture())

  expect_true(all(grepl("AUC", out$label)))
  expect_match(out$label[out$factor == "Factor2"], "0.888", fixed = TRUE)
  expect_match(out$label[out$factor == "Factor2"], FACTOR_VERDICT_CONFOUNDED,
               fixed = TRUE)
})

test_that("fn_factor_platform_labels refuses to trust a degenerate AUC row", {
  # Factor6's recorded p = q = 0 underflowed; the row's AUC is unusable, and
  # calling it "clean" or "confounded" would both be claims the data cannot
  # support.
  out <- fn_factor_platform_labels("Factor6", fp_fixture())

  expect_equal(out$verdict, FACTOR_VERDICT_DEGENERATE)
})

test_that("fn_factor_platform_labels marks an unscored factor as unscored", {
  # Never silently defaults to clean: an unscored factor is an unknown, and an
  # unknown shown as clean is exactly the laundering this page must prevent.
  out <- fn_factor_platform_labels("Factor99", fp_fixture())

  expect_equal(out$verdict, FACTOR_VERDICT_UNSCORED)
  expect_true(is.na(out$auc))
})

test_that("fn_factor_platform_labels preserves the requested factor order", {
  out <- fn_factor_platform_labels(c("Factor2", "Factor1"), fp_fixture())

  expect_equal(out$factor, c("Factor2", "Factor1"))
})

test_that("fn_factor_platform_labels rejects a table missing its q column", {
  bad <- fp_fixture()[, c("factor", "auc")]

  expect_error(fn_factor_platform_labels("Factor1", bad))
})

# --- Task 5.6: value boxes and the KM frame -----------------------------------
#
# The four value boxes on dashboard.qmd are the most quotable numbers on the
# site. Three of them (`survival_metrics$cindex$cox`, `bap1_auroc$heldout_auroc`,
# `concordance$ari`) are produced by the container and are NOT available in a
# render without the restored store, so the ONLY two acceptable behaviours are
# "print the measured value" and "say there is no value". The plan's verbatim
# `sprintf("%.3f", survival_metrics$cindex$cox)` has a third: on a NULL it
# yields `character(0)`, which Quarto renders as an EMPTY value box — a blank
# where a reader expects a metric, indistinguishable from a bad number.

test_that("fn_dashboard_valuebox formats a measured scalar", {
  # Arrange
  measured <- 0.6789

  # Act
  vb <- fn_dashboard_valuebox(measured, icon = "activity", colour = "info")

  # Assert
  expect_equal(vb$value, "0.679")
  expect_equal(vb$color, "info")
  expect_equal(vb$icon, "activity")
})

test_that("fn_dashboard_valuebox honours an explicit format", {
  vb <- fn_dashboard_valuebox(524, icon = "people", colour = "primary",
                              fmt = "%d")

  expect_equal(vb$value, "524")
})

test_that("fn_dashboard_valuebox degrades a missing target to a stated gap", {
  # Arrange: exactly what fn_dashboard_read() returns with no store restored.
  absent <- NULL

  # Act
  vb <- fn_dashboard_valuebox(absent, icon = "cpu", colour = "success")

  # Assert: a sentence, not a blank and not a number.
  expect_equal(vb$value, DASHBOARD_PENDING_VALUE)
  expect_equal(vb$color, DASHBOARD_PENDING_COLOUR)
  expect_false(vb$color == "success")
})

test_that("fn_dashboard_valuebox degrades NA and non-scalars too", {
  # An errored or half-restored target can yield NA or a length-0 slot; both
  # must read as absent rather than as `NA` or `character(0)` on the page.
  expect_equal(fn_dashboard_valuebox(NA_real_, "x", "info")$value,
               DASHBOARD_PENDING_VALUE)
  expect_equal(fn_dashboard_valuebox(numeric(0), "x", "info")$value,
               DASHBOARD_PENDING_VALUE)
  expect_equal(fn_dashboard_valuebox(c(1, 2), "x", "info")$value,
               DASHBOARD_PENDING_VALUE)
})

test_that("fn_km_curve_df expands a survfit into a step-plot frame", {
  # Arrange: two subtypes, no censoring subtleties needed for the shape test.
  km_src <- data.frame(
    time    = c(10, 20, 30, 40, 50, 60),
    status  = c(1L, 0L, 1L, 1L, 0L, 1L),
    subtype = c("S1", "S1", "S1", "S2", "S2", "S2"),
    stringsAsFactors = FALSE
  )
  fit <- survival::survfit(survival::Surv(time, status) ~ subtype,
                           data = km_src)

  # Act
  out <- fn_km_curve_df(fit)

  # Assert
  expect_named(out, c("time", "surv", "strata"))
  expect_equal(nrow(out), length(fit$time))
  expect_setequal(unique(out$strata), c("S1", "S2"))
  expect_true(all(out$surv >= 0 & out$surv <= 1))
})

test_that("fn_km_curve_df strips the `subtype=` prefix survfit adds", {
  km_src <- data.frame(time = c(10, 20), status = c(1L, 1L),
                       subtype = c("S1", "S2"), stringsAsFactors = FALSE)
  fit <- survival::survfit(survival::Surv(time, status) ~ subtype,
                           data = km_src)

  out <- fn_km_curve_df(fit)

  expect_false(any(grepl("=", out$strata, fixed = TRUE)))
})

# --- fn_metric_text: the landing-page headline table --------------------------
#
# Task 5.8 builds that table with `sprintf("%.3f", sm$cindex$cox)` inside a
# `data.frame(...)` call. Without the frozen store `sm` is NULL, the sprintf
# returns `character(0)`, and `data.frame()` then dies with
# "arguments imply differing number of rows" -- the landing page, the one page
# a reader reaches first, is the ONLY page in the site that would fail to
# render at all rather than degrade. Every cell therefore goes through this
# helper, which returns a length-1 string in every case.

test_that("fn_metric_text formats a measured value", {
  expect_equal(fn_metric_text(0.74861, "%.3f"), "0.749")
  expect_equal(fn_metric_text(524, "%d"), "524")
})

test_that("fn_metric_text returns the stated gap, never an empty cell", {
  for (absent in list(NULL, NA_real_, numeric(0), c(1, 2))) {
    out <- fn_metric_text(absent, "%.3f")
    expect_type(out, "character")
    expect_length(out, 1L)
    expect_equal(out, DASHBOARD_PENDING_VALUE)
  }
})

test_that("fn_metric_text never invents a value for an absent target", {
  # The failure this guards: a plausible-looking default (0, 0.5, "—") in the
  # cell a reader takes as the headline result.
  expect_false(grepl("[0-9]", fn_metric_text(NULL, "%.3f")))
})

# --- Task 5.6: the driver-gene panel, and the feature filter that can empty it -
#
# `rna_mat` is fn_top_variable(rna_full, N_TOP_GENES) -- the 5000 most-variable
# genes of ~20500. NOTHING guarantees VHL/PBRM1/SETD2/BAP1/MTOR/KDM5C survive
# that filter. Two consequences, both of which used to be unhandled on
# dashboard.qmd:
#   * zero survivors -> facet_wrap() on a 0-row frame ABORTS the render with
#     "Faceting variables must have at least one value.", taking down
#     `quarto render dashboard`, `tar_make(dashboard_site)` and the deploy;
#   * some survivors -> the card silently showed a subset and said nothing
#     about the rest, which is the under-reporting R/functions_dashboard.R
#     exists to prevent.
# Both fire ONLY when the store is present, i.e. on the one path a storeless
# render can never reach.

driver_rna <- function(genes, samples) {
  m <- matrix(seq_len(length(genes) * length(samples)),
              nrow = length(genes), ncol = length(samples))
  dimnames(m) <- list(genes, samples)
  m
}

driver_subtypes <- function(samples) {
  data.frame(sample_id = samples,
             subtype = rep_len(c("S1", "S2"), length(samples)),
             stringsAsFactors = FALSE)
}

test_that("fn_driver_gene_frame keeps every driver gene when all survive", {
  samples <- c("A", "B", "C")
  rna <- driver_rna(DRIVER_GENES, samples)

  gd <- fn_driver_gene_frame(rna, driver_subtypes(samples))

  expect_equal(gd$present, DRIVER_GENES)
  expect_length(gd$dropped, 0L)
  expect_equal(nrow(gd$frame), length(DRIVER_GENES) * length(samples))
  expect_named(gd$frame, c("gene", "expr", "subtype"))
  # Expression must stay aligned to its (gene, sample) cell, not be recycled.
  expect_equal(gd$frame$expr[gd$frame$gene == "BAP1"],
               as.numeric(rna["BAP1", samples]))
})

test_that("fn_driver_gene_frame reports the genes the variance filter dropped", {
  samples <- c("A", "B")
  rna <- driver_rna(c("VHL", "BAP1"), samples)

  gd <- fn_driver_gene_frame(rna, driver_subtypes(samples))

  expect_equal(gd$present, c("VHL", "BAP1"))
  expect_equal(gd$dropped, c("PBRM1", "SETD2", "MTOR", "KDM5C"))
  expect_equal(nrow(gd$frame), 2L * length(samples))
})

test_that("fn_driver_gene_frame returns an empty frame when none survive", {
  samples <- c("A", "B")
  rna <- driver_rna(c("GAPDH", "ACTB"), samples)

  gd <- fn_driver_gene_frame(rna, driver_subtypes(samples))

  expect_length(gd$present, 0L)
  expect_equal(gd$dropped, DRIVER_GENES)
  expect_equal(nrow(gd$frame), 0L)
})

test_that("fn_driver_gene_frame returns an empty frame with no shared samples", {
  rna <- driver_rna(DRIVER_GENES, c("A", "B"))

  gd <- fn_driver_gene_frame(rna, driver_subtypes(c("Y", "Z")))

  expect_equal(gd$present, DRIVER_GENES)
  expect_equal(nrow(gd$frame), 0L)
})

test_that("fn_driver_gene_frame carries each sample's own subtype", {
  samples <- c("A", "B", "C")
  st <- data.frame(sample_id = c("C", "A", "B"),
                   subtype = c("S3", "S1", "S2"), stringsAsFactors = FALSE)

  gd <- fn_driver_gene_frame(driver_rna(c("VHL"), samples), st)

  expect_equal(gd$frame$subtype, c("S1", "S2", "S3"))
})

test_that("fn_driver_gene_gap blames the variance filter, never the store", {
  gap <- fn_driver_gene_gap(list(present = character(0),
                                 dropped = DRIVER_GENES,
                                 frame = data.frame()))

  expect_type(gap, "character")
  expect_length(gap, 1L)
  expect_match(gap, "variance filter")
  # The store IS present on this path, so the standing "restore the release
  # asset" sentence would state a false reason.
  expect_false(grepl("release asset", gap, fixed = TRUE))
})

test_that("fn_driver_gene_gap names the sample overlap when that is the gap", {
  gap <- fn_driver_gene_gap(list(present = DRIVER_GENES,
                                 dropped = character(0),
                                 frame = data.frame()))

  expect_match(gap, "no sample")
})

test_that("fn_driver_gene_note names every dropped gene", {
  note <- fn_driver_gene_note(c("MTOR", "KDM5C"))

  expect_match(note, "MTOR")
  expect_match(note, "KDM5C")
  expect_match(note, "4 of 6|2 of 6")
})

test_that("fn_driver_gene_note is silent when nothing was dropped", {
  expect_null(fn_driver_gene_note(character(0)))
})

# --- Task 5.6 / 5.4: a failing verdict must LOOK failing ----------------------
#
# dashboard/styles.css states the rule: "A failing positive control ... must
# LOOK failing. The m1-m4 anchor is RED and stays red on the site; styling it
# the same as a passing row is how a negative result gets laundered into a
# positive one." index.qmd honoured it and dashboard.qmd did not, because each
# page built its own verdict cell. These two helpers exist so there is one
# implementation to diverge from.

test_that("fn_pass_span marks a failed check red and a passed check green", {
  spans <- fn_pass_span(c(TRUE, FALSE))

  expect_match(spans[1], "verdict-green")
  expect_match(spans[1], "PASS")
  expect_match(spans[2], "verdict-red")
  expect_match(spans[2], "FAIL (RED)", fixed = TRUE)
})

test_that("fn_pass_span never emits a failing row without the red class", {
  # The laundering failure mode, asserted directly.
  expect_false(grepl("verdict-green", fn_pass_span(FALSE), fixed = TRUE))
})

test_that("fn_verdict_span colours every factor verdict the pipeline emits", {
  spans <- fn_verdict_span(c(FACTOR_VERDICT_CLEAN, FACTOR_VERDICT_CONFOUNDED,
                             FACTOR_VERDICT_DEGENERATE,
                             FACTOR_VERDICT_UNSCORED))

  expect_match(spans[1], "verdict-green")
  expect_match(spans[2], "verdict-red")
  # An unusable AUC is not a clean bill of health, so it is not styled as one.
  expect_match(spans[3], "verdict-red")
  expect_match(spans[4], "verdict-pending")
  expect_true(all(grepl(">", spans, fixed = TRUE)))
})

test_that("fn_verdict_span keeps the verdict text the pipeline wrote", {
  span <- fn_verdict_span(FACTOR_VERDICT_CONFOUNDED)

  expect_true(grepl(FACTOR_VERDICT_CONFOUNDED, span, fixed = TRUE))
})

test_that("fn_verdict_span refuses a verdict it has no styling rule for", {
  expect_error(fn_verdict_span("something new"), "no styling rule")
})
