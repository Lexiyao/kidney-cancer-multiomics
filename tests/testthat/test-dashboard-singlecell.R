# Task 6.9: the v1.1 single-cell page, its store reader, and its verdict text.
#
# The page has one property no other page has: its targets are CONDITIONALLY
# DECLARED (`run_singlecell`, Task 6.8), so "the store is populated and does not
# contain this target" is the NORMAL release/CI state here, not a broken
# restore. fn_dashboard_read() is right to throw in that situation for every
# Module 0-5 target and would be wrong to throw for these, which is why the
# optional reader exists and is tested against the same three store states.
#
# The rest is the gate text. Both arms of the purity verdict must be
# distinguishable at a glance on the rendered page, so both are asserted here --
# a formatter that prints the same thing either way would let a FAILED gate be
# read as a passing one, which is the exact failure spec section 10 exists to
# prevent.

# The `purity_check` producer shape (fn_subtype_purity_test, Task 6.3).
synthetic_purity_check <- function(is_purity_proxy,
                                   gate_arm = if (isTRUE(is_purity_proxy)) {
                                     "purity"
                                   } else {
                                     "none"
                                   }) {
  list(purity_p = 1.2e-08, immune_p = 3.4e-05, purity_eta2 = 0.31,
       immune_eta2 = 0.22, is_purity_proxy = is_purity_proxy,
       gate_arm = gate_arm, n = 60L)
}

# The `sc_mapping` producer shape (map_bulk_signature, Task 6.7), including
# the per-signature coverage the score table refuses to render without.
synthetic_sc_mapping <- function() {
  list(
    scores_by_celltype = list(
      S1 = list(Tumor_epithelial = 0.42, T_cell = -0.11),
      S2 = list(Tumor_epithelial = -0.07, T_cell = 0.38)
    ),
    coverage = list(
      S1 = list(n_genes_requested = 20L, n_genes_used = 18L,
                low_coverage = FALSE),
      S2 = list(n_genes_requested = 20L, n_genes_used = 20L,
                low_coverage = FALSE)
    ),
    purity_confounded = TRUE,
    status = "exploratory",
    interpretation = "stand-in interpretation string"
  )
}

# --- fn_dashboard_read_optional: absence is expected, breakage is not ---------

test_that("fn_dashboard_read_optional returns NULL when there is no store", {
  empty <- tempfile("nostore")
  dir.create(empty)

  expect_null(fn_dashboard_read_optional("sc_mapping", store = empty))
})

test_that("fn_dashboard_read_optional returns NULL for a target a populated store lacks", {
  # The contrast that justifies a second reader: fn_dashboard_read() ERRORS
  # here (test-dashboard.R asserts that), because for a Module 0-5 target a
  # populated store without it is a broken restore. For a flag-gated Module 6
  # target it is the frozen release built with `run_singlecell: false`.
  store <- build_tiny_store()

  expect_null(fn_dashboard_read_optional("sc_mapping", store = store))
  expect_error(fn_dashboard_read("sc_mapping", store = store), "records no")
})

test_that("fn_dashboard_read_optional reads a target the store does hold", {
  store <- build_tiny_store()

  expect_identical(fn_dashboard_read_optional("alpha", store = store), 42L)
})

test_that("fn_dashboard_read_optional still errors on a target recorded as FAILED", {
  # Absent and failed are different: a target that ran and errored is a real
  # problem upstream, and returning NULL would render it as "not built yet".
  store <- build_tiny_store()

  expect_error(fn_dashboard_read_optional("broken", store = store), "ERRORED")
})

# --- fn_purity_gate_line: the verdict a reader actually sees ------------------

test_that("a FAILED purity gate is stated as failed, in the failing colour", {
  line <- fn_purity_gate_line(synthetic_purity_check(TRUE))

  expect_match(line, "FAILED")
  expect_match(line, "verdict-red", fixed = TRUE)
  expect_no_match(line, "verdict-green", fixed = TRUE)
  # The evidence, not just the verdict.
  expect_match(line, "0.31")
  expect_match(line, "60")
})

test_that("a passed purity gate is stated as passed, and is not called failed", {
  line <- fn_purity_gate_line(synthetic_purity_check(FALSE))

  expect_match(line, "verdict-green", fixed = TRUE)
  expect_no_match(line, "FAILED")
})

test_that("both arms of the gate reach the reader, not just the purity one", {
  # The gate is named for two confounds and fires on either, so a line that
  # shows only the purity arm lets an immune-driven FAILURE be read as though
  # purity had produced it -- and, before the OR verdict landed, let an
  # immune-only proxy be reported as a pass with no immune number in sight.
  line <- fn_purity_gate_line(synthetic_purity_check(TRUE, gate_arm = "immune"))

  expect_match(line, "ImmuneScore", fixed = TRUE)
  expect_match(line, "TumorPurity", fixed = TRUE)
  expect_match(line, "3.4e-05", fixed = TRUE)   # immune p
  expect_match(line, "0.22", fixed = TRUE)      # immune eta-squared
  expect_match(line, "immune", fixed = TRUE)    # which arm fired
})

test_that("the stated decision rule uses the operators the gate actually applies", {
  # A published verdict is quotable only if the rule printed beside it is the
  # rule that produced it. fn_subtype_purity_test compares p with a STRICT `<`
  # and eta-squared with an inclusive `>=`; the line said "p <=".
  # The two source assertions are deliberately keyed on the OPERATOR + CONSTANT
  # and not on the name of the variable being compared. The immune-arm refactor
  # moved the comparison into a local `fires(p_value, eta2)` helper, which left
  # the operators untouched but broke an assertion pinned to `p.value` -- a red
  # test that named no real defect. Matching " < PURITY_PROXY_ALPHA" still
  # falsifies the drift this test exists for: "<=" puts "=" where the space is,
  # so a `p_value <= PURITY_PROXY_ALPHA` regression fails here as it should.
  line <- fn_purity_gate_line(synthetic_purity_check(TRUE))
  src  <- paste(readLines(testthat::test_path("..", "..", "R",
                                              "functions_purity.R")),
                collapse = " ")

  expect_match(src, " < PURITY_PROXY_ALPHA", fixed = TRUE)
  expect_no_match(src, "<= PURITY_PROXY_ALPHA", fixed = TRUE)
  expect_match(src, " >= PURITY_PROXY_ETA2", fixed = TRUE)
  expect_match(line, "p &lt;", fixed = TRUE)
  expect_no_match(line, "p &le;", fixed = TRUE)
  expect_match(line, "&ge;", fixed = TRUE)
})

test_that("fn_purity_gate_line refuses an NA verdict instead of reading it as a pass", {
  # Same rule as `_validate_gate_verdict` on the Python side: the unsafe
  # direction is silently presenting a missing verdict as "gate passed", so it
  # stops rather than defaulting.
  expect_error(fn_purity_gate_line(synthetic_purity_check(NA)), "TRUE or FALSE")
})

test_that("fn_purity_gate_line refuses a purity_check missing its evidence fields", {
  pc <- synthetic_purity_check(TRUE)
  pc$purity_eta2 <- NULL

  expect_error(fn_purity_gate_line(pc), "purity_eta2")
})

# --- fn_sc_score_table: the mapping, reshaped for display ---------------------

test_that("fn_sc_score_table gives one row per cell type and one coverage-labelled column per subtype", {
  tbl <- fn_sc_score_table(synthetic_sc_mapping())

  # The gene counts are IN the column headers, so a score can never be
  # separated from the coverage that produced it.
  expect_equal(names(tbl),
               c("cell_type", "S1 (18/20 genes)", "S2 (20/20 genes)"))
  expect_setequal(tbl$cell_type, c("Tumor_epithelial", "T_cell"))
  expect_equal(tbl[["S1 (18/20 genes)"]][tbl$cell_type == "Tumor_epithelial"],
               0.42)
  expect_equal(tbl[["S2 (20/20 genes)"]][tbl$cell_type == "T_cell"], 0.38)
})

test_that("fn_sc_score_table refuses a mapping with no scores rather than showing an empty table", {
  m <- synthetic_sc_mapping()
  m$scores_by_celltype <- list()

  expect_error(fn_sc_score_table(m), "scores_by_celltype")
})

test_that("fn_sc_score_table refuses a mapping without per-signature coverage", {
  # A score without the gene count behind it prints a 2-of-20-gene value with
  # the same authority as a 20-of-20 one; the table must not be renderable in
  # that state.
  m <- synthetic_sc_mapping()
  m$coverage <- NULL

  expect_error(fn_sc_score_table(m), "coverage")
})

test_that("fn_sc_score_table renders an unmeasurable signature as NaN, not 0", {
  # map_bulk_signature emits NaN when no signature gene is present;
  # a 0.000 in the rendered table would be indistinguishable from a measured
  # absence of enrichment.
  m <- synthetic_sc_mapping()
  m$scores_by_celltype$S1 <- list(Tumor_epithelial = NaN, T_cell = NaN)
  m$coverage$S1 <- list(n_genes_requested = 20L, n_genes_used = 0L,
                        low_coverage = TRUE)

  tbl <- fn_sc_score_table(m)

  expect_true(all(is.nan(tbl[["S1 (0/20 genes)"]])))
  expect_false(any(is.nan(tbl[["S2 (20/20 genes)"]])))
})

# --- The gap sentence: the right cause, not the nearest one ------------------

test_that("the single-cell gap sentence blames the flag, not an unrestored asset", {
  # fn_dashboard_pending() says "restore the targets-store release asset", which
  # is FALSE here: the release store is complete and simply has no single-cell
  # targets, because it was built with run_singlecell off. A true gap explained
  # by a false cause is its own kind of dishonesty (cf. fn_driver_gene_gap).
  gap <- fn_singlecell_gap("the bulk to single-cell mapping")

  expect_match(gap, "run_singlecell", fixed = TRUE)
  expect_no_match(gap, "release asset", fixed = TRUE)
})

# --- The page itself ---------------------------------------------------------

test_that("singlecell.qmd is registered in the dashboard navbar", {
  cfg <- yaml::read_yaml(testthat::test_path("..", "..", "dashboard",
                                             "_quarto.yml"))
  left <- cfg$website$navbar$left

  expect_true(any(vapply(left, function(item) {
    identical(item$href, "singlecell.qmd")
  }, logical(1))))
})

test_that("the landing page's page count matches what _quarto.yml renders", {
  # The pre-publication audit removed exactly this class of self-contradiction
  # and the singlecell branch reintroduced it: index.qmd said 'Five pages' and
  # 'no code for it exists' while the same CI job rendered six pages, one of
  # them singlecell.html. Pin the count on both sides so they cannot drift
  # apart again.
  cfg <- yaml::read_yaml(testthat::test_path("..", "..", "dashboard",
                                             "_quarto.yml"))
  index <- paste(readLines(testthat::test_path("..", "..", "dashboard",
                                               "index.qmd")),
                 collapse = "\n")

  expect_length(cfg$website$navbar$left, 6L)
  expect_match(index, "Six pages", fixed = TRUE)
  expect_no_match(index, "Five pages", fixed = TRUE)
  expect_no_match(index, "NOT BUILT.** Module 6", fixed = TRUE)
})

test_that("singlecell.qmd reads the store, not the manifest, and not from the repo root", {
  lines <- readLines(testthat::test_path("..", "..", "dashboard",
                                         "singlecell.qmd"))
  page <- paste(lines, collapse = "\n")
  # Comment and heading lines are dropped before the forbidden-call check: the
  # page's own chunk comment NAMES these two functions to explain why it does
  # not call them, and a match there would be the opposite of the defect.
  code <- paste(lines[!grepl("^\\s*#", lines)], collapse = "\n")

  # `_quarto.yml` sets `execute-dir: file`, so a page's working directory is
  # dashboard/. The plan's version of this page calls `tar_manifest()` (which
  # needs `_targets.R`, absent from dashboard/) and a bare `tar_read()` (which
  # looks for `_targets` under dashboard/). Both take the whole site render
  # down. Every page reads through the fn_dashboard_read* guards instead.
  expect_match(page, "fn_dashboard_read_optional", fixed = TRUE)
  expect_no_match(code, "tar_manifest(", fixed = TRUE)
  expect_no_match(code, "tar_read(", fixed = TRUE)
})
