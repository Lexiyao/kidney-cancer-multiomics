# Auto-sourced by testthat before every test file.
# Loads project constants + any pure-function modules, and exposes the
# subsampled-fixture directory used by all R tests (see tests/fixtures/).
root <- testthat::test_path("..", "..")

source(file.path(root, "R", "constants.R"))

# This glob already covers every module's function file, INCLUDING Module 3's
# R/functions_sanity.R. The plan's Task 3.2 asks for an explicit
# `source(test_path("..", "..", "R", "functions_sanity.R"))` here, but that
# would source the file a SECOND time and re-create the very duplicate-
# definition churn the plan's own note warns against. The glob is the single
# source of truth; new R/functions_*.R files need no helper edit.
fn_files <- list.files(
  file.path(root, "R"),
  pattern = "^functions_.*\\.R$",
  full.names = TRUE
)
invisible(lapply(fn_files, source))

FIXTURE_DIR <- file.path(root, "tests", "fixtures")

# Loads committed subsampled fixtures for the R test suite.
# NOTE: R/constants.R, all R/functions_*.R, and FIXTURE_DIR are already
# sourced/defined above by the Phase 0 (Task 0.2) portion of this helper.
load_fixture <- function(name) {
  readRDS(testthat::test_path("..", "fixtures", name))
}

# --- Module 2 (integration) test helpers ---
# Named list of features x samples matrices on the shared sample columns,
# in the exact view shape fn_run_mofa / fn_run_snf consume. Mutation is
# deliberately absent: it is an external label, never a MOFA view.
load_view_list <- function() {
  rna    <- readRDS(testthat::test_path("..", "fixtures", "rna_subset.rds"))
  methyl <- readRDS(testthat::test_path("..", "fixtures", "methyl_subset.rds"))
  cnv    <- readRDS(testthat::test_path("..", "fixtures", "cnv_subset.rds"))
  common <- Reduce(intersect, list(colnames(rna), colnames(methyl), colnames(cnv)))
  list(
    RNA         = rna[, common, drop = FALSE],
    Methylation = methyl[, common, drop = FALSE],
    CNV         = cnv[, common, drop = FALSE]
  )
}

# MOFA2 trains through reticulate + the mofapy2 Python module. On a bare host
# neither is present, so MOFA-dependent tests SKIP; they run for real in the
# container (RETICULATE_PYTHON + mofapy2, use_basilisk = FALSE).
skip_if_no_mofapy2 <- function() {
  testthat::skip_if_not(
    requireNamespace("reticulate", quietly = TRUE) &&
      reticulate::py_module_available("mofapy2"),
    "mofapy2 not available in reticulate Python"
  )
}
