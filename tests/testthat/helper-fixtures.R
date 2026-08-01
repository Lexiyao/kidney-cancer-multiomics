# Auto-sourced by testthat before every test file.
# Loads project constants + any pure-function modules, and exposes the
# subsampled-fixture directory used by all R tests (see tests/fixtures/).
root <- testthat::test_path("..", "..")

source(file.path(root, "R", "constants.R"))

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
