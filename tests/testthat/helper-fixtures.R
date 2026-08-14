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

# --- Reading frozen targets out of the real `_targets` store ---------------
# MOVED HERE from test-clinical.R, unchanged: it was file-local, so the Module-4
# anchors in test-survival.R could not reuse it and the Phase-4 targets ended up
# with no assertion anywhere.
#
# Skipping is reserved for the ONE honest case — a store with no pipeline
# metadata at all (a fresh clone, or the dev machine, where Modules 1-2 never
# run). A populated store that lacks the target, or records it as errored, is a
# FAILURE: an anchor that stays green by silently skipping is the false
# confidence this suite exists to prevent.
read_pipeline_target <- function(name) {
  store <- testthat::test_path("..", "..", "_targets")
  if (!file.exists(file.path(store, "meta", "meta"))) {
    return(NULL)
  }
  meta <- targets::tar_meta(store = store)
  if (!name %in% meta$name) {
    stop("the _targets store at ", store, " has pipeline metadata but records ",
         "no `", name, "` target: the restore predates the platform covariate, ",
         "or an upstream target failed before it.")
  }
  err <- meta$error[meta$name == name]
  if (!is.na(err[[1]])) {
    stop(name, " is recorded in the _targets store but ERRORED: ", err[[1]])
  }
  targets::tar_read_raw(name, store = store)
}

# --- Reading target EXPRESSIONS out of the committed `_targets.R` -----------
# PROMOTED HERE from test-dashboard-targets.R, unchanged. They were file-local,
# so the Module 6 wiring tests (Task 6.5) could not reuse them and would have
# had to copy them a second time. `_targets.R` sources `R/` by relative path,
# so both helpers must run from the repo root, not from tests/testthat.
#
# The point of reading the manifest instead of retyping the target body: these
# tests then exercise the SHIPPED wiring. A target that was edited in
# `_targets.R` but not here fails, rather than passing against a stale copy.
# `fields` is a LITERAL, not a parameter: tar_manifest selects with tidyselect,
# and forwarding a variable there raises the "external vector in selections"
# deprecation on every call. Nothing needs a second field.
#
# `env` exists for Module 6 (Task 6.8). The single-cell targets are declared
# inside `if (isTRUE(config$singlecell$run_singlecell))`, and `_targets.R` reads
# that config at PARSE time -- there is no argument to hand it a different flag.
# So the flag-ON branch is reachable only by overriding the config path in the
# environment for the duration of the call, which is exactly what the plan's own
# Task 6.8 Step 3 dry-run does (`R_CONFIG_FILE=config/params.local.yml`).
# `tar_manifest`/`tar_network` run in a `callr` child that inherits this
# process's environment, so `Sys.setenv` here is visible to the parse.
in_repo_root <- function(expr, env = NULL) {
  # `env` is forced BEFORE the setwd: a caller that builds it from
  # `testthat::test_path()` resolves those paths relative to tests/testthat, and
  # forcing it after the move would look them up under the repo root instead.
  force(env)
  old <- setwd(testthat::test_path("..", ".."))
  on.exit(setwd(old), add = TRUE)
  if (length(env) > 0L) {
    previous <- Sys.getenv(names(env), unset = NA_character_, names = TRUE)
    do.call(Sys.setenv, as.list(env))
    on.exit(restore_env(previous), add = TRUE)
  }
  expr
}

# Puts the environment back exactly as it was: variables that were unset before
# are UNSET again, not set to "". A leaked R_CONFIG_FILE would silently flip
# every later test in the session onto the single-cell config.
restore_env <- function(previous) {
  had <- previous[!is.na(previous)]
  if (length(had) > 0L) do.call(Sys.setenv, as.list(had))
  lacked <- names(previous)[is.na(previous)]
  if (length(lacked) > 0L) Sys.unsetenv(lacked)
}

# A params file with the Module 6 flag flipped ON, exposed the way `_targets.R`
# accepts a local override. Nothing here touches the committed config.
# PROMOTED from test-singlecell-targets.R: test-purity-targets.R needs the same
# override now that the ESTIMATE gate targets are declared behind the same flag,
# and testthat gives each test FILE its own environment, so a second copy would
# be the only alternative.
singlecell_config_env <- function(h5_path = "tests/fixtures/gse159115_subset.h5") {
  cfg <- yaml::read_yaml(testthat::test_path("..", "..", "config", "params.yml"))
  cfg$singlecell$run_singlecell <- TRUE
  cfg$singlecell$h5_path <- h5_path
  path <- tempfile("params-singlecell-", fileext = ".yml")
  yaml::write_yaml(cfg, path)
  c(R_CONFIG_FILE = path)
}

targets_manifest <- function(env = NULL) {
  in_repo_root(targets::tar_manifest(fields = "command"), env)
}

targets_edges <- function(env = NULL) {
  in_repo_root(targets::tar_network(targets_only = TRUE)$edges, env)
}

# --- A real, tiny `_targets` store for the dashboard reader guards ------------
# PROMOTED HERE from test-dashboard.R, unchanged. The store guards have three
# states to distinguish -- absent store, populated store missing the target,
# populated store recording an ERRORED target -- and the Module 6 page reader
# (fn_dashboard_read_optional) must be asserted against all three too.
# `alpha` builds; `broken` is recorded with an error rather than absent.
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

eval_target <- function(man, name, env) {
  cmd <- man$command[man$name == name]
  testthat::expect_length(cmd, 1L)
  eval(parse(text = cmd), envir = env)
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
