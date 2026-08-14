# Task 6.5 wiring, tested WITHOUT a _targets store and WITHOUT the network.
#
# The plan verifies `purity_bulk` / `purity_check` with `tar_make(purity_check)`
# and an assertion on the built object. That cannot run here or in CI: the
# target's ancestors reach `mae_raw` (curatedTCGAData) and `mofa_model` (MOFA2 +
# reticulate), and `purity_bulk` needs the R-Forge `estimate` package — none of
# which are installed on a dev box or a CI runner. So that check would exist
# nowhere. Following tests/testthat/test-dashboard-targets.R, this file reads the
# ACTUAL committed target expressions out of `_targets.R` via `tar_manifest` and
# evaluates them against synthetic upstream objects in the documented producer
# shapes, so it tests the shipped wiring rather than a retyped copy of it.
#
# DECLARATION IS FLAG-GATED (Phase 6 review, findings 1/13). purity_bulk /
# purity_check are declared behind the same `run_singlecell` conditional as the
# sc_* targets, because `estimate` is installed nowhere the frozen core builds
# (renv.lock, freeze-store.yml, ci.yml) and an unconditional declaration made a
# bare `tar_make()` — the freeze-store/heavy-pull/verify-module2 path that
# produces the LIVE `targets-store` release asset — halt at purity_bulk. So the
# flag-ON manifest (singlecell_config_env, helper-fixtures.R) is where the
# wiring is asserted, and the committed default is asserted NOT to carry the
# targets at all.
#
# The synthetic values below are STAND-INS chosen to exercise shapes and both
# gate verdicts. They are not measurements and none may be quoted anywhere.
# `purity_bulk` itself stays container-deferred: see the note on
# fn_estimate_purity in R/functions_purity.R.

# purity_bulk's documented output shape (Task 6.4), with `mean_purity` per
# subtype under the caller's control so both gate verdicts are reachable.
synthetic_purity_bulk <- function(purity_means, sd = 0.02, n_each = 20L) {
  subs <- rep(paste0("S", seq_along(purity_means)), each = n_each)
  ids  <- sprintf("TCGA-B0-%04d", seq_along(subs))
  estimate_score <- stats::rnorm(length(ids), 0, 500)
  data.frame(
    sample_id     = ids,
    StromalScore  = stats::rnorm(length(ids), 0, 500),
    ImmuneScore   = stats::rnorm(length(ids), 0, 500),
    ESTIMATEScore = estimate_score,
    TumorPurity   = unlist(lapply(purity_means, function(m) {
      stats::rnorm(n_each, m, sd)
    })),
    stringsAsFactors = FALSE
  )
}

# subtypes_mofa's documented shape: a NAMED FACTOR (fn_assign_subtypes).
# Deliberately not a character vector -- letting the committed `subtypes_df`
# target do the conversion is what proves the two contracts actually meet.
synthetic_subtypes_mofa <- function(purity_means, n_each = 20L) {
  subs <- rep(paste0("S", seq_along(purity_means)), each = n_each)
  stats::setNames(factor(subs), sprintf("TCGA-B0-%04d", seq_along(subs)))
}

# Runs the committed subtypes_df + purity_check expressions end to end.
run_gate_targets <- function(purity_means, seed) {
  set.seed(seed)
  man <- targets_manifest(env = singlecell_config_env())
  env <- new.env(parent = globalenv())
  env$subtypes_mofa <- synthetic_subtypes_mofa(purity_means)
  env$purity_bulk   <- synthetic_purity_bulk(purity_means)
  env$subtypes_df   <- eval_target(man, "subtypes_df", env)
  eval_target(man, "purity_check", env)
}

test_that("the gate targets are NOT declared under the committed default config", {
  # The critical property (Phase 6 findings 1/13): a bare `tar_make()` on a
  # host without the R-Forge `estimate` package — every CI runner, the
  # freeze-store job that produces the live `targets-store` release asset —
  # must never have purity_bulk in its DAG, or the whole Modules 0-5 release
  # is blocked by a v1.1 non-blocking increment.
  man <- targets_manifest()

  expect_equal(intersect(man$name, c("purity_bulk", "purity_check")),
               character(0))
})

test_that("the Module 6 gate targets are declared when run_singlecell is on", {
  man <- targets_manifest(env = singlecell_config_env())

  expect_true(all(c("purity_bulk", "purity_check") %in% man$name))
  # The plan's Task 6.5 also asks for a second `subtypes_df` target. `targets`
  # rejects duplicate names, so this is the guard against re-adding it.
  expect_equal(sum(man$name == "subtypes_df"), 1L)
})

test_that("purity_check consumes purity_bulk and subtypes_df", {
  edges <- targets_edges(env = singlecell_config_env())

  expect_setequal(edges$from[edges$to == "purity_check"],
                  c("purity_bulk", "subtypes_df"))
  # rna_full, NOT rna_mat: ESTIMATE is a published panel scored against a
  # 10412 common-gene background, and the top-5000-variance filter both drops
  # signature genes and changes the within-sample ranks (finding 2). The same
  # rule that keeps the ccA/ccB panel on rna_full.
  expect_equal(edges$from[edges$to == "purity_bulk"], "rna_full")
})

test_that("the purity gate is computed on bulk data alone, upstream of any single cell", {
  edges <- targets_edges(env = singlecell_config_env())

  # Spec section 10: the confound check must run BEFORE the bulk->single-cell
  # mapping. A single-cell ancestor would mean the gate could only be reached
  # after the mapping had been seen, which is when a verdict stops being a gate
  # and starts being a rationalisation.
  ancestors <- character(0)
  frontier  <- "purity_check"
  repeat {
    parents <- edges$from[edges$to %in% frontier]
    frontier <- setdiff(parents, ancestors)
    if (!length(frontier)) break
    ancestors <- c(ancestors, frontier)
  }

  expect_true("rna_full" %in% ancestors)
  expect_false(any(grepl("^sc_|singlecell", ancestors)))
})

test_that("the purity_check target returns the documented gate structure", {
  pc <- run_gate_targets(purity_means = c(0.9, 0.6, 0.3), seed = 11)

  expect_type(pc, "list")
  expect_setequal(names(pc),
                  c("purity_p", "immune_p", "purity_eta2", "immune_eta2",
                    "is_purity_proxy", "gate_arm", "n"))
  expect_type(pc$is_purity_proxy, "logical")
  expect_length(pc$is_purity_proxy, 1L)
  expect_equal(pc$n, 60L)
})

test_that("the wired gate FLAGS a proxy when purity tracks subtype", {
  pc <- run_gate_targets(purity_means = c(0.9, 0.6, 0.3), seed = 11)

  expect_true(pc$is_purity_proxy)
})

test_that("the wired gate CLEARS the subtypes when purity is independent of them", {
  # The other half of the falsifiability requirement, asserted through the
  # committed targets and not only through the function: a gate that cannot
  # return FALSE licenses the re-framing unconditionally, which is the same
  # defect as having no gate.
  pc <- run_gate_targets(purity_means = c(0.6, 0.6, 0.6), seed = 12)

  expect_false(pc$is_purity_proxy)
})
