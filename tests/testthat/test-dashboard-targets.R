# Task 5.2 wiring, tested WITHOUT a _targets store and WITHOUT the network.
#
# The plan verifies these three targets by running `tar_make()` on them, which
# needs a full real store (Modules 1-4). There is no such store on a dev machine
# or in CI before the release asset is restored, so that check exists nowhere.
# Instead we read the ACTUAL committed target expressions out of `_targets.R`
# via `tar_manifest` and evaluate them against synthetic upstream objects in the
# documented producer shapes. This tests the shipped code, not a copy of it.
#
# The synthetic values below are STAND-INS, chosen to exercise shapes and edge
# cases. They are not measurements and none of them may be quoted anywhere.

# `targets_manifest`, `targets_edges` and `eval_target` used to be defined here.
# They are generic target-introspection helpers, not Module 5 ones, so Task 6.5
# promoted them to tests/testthat/helper-fixtures.R rather than copy them.

upstream_env <- function() {
  env <- new.env(parent = globalenv())
  env$subtypes_mofa <- stats::setNames(
    factor(c("S1", "S2", "S2", "S3")),
    c("TCGA-AA-0001", "TCGA-AA-0002", "TCGA-AA-0003", "TCGA-AA-0004")
  )
  env$survival_df <- data.frame(
    sample_id = c("TCGA-AA-0001", "TCGA-AA-0002", "TCGA-AA-0009"),
    time = c(100, 200, 300), status = c(1L, 0L, 1L),
    Factor1 = 0, Factor4 = 0, age_years = 60, stage_num = 2,
    platform = "HM27", stringsAsFactors = FALSE
  )
  env$sanity_results <- list(
    mutation_freq = list(label = "mutation frequencies", pass = TRUE, n = 417),
    bap1_survival = list(label = "BAP1 worse OS", pass = TRUE,
                         underpowered = TRUE, hr = 1.584, p_value = 0.0677,
                         events_required = 470),
    methyl_strata = list(label = "m1-m4", pass = FALSE,
                         message = "RED: platform ARI 0.583 exceeds the veto",
                         platform_ari = 0.583)
  )
  env
}

test_that("the Module 5 targets are declared in _targets.R", {
  man <- targets_manifest()

  expect_true(all(c("gdc_live_panel", "subtypes_df", "km_subtype_df",
                    "sanity_table") %in% man$name))
})

test_that("gdc_live_panel depends on no research target", {
  root <- testthat::test_path("..", "..")
  old  <- setwd(root)
  on.exit(setwd(old), add = TRUE)

  edges <- targets::tar_network(targets_only = TRUE)$edges

  # The live panel is the ONLY updating thing in the DAG (spec section 9). An
  # edge from any frozen target would mean a cron run could perturb the
  # research core, and would blur the frozen/live line the site depends on.
  expect_equal(nrow(edges[edges$to == "gdc_live_panel", ]), 0L)
})

test_that("subtypes_df adapts the named factor into a two-column frame", {
  env <- upstream_env()

  sd <- eval_target(targets_manifest(), "subtypes_df", env)

  expect_s3_class(sd, "data.frame")
  expect_named(sd, c("sample_id", "subtype"))
  expect_equal(nrow(sd), 4L)
  expect_type(sd$subtype, "character")   # not a factor, not a factor's integer codes
})

test_that("km_subtype_df inner-joins survival to subtypes", {
  env <- upstream_env()
  man <- targets_manifest()
  env$subtypes_df <- eval_target(man, "subtypes_df", env)

  km <- eval_target(man, "km_subtype_df", env)

  expect_true(all(c("sample_id", "time", "status", "subtype") %in% names(km)))
  # TCGA-AA-0009 has survival but no subtype; TCGA-AA-0003/0004 the reverse.
  expect_equal(nrow(km), 2L)
  expect_false("TCGA-AA-0009" %in% km$sample_id)
})

test_that("sanity_table keeps a failing check failing", {
  env <- upstream_env()

  st <- eval_target(targets_manifest(), "sanity_table", env)

  expect_named(st, c("check", "passed", "detail"))
  expect_type(st$passed, "logical")
  expect_false(st$passed[st$check == "m1-m4"])
  expect_match(st$detail[st$check == "m1-m4"], "RED")
})

test_that("sanity_table never ships a blank detail cell", {
  env <- upstream_env()

  st <- eval_target(targets_manifest(), "sanity_table", env)

  # The plan's `else ""` branch would make this fail on all three rows.
  expect_true(all(nzchar(st$detail)))
})

test_that("sanity_table carries the underpowered caveat on a PASSING check", {
  env <- upstream_env()

  st <- eval_target(targets_manifest(), "sanity_table", env)
  row <- st[st$check == "BAP1 worse OS", ]

  expect_true(row$passed)
  expect_match(row$detail, "UNDERPOWERED")
})

# --- Task 5.9: the dashboard_site render target -------------------------------
#
# The plan verifies this target by RUNNING `tar_make(dashboard_site)`, which
# needs a full Modules 1-4 store plus Quarto plus every render-time package.
# That check therefore runs nowhere: not on a dev machine, not in CI before the
# release asset is restored. What can be checked anywhere, and is what actually
# breaks, is the WIRING: the target must exist, must be a file target pointing
# at the rendered site, and must carry an edge from every upstream target whose
# change should re-render the site. A missing edge does not fail loudly -- it
# silently publishes a stale page next to a fresh one.

test_that("dashboard_site is declared and is a file target", {
  root <- testthat::test_path("..", "..")
  old  <- setwd(root)
  on.exit(setwd(old), add = TRUE)

  man <- targets::tar_manifest(fields = c("command", "format"))
  row <- man[man$name == "dashboard_site", ]

  expect_equal(nrow(row), 1L)
  expect_equal(row$format, "file")
  expect_match(row$command, "quarto_render", fixed = TRUE)
})

test_that("dashboard_site depends on every target the pages render", {
  root <- testthat::test_path("..", "..")
  old  <- setwd(root)
  on.exit(setwd(old), add = TRUE)

  edges <- targets::tar_network(targets_only = TRUE)$edges
  parents <- edges$from[edges$to == "dashboard_site"]

  # Each of these is read by at least one .qmd via fn_dashboard_read(). Without
  # the edge, `targets` considers the site up to date after the upstream target
  # changed, and the published page keeps showing the previous number.
  expect_true(all(c("mofa_factors", "mofa_varexp", "subtypes_df",
                    "km_subtype_df", "concordance", "survival_metrics",
                    "bap1_auroc", "sanity_table", "rna_mat", "cox_fit",
                    "factor_platform", "methyl_platform",
                    "mutation_factor_annot", "gdc_live_panel") %in% parents))
})

test_that("no research target depends on dashboard_site", {
  root <- testthat::test_path("..", "..")
  old  <- setwd(root)
  on.exit(setwd(old), add = TRUE)

  edges <- targets::tar_network(targets_only = TRUE)$edges

  # The site is a leaf. If anything upstream read it, the always-refreshed
  # gdc_live_panel would reach the frozen core through it.
  expect_equal(nrow(edges[edges$from == "dashboard_site", ]), 0L)
})

# --- Task 5.11: the weekly cron must be able to build the panel alone ---------

test_that("gdc_live_panel declares only the packages the query needs", {
  root <- testthat::test_path("..", "..")
  old  <- setwd(root)
  on.exit(setwd(old), add = TRUE)

  man <- targets::tar_manifest(fields = "packages")
  pkgs <- man$packages[man$name == "gdc_live_panel"][[1]]

  # The weekly cron builds THIS TARGET AND NOTHING ELSE on a bare runner. If it
  # inherits the global default it also attaches MultiAssayExperiment, and the
  # live-panel refresh starts failing for a Bioconductor reason unrelated to the
  # panel. Equally, it must not quietly stop declaring httr.
  expect_true("httr" %in% pkgs)
  expect_false("MultiAssayExperiment" %in% pkgs)
})

# --- Task 5.10: publishing stays a decision, not a side effect of a merge -----
#
# The single most consequential choice in Phase 5 -- that deploying the site is
# manual and reserved for the repo owner -- rested entirely on a YAML comment.
# Nothing executable asserted it, while the narrower claims around it (cron.yml's
# "only the live panel rebuilt", every target-wiring claim above) were all made
# falsifiable. Meanwhile the plan's own Task 5.10 literal block contains
#   on:
#     push:
#       branches: [main]
# so an executor replaying the plan, or an agent asked to "match the plan",
# reintroduces an automatic outward-facing publish on every merge to main and
# nothing goes red.

test_that("pages.yml stays workflow_dispatch-only until publication is decided", {
  wf <- yaml::read_yaml(testthat::test_path("..", "..",
                                            ".github/workflows/pages.yml"))

  # YAML 1.1 parses the bare `on:` key as the logical TRUE, hence the name.
  expect_identical(names(wf[["TRUE"]]), "workflow_dispatch")
})

test_that("no workflow deploys Pages on a push", {
  root <- testthat::test_path("..", "..")
  files <- list.files(file.path(root, ".github/workflows"), pattern = "\\.yml$",
                      full.names = TRUE)

  # Guards the same decision from the other side: a NEW workflow that runs
  # actions/deploy-pages on push would publish the site without touching
  # pages.yml at all.
  for (f in files) {
    wf <- yaml::read_yaml(f)
    deploys <- grepl("actions/deploy-pages", paste(readLines(f), collapse = "\n"),
                     fixed = TRUE)
    if (!deploys) next
    expect_false("push" %in% names(wf[["TRUE"]]),
                 info = paste(basename(f), "deploys Pages on a push trigger"))
  }
})
