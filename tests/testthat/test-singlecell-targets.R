# Task 6.8 wiring: the flag-gated single-cell targets and the gate they consume.
#
# The plan verifies this with `tar_make(sc_mapping)` on a local dry-run. That
# cannot run here or in CI. `sc_object` needs `reticulate` + `scanpy` + a real
# GSE159115 H5, and `sc_mapping`'s ancestors reach `mae_raw` (curatedTCGAData)
# and `mofa_model` (MOFA2), none of which are installed on a dev box or a CI
# runner. Execution is therefore container-deferred, and what is asserted here
# instead is the WIRING, read out of the committed `_targets.R` through
# `tar_manifest`/`tar_network` -- the same approach as
# tests/testthat/test-purity-targets.R, so these tests exercise the shipped
# pipeline rather than a retyped copy of it.
#
# Two properties carry the whole point of spec section 10 and are asserted
# separately below, because either one alone is satisfiable by a broken gate:
#   1. `purity_check` is an INPUT of `sc_mapping` (the edge), and
#   2. the mapping actually READS `is_purity_proxy` off it (the command text).
# A mapping that takes the edge and ignores the verdict is the defect the spec
# is written against.

# The flag-ON config override (singlecell_config_env) lives in
# helper-fixtures.R, shared with test-purity-targets.R.

SC_TARGET_NAMES <- c("sc_h5_path", "sc_object", "bulk_signature_sets",
                     "sc_mapping")

# Synthetic upstream objects in the documented producer shapes. STAND-INS to
# exercise wiring; not measurements, and none may be quoted anywhere.
synthetic_rna_mat <- function(n_genes = 50L, n_samples = 12L, seed = 7L) {
  set.seed(seed)
  m <- matrix(stats::rnorm(n_genes * n_samples), nrow = n_genes)
  rownames(m) <- sprintf("GENE%03d", seq_len(n_genes))
  colnames(m) <- sprintf("TCGA-B0-%04d", seq_len(n_samples))
  m
}

synthetic_subtypes_df <- function(rna_mat, n_groups = 3L) {
  ids <- colnames(rna_mat)
  data.frame(
    sample_id = ids,
    subtype   = paste0("S", rep_len(seq_len(n_groups), length(ids))),
    stringsAsFactors = FALSE
  )
}

test_that("the single-cell targets are ABSENT when run_singlecell is off", {
  # The CI / frozen-release default. The GSE159115 pull is the one heavy step
  # outside the frozen research core, so with the committed config the targets
  # must not merely be skipped -- they must not be DECLARED, or `tar_make()` on
  # a CI runner would try to build them.
  man <- targets_manifest()

  expect_equal(intersect(man$name, SC_TARGET_NAMES), character(0))
})

test_that("the single-cell targets are declared when run_singlecell is on", {
  man <- targets_manifest(env = singlecell_config_env())

  expect_true(all(SC_TARGET_NAMES %in% man$name))
})

test_that("sc_object is built from the configured H5 path, not a literal", {
  edges <- targets_edges(env = singlecell_config_env())

  expect_equal(edges$from[edges$to == "sc_object"], "sc_h5_path")
})

test_that("sc_h5_path REFUSES a null h5_path instead of inventing a layout", {
  # GEO ships GSE159115 as GSE159115_RAW.tar (per-sample 10x H5s), NOT as a
  # single filtered_feature_bc_matrix.h5 -- the path the config used to carry
  # described a file that does not exist on GEO. The committed config now
  # carries no path at all, and a flag-ON run without the manual download must
  # fail here with instructions, not downstream with a missing-file error.
  man <- targets_manifest(env = singlecell_config_env())
  env <- new.env(parent = globalenv())
  env$config <- list(singlecell = list(h5_path = NULL))

  expect_error(eval_target(man, "sc_h5_path", env), "GSE159115_RAW.tar")

  # And the committed default must not resurrect an invented layout.
  committed <- yaml::read_yaml(testthat::test_path("..", "..", "config",
                                                   "params.yml"))
  expect_null(committed$singlecell$h5_path)
})

test_that("sc_mapping consumes the purity gate as an input", {
  edges <- targets_edges(env = singlecell_config_env())

  expect_setequal(edges$from[edges$to == "sc_mapping"],
                  c("sc_object", "bulk_signature_sets", "purity_check"))
})

test_that("sc_mapping reads the gate VERDICT, not just the gate target", {
  # The edge above is necessary and not sufficient: a command that mentions
  # `purity_check` anywhere gets the edge, including one that passes a hard-
  # coded FALSE to `purity_confounded`. This asserts the verdict field is what
  # reaches the Python mapping, so the re-framing cannot be switched off
  # without this test going red.
  man <- targets_manifest(env = singlecell_config_env())
  cmd <- man$command[man$name == "sc_mapping"]

  expect_match(cmd, "purity_confounded\\s*=\\s*purity_check\\$is_purity_proxy")
})

test_that("bulk_signature_sets is derived from bulk expression and subtypes", {
  edges <- targets_edges(env = singlecell_config_env())

  expect_setequal(edges$from[edges$to == "bulk_signature_sets"],
                  c("rna_mat", "subtypes_df"))
})

test_that("the wired bulk_signature_sets returns one gene set per subtype", {
  man <- targets_manifest(env = singlecell_config_env())
  env <- new.env(parent = globalenv())
  env$rna_mat     <- synthetic_rna_mat()
  env$subtypes_df <- synthetic_subtypes_df(env$rna_mat)

  sets <- eval_target(man, "bulk_signature_sets", env)

  expect_setequal(names(sets), c("S1", "S2", "S3"))
  expect_true(all(vapply(sets, length, integer(1)) == SC_SIGNATURE_N_GENES))
  expect_true(all(unlist(sets) %in% rownames(env$rna_mat)))
})

test_that("bulk_signature_sets ranks by the subtype-vs-rest mean difference", {
  man <- targets_manifest(env = singlecell_config_env())
  env <- new.env(parent = globalenv())
  rna <- synthetic_rna_mat()
  sub <- synthetic_subtypes_df(rna)
  # Plant one unambiguous winner for S1 so the ordering is checkable rather
  # than merely well-shaped.
  s1 <- sub$sample_id[sub$subtype == "S1"]
  rna["GENE042", s1] <- rna["GENE042", s1] + 100
  env$rna_mat     <- rna
  env$subtypes_df <- sub

  sets <- eval_target(man, "bulk_signature_sets", env)

  expect_equal(sets$S1[1], "GENE042")
})

test_that("bulk_signature_sets STOPS when a subtype has no expression columns", {
  # rowMeans over a zero-column matrix is NaN, which would silently yield a
  # signature of arbitrary genes (or NAs) and a mapping scored against it. An
  # empty side of the contrast is a broken input, not a result.
  man <- targets_manifest(env = singlecell_config_env())
  env <- new.env(parent = globalenv())
  env$rna_mat     <- synthetic_rna_mat()
  env$subtypes_df <- rbind(
    synthetic_subtypes_df(env$rna_mat),
    data.frame(sample_id = "TCGA-ZZ-9999", subtype = "S9",
               stringsAsFactors = FALSE)
  )

  expect_error(eval_target(man, "bulk_signature_sets", env), "S9")
})
