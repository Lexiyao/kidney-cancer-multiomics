# Module 1 fills this file with fn_harmonise_ids / fn_align_samples /
# fn_intersect_cases tests. Phase 0 seeds a constants smoke test so the
# testthat harness is provably wired.

suppressPackageStartupMessages(library(MultiAssayExperiment))

test_that("constants define the ccRCC driver-gene panel", {
  # Assert (constants sourced by helper-fixtures.R)
  expect_true("BAP1" %in% DRIVER_GENES)
  expect_setequal(
    DRIVER_GENES,
    c("VHL", "PBRM1", "SETD2", "BAP1", "MTOR", "KDM5C")
  )
})

test_that("EPV cap is 10 and methylation is never HM450-only", {
  expect_identical(EPV_CAP, 10L)
  expect_setequal(METHYL_PLATFORMS, c("HM27", "HM450"))
})

test_that("fn_harmonise_ids truncates barcodes to patient level by default", {
  # Arrange
  barcodes <- c("TCGA-B0-4700-01A-01R-1289-07", "TCGA-B0-4700-11A-01R-1289-07")

  # Act
  ids <- fn_harmonise_ids(barcodes)

  # Assert
  expect_equal(ids, c("TCGA-B0-4700", "TCGA-B0-4700"))
})

test_that("fn_harmonise_ids keeps the sample-type code at sample level", {
  # Arrange
  barcodes <- c("TCGA-B0-4700-01A-01R-1289-07", "TCGA-B0-4700-11A-01R-1289-07")

  # Act
  ids <- fn_harmonise_ids(barcodes, level = "sample")

  # Assert
  expect_equal(ids, c("TCGA-B0-4700-01", "TCGA-B0-4700-11"))
})

test_that("fn_intersect_cases returns cases present in every modality", {
  # Arrange
  rna <- c("TCGA-A", "TCGA-B", "TCGA-C", "TCGA-D")
  met <- c("TCGA-B", "TCGA-C", "TCGA-D", "TCGA-E")
  cnv <- c("TCGA-C", "TCGA-D", "TCGA-F")

  # Act
  common <- fn_intersect_cases(list(rna, met, cnv))

  # Assert
  expect_setequal(common, c("TCGA-C", "TCGA-D"))
})

test_that("fn_intersect_cases throws when given fewer than two modalities", {
  # Arrange / Act / Assert
  expect_error(fn_intersect_cases(list(c("TCGA-A"))), "at least two")
})

test_that("fn_align_samples dedupes by patient and reorders to common_ids", {
  # Arrange: two aliquots of TCGA-B0-4700, columns out of cohort order
  mat <- matrix(1:12, nrow = 3,
                dimnames = list(c("g1", "g2", "g3"),
                                c("TCGA-B0-4700-01A-01R-1289-07",
                                  "TCGA-B0-5075-01A-01R-1289-07",
                                  "TCGA-B0-4700-01B-02R-1289-07",
                                  "TCGA-CZ-9999-01A-01R-1289-07")))
  common <- c("TCGA-B0-5075", "TCGA-B0-4700")

  # Act
  aligned <- fn_align_samples(mat, common)

  # Assert
  expect_equal(colnames(aligned), c("TCGA-B0-5075", "TCGA-B0-4700"))
  expect_equal(ncol(aligned), 2L)
  expect_equal(unname(aligned[, "TCGA-B0-4700"]), c(1, 2, 3))  # first aliquot kept
})

test_that("fn_experiment resolves the curatedTCGAData '<disease>_<stub>-<date>' name", {
  # Arrange
  mae <- load_fixture("kirc_mae_subset.rds")

  # Act
  rna <- fn_experiment(mae, "RNASeq2GeneNorm")

  # Assert
  expect_true(
    paste0(CURATED_CANCER, "_RNASeq2GeneNorm-", SNAPSHOT_DATE) %in% names(mae)
  )
  expect_equal(nrow(rna), nrow(mae[[1L]]))
})

test_that("fn_experiment errors loudly instead of returning NULL for a bad stub", {
  # Arrange
  mae <- load_fixture("kirc_mae_subset.rds")

  # Act / Assert
  expect_error(fn_experiment(mae, "NotAnAssay"), "experiment not found in MAE")
})
