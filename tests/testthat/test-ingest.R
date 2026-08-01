# S4 `[[` / `[` dispatch on MultiAssayExperiment requires the package to be
# attached, not merely loaded (the _targets DAG gets this via tar_option_set).
suppressPackageStartupMessages(library(MultiAssayExperiment))

test_that("fn_qc_mae drops non-primary-tumour samples from every assay", {
  # Arrange
  mae <- load_fixture("kirc_mae_subset.rds")
  rna_before <- colnames(fn_experiment(mae, "RNASeq2GeneNorm"))
  expect_true(any(substr(rna_before, 14L, 15L) == "11"))  # fixture has normals

  # Act
  qc <- fn_qc_mae(mae)

  # Assert
  rna_after <- colnames(fn_experiment(qc, "RNASeq2GeneNorm"))
  expect_true(all(substr(rna_after, 14L, 15L) == PRIMARY_TUMOUR_CODE))
  expect_lt(length(rna_after), length(rna_before))
})

test_that("fn_extract_mutation_status returns a sample_id data.frame with 0/1 gene columns", {
  # Arrange
  mae <- load_fixture("kirc_mae_subset.rds")

  # Act
  status <- fn_extract_mutation_status(mae, genes = c("BAP1", "PBRM1"))

  # Assert
  expect_true(is.data.frame(status))
  expect_true("sample_id" %in% names(status))
  expect_setequal(setdiff(names(status), "sample_id"), c("BAP1", "PBRM1"))
  # `%in%` coerces, so type and value have to be asserted separately
  for (g in c("BAP1", "PBRM1")) {
    expect_type(status[[g]], "integer")
    expect_true(all(status[[g]] %in% c(0L, 1L)))
  }
  expect_type(status$sample_id, "character")
  expect_true(all(startsWith(status$sample_id, "TCGA-")))
  expect_true(!anyDuplicated(status$sample_id))
  expect_true(all(nchar(status$sample_id) == PATIENT_BARCODE_LEN))  # harmonised IDs
})

# Hand-built deterministic Mutation MAE (curatedTCGAData experiment naming).
mut_mae <- function(mut) {
  exps <- list(mut)
  names(exps) <- paste0(CURATED_CANCER, "_Mutation-", SNAPSHOT_DATE)
  MultiAssayExperiment(experiments = ExperimentList(exps))
}

test_that("fn_extract_mutation_status emits one column per requested gene, absent genes all-zero", {
  # Arrange: only BAP1/PBRM1 are in the Mutation matrix
  mut <- matrix(c("Missense_Mutation", "Silent"), nrow = 2,
                dimnames = list(c("BAP1", "PBRM1"),
                                "TCGA-B0-0001-01A-01D-1-01"))

  # Act
  status <- fn_extract_mutation_status(mut_mae(mut), genes = DRIVER_GENES)

  # Assert
  expect_named(status, c("sample_id", DRIVER_GENES))
  expect_equal(status$MTOR, 0L)      # absent gene -> present column, all zero
  expect_equal(status$BAP1, 1L)
  expect_equal(status$PBRM1, 0L)
})

test_that("fn_mutation_status_dense returns a full gene panel when no gene matches", {
  # Arrange
  mut <- matrix(rep("Missense_Mutation", 2L), nrow = 1,
                dimnames = list("BAP1", c("S1", "S2")))

  # Act
  out <- fn_mutation_status_dense(mut, c("NOPE1", "NOPE2"))

  # Assert
  expect_equal(dim(out), c(2L, 2L))
  expect_equal(rownames(out), c("NOPE1", "NOPE2"))
  expect_true(all(out == 0L))
})

test_that("fn_extract_mutation_status keeps sample_id as rownames for Module 2", {
  # Arrange: Module 2's fn_annotate_mutation does
  # intersect(rownames(factor_matrix), rownames(mut_annot))
  mut <- matrix(c("Missense_Mutation", "Silent", "Silent", "Nonsense_Mutation"),
                nrow = 2,
                dimnames = list(c("BAP1", "PBRM1"),
                                c("TCGA-B0-0001-01A-01D-1-01",
                                  "TCGA-B0-0002-01A-01D-1-01")))

  # Act
  status <- fn_extract_mutation_status(mut_mae(mut), genes = c("BAP1", "PBRM1"))

  # Assert
  expect_equal(rownames(status), c("TCGA-B0-0001", "TCGA-B0-0002"))
  expect_equal(rownames(status), status$sample_id)
})

test_that("fn_extract_mutation_status scores non-silent variants only", {
  # Arrange: deterministic classes, not the random fixture
  mut <- matrix(c("Missense_Mutation", "Silent", "3'UTR", NA_character_),
                nrow = 2,
                dimnames = list(c("BAP1", "PBRM1"),
                                c("TCGA-B0-0001-01A-01D-1-01",
                                  "TCGA-B0-0002-01A-01D-1-01")))

  # Act
  st <- fn_extract_mutation_status(mut_mae(mut), genes = c("BAP1", "PBRM1"))

  # Assert
  expect_equal(st$BAP1,  c(1L, 0L))   # Missense -> 1, 3'UTR (SILENT_CLASSES) -> 0
  expect_equal(st$PBRM1, c(0L, 0L))   # Silent -> 0, NA -> 0
})

test_that("fn_mutation_status_ragged aggregates non-silent hits by Hugo_Symbol", {
  # Arrange: the real curatedTCGAData Mutation is a RaggedExperiment, so this
  # is the branch that actually runs in production.
  mk <- function(sym, vc, starts) {
    GenomicRanges::GRanges(
      "chr3", IRanges::IRanges(starts, width = 1),
      Hugo_Symbol = sym, Variant_Classification = vc
    )
  }
  re <- RaggedExperiment::RaggedExperiment(
    S1 = mk(c("VHL", "BAP1"), c("Missense_Mutation", "Silent"), c(10, 20)),
    S2 = mk("BAP1", "Nonsense_Mutation", 20),
    S3 = mk("PBRM1", "Silent", 30)
  )

  # Act
  st <- fn_mutation_status_ragged(re, c("VHL", "BAP1", "PBRM1", "MTOR"))

  # Assert
  expect_equal(rownames(st), c("VHL", "BAP1", "PBRM1", "MTOR"))
  expect_equal(unname(st["VHL", ]),   c(1L, 0L, 0L))
  expect_equal(unname(st["BAP1", ]),  c(0L, 1L, 0L))
  expect_equal(unname(st["PBRM1", ]), c(0L, 0L, 0L))
  expect_equal(unname(st["MTOR", ]),  c(0L, 0L, 0L))
})
