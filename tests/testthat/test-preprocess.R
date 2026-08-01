test_that("fn_log2_normalise_rna applies log2(x+1) to RSEM values", {
  # Arrange
  rsem <- matrix(c(0, 1, 3, 7), nrow = 2,
                 dimnames = list(c("g1", "g2"), c("s1", "s2")))

  # Act
  out <- fn_log2_normalise_rna(rsem)

  # Assert
  expect_equal(out, matrix(c(0, 1, 2, 3), nrow = 2,
                           dimnames = list(c("g1", "g2"), c("s1", "s2"))))
})

test_that("fn_log2_normalise_rna rejects negative (non-RSEM) input", {
  # Arrange
  bad <- matrix(c(-1, 2), nrow = 1)

  # Act / Assert
  expect_error(fn_log2_normalise_rna(bad), "negative")
})

test_that("fn_beta_to_mvalue maps beta to the logit (M-value) scale", {
  # Arrange
  beta <- matrix(c(0.5, 0.8), nrow = 1, dimnames = list("cg1", c("s1", "s2")))

  # Act
  m <- fn_beta_to_mvalue(beta)

  # Assert
  expect_equal(m[1, "s1"], 0, tolerance = 1e-8)             # log2(0.5/0.5)
  expect_equal(m[1, "s2"], log2(0.8 / 0.2), tolerance = 1e-8)
})

test_that("fn_beta_to_mvalue clamps 0 and 1 to avoid +/-Inf", {
  # Arrange
  beta <- matrix(c(0, 1), nrow = 1, dimnames = list("cg1", c("s1", "s2")))

  # Act
  m <- fn_beta_to_mvalue(beta)

  # Assert
  expect_true(all(is.finite(m)))
})

test_that("fn_drop_bad_probes removes sex-chromosome and SNP-adjacent CpGs", {
  # Arrange
  mval <- matrix(1:8, nrow = 4,
                 dimnames = list(c("cgA", "cgB", "cgC", "cgD"), c("s1", "s2")))
  anno <- data.frame(
    chr    = c("chr1", "chrX", "chr7", "chr2"),
    is_snp = c(FALSE,  FALSE,  TRUE,   FALSE),
    row.names = c("cgA", "cgB", "cgC", "cgD")
  )

  # Act
  kept <- fn_drop_bad_probes(mval, anno)

  # Assert
  expect_equal(rownames(kept), c("cgA", "cgD"))
})

test_that("fn_drop_bad_probes drops CpGs missing from the annotation, never NA rows", {
  # Arrange: cgX1/cgX2 are HM27-only probes absent from the 450k manifest
  mval <- matrix(1:8, nrow = 4,
                 dimnames = list(c("cgA", "cgX1", "cgX2", "cgB"), c("s1", "s2")))
  anno <- data.frame(
    chr    = c("chr1", "chr2"),
    is_snp = c(FALSE,  FALSE),
    row.names = c("cgA", "cgB")
  )

  # Act
  kept <- fn_drop_bad_probes(mval, anno)

  # Assert
  expect_equal(rownames(kept), c("cgA", "cgB"))
  expect_false(anyNA(rownames(kept)))
  expect_false(anyNA(kept))
})

test_that("fn_merge_methyl_platforms tolerates probes absent from the annotation", {
  # Arrange: cgU is on both platforms but on neither manifest row
  anno <- data.frame(chr = c("chr1", "chr2"), is_snp = c(FALSE, FALSE),
                     row.names = c("cg1", "cg2"))
  hm27  <- matrix(1:6, nrow = 3, dimnames = list(c("cg1", "cg2", "cgU"), c("a", "b")))
  hm450 <- matrix(7:12, nrow = 3, dimnames = list(c("cg1", "cg2", "cgU"), c("c", "d")))

  # Act
  merged <- fn_merge_methyl_platforms(fn_drop_bad_probes(hm27, anno),
                                      fn_drop_bad_probes(hm450, anno))

  # Assert
  expect_equal(rownames(merged), c("cg1", "cg2"))
  expect_false(anyNA(merged))
})

test_that("fn_merge_methyl_platforms keeps common CpGs and all samples", {
  # Arrange
  hm27  <- matrix(1:6,  nrow = 3, dimnames = list(c("cg1","cg2","cg3"), c("a","b")))
  hm450 <- matrix(7:12, nrow = 3, dimnames = list(c("cg2","cg3","cg4"), c("c","d")))

  # Act
  merged <- fn_merge_methyl_platforms(hm27, hm450)

  # Assert
  expect_equal(rownames(merged), c("cg2", "cg3"))
  expect_equal(colnames(merged), c("a", "b", "c", "d"))
})

test_that("fn_merge_methyl_platforms errors when platforms share no CpGs", {
  # Arrange
  hm27  <- matrix(1, nrow = 1, dimnames = list("cg1", "a"))
  hm450 <- matrix(1, nrow = 1, dimnames = list("cg9", "b"))

  # Act / Assert
  expect_error(fn_merge_methyl_platforms(hm27, hm450), "no common CpGs")
})

test_that("fn_prep_cnv coerces to numeric and drops rows with NA", {
  # Arrange
  cnv <- matrix(c(-2, 0, 1, NA, 2, -1), nrow = 3, byrow = TRUE,
                dimnames = list(c("g1", "g2", "g3"), c("s1", "s2")))

  # Act
  out <- fn_prep_cnv(cnv)

  # Assert
  expect_equal(rownames(out), c("g1", "g3"))
  expect_true(is.numeric(out))
})

test_that("fn_top_variable keeps the highest-variance rows", {
  # Arrange: g2 is constant (var 0), g1 and g3 vary
  mat <- matrix(c(1, 100, 5, 5, 2, 300), nrow = 3, byrow = TRUE,
                dimnames = list(c("g1", "g2", "g3"), c("s1", "s2")))

  # Act
  out <- fn_top_variable(mat, 2L)

  # Assert
  expect_equal(nrow(out), 2L)
  expect_false("g2" %in% rownames(out))
})

test_that("fn_top_variable caps n_top at the number of rows", {
  # Arrange
  mat <- matrix(1:4, nrow = 2, dimnames = list(c("g1", "g2"), c("s1", "s2")))

  # Act / Assert
  expect_equal(nrow(fn_top_variable(mat, 99L)), 2L)
})

test_that("probe filtering and common-CpG restriction commute with the M-value transform", {
  # Arrange: licenses _targets.R doing drop/subset BEFORE fn_beta_to_mvalue,
  # so the HDF5-backed HM450 matrix is never realised at full size.
  set.seed(11)
  beta27  <- matrix(runif(20), nrow = 5,
                    dimnames = list(paste0("cg", 1:5), paste0("s", 1:4)))
  beta450 <- matrix(runif(20), nrow = 5,
                    dimnames = list(paste0("cg", 3:7), paste0("t", 1:4)))
  anno <- data.frame(
    chr    = c("chr1", "chrX", "chr2", "chr3", "chr4", "chr5", "chr6"),
    is_snp = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
    row.names = paste0("cg", 1:7)
  )
  common <- intersect(rownames(beta27), rownames(beta450))

  # Act
  transform_last <- fn_merge_methyl_platforms(
    fn_beta_to_mvalue(fn_drop_bad_probes(beta27[common, , drop = FALSE],  anno)),
    fn_beta_to_mvalue(fn_drop_bad_probes(beta450[common, , drop = FALSE], anno))
  )
  transform_first <- fn_merge_methyl_platforms(
    fn_drop_bad_probes(fn_beta_to_mvalue(beta27),  anno),
    fn_drop_bad_probes(fn_beta_to_mvalue(beta450), anno)
  )

  # Assert
  expect_equal(transform_last, transform_first)
  expect_equal(rownames(transform_last), c("cg4", "cg5"))  # cg3 is SNP-adjacent
})
