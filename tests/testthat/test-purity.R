test_that("fn_purity_from_estimate_score matches the published Yoshihara formula", {
  # Arrange: an ESTIMATE score of 0 gives cos(intercept)
  estimate_score <- 0

  # Act
  purity <- fn_purity_from_estimate_score(estimate_score)

  # Assert
  expect_equal(purity, cos(ESTIMATE_PURITY_INTERCEPT), tolerance = 1e-9)
})

test_that("fn_purity_from_estimate_score is vectorised and monotone-decreasing in score", {
  # Arrange
  scores <- c(-1000, 0, 1000, 3000)

  # Act
  purity <- fn_purity_from_estimate_score(scores)

  # Assert: higher ESTIMATE score => lower tumour purity over this range
  expect_length(purity, 4)
  expect_true(all(diff(purity) < 0))
})

# ADDED BEYOND THE PLAN. The zero crossing of the purity transform was carried
# in prose only, in two files, and the number was wrong by ~63% in the UNSAFE
# direction: it claimed the transform stays positive far above the range TCGA
# ESTIMATEScores actually reach, when in fact it goes negative inside it. Pin
# the boundary as a test so prose cannot drift away from arithmetic again.
test_that("the purity transform crosses zero where the source says it does", {
  # Arrange: the analytic crossing of cos(intercept + slope * S)
  crossing <- (pi / 2 - ESTIMATE_PURITY_INTERCEPT) / ESTIMATE_PURITY_SLOPE

  # Act / Assert: behaviour at and around the boundary
  expect_equal(crossing, 6579.6, tolerance = 1e-4)
  expect_equal(fn_purity_from_estimate_score(crossing), 0, tolerance = 1e-9)
  expect_lt(fn_purity_from_estimate_score(7000), 0)
  expect_gt(fn_purity_from_estimate_score(6000), 0)
})

test_that("the crossing documented in the source matches the computed crossing", {
  # Arrange: every source line that states the crossing in prose
  crossing <- (pi / 2 - ESTIMATE_PURITY_INTERCEPT) / ESTIMATE_PURITY_SLOPE
  # Comment markers are stripped and the file flattened to one string: the
  # sentence is line-wrapped differently in the two files, so a per-line regex
  # would silently miss one of them.
  flatten <- function(path) {
    txt <- readLines(testthat::test_path("..", "..", "R", path))
    paste(gsub("^[[:space:]]*#'?", "", txt), collapse = " ")
  }

  # Act
  claims <- vapply(c("constants.R", "functions_purity.R"), function(path) {
    # perl = TRUE is load-bearing: with the default TRE engine the lazy `.*?`
    # makes the trailing `[0-9]+` match a single digit, so "~6580" reads as 6.
    m <- regmatches(flatten(path), regexpr("crosses zero.*?exceeds ~[0-9]+",
                                           flatten(path), perl = TRUE))
    if (length(m) == 0L) NA_character_ else m
  }, character(1))
  documented <- as.numeric(sub(".*~", "", claims))

  # Assert: the prose is stated in both files, and both agree with arithmetic
  expect_length(documented, 2L)
  expect_true(all(abs(documented - crossing) / crossing < 0.01))
})

test_that("fn_kruskal_eta_squared returns the epsilon-corrected eta-squared", {
  # Arrange: H, k groups, n observations
  h_stat <- 12; n_levels <- 3L; n <- 30L

  # Act
  eta2 <- fn_kruskal_eta_squared(h_stat, n_levels, n)

  # Assert: (H - k + 1) / (n - k) = (12 - 3 + 1) / (30 - 3) = 10/27
  expect_equal(eta2, 10 / 27, tolerance = 1e-9)
})

test_that("fn_subtype_purity_test flags subtypes as a purity proxy when purity tracks subtype", {
  # Arrange: 3 subtypes with strongly separated purity/immune means
  set.seed(1)
  subs <- rep(c("S1", "S2", "S3"), each = 20)
  ids  <- paste0("case", seq_along(subs))
  purity_df <- data.frame(
    sample_id    = ids,
    ImmuneScore  = c(rnorm(20, -1000), rnorm(20, 0), rnorm(20, 1000)),
    TumorPurity  = c(rnorm(20, 0.9, .02), rnorm(20, 0.6, .02), rnorm(20, 0.3, .02)),
    stringsAsFactors = FALSE
  )
  subtypes <- data.frame(sample_id = ids, subtype = subs, stringsAsFactors = FALSE)

  # Act
  res <- fn_subtype_purity_test(purity_df, subtypes)

  # Assert
  expect_true(res$is_purity_proxy)
  expect_lt(res$purity_p, PURITY_PROXY_ALPHA)
  expect_gte(res$purity_eta2, PURITY_PROXY_ETA2)
  expect_equal(res$n, 60L)
})

test_that("fn_subtype_purity_test does NOT flag when purity is independent of subtype", {
  # Arrange: purity drawn from the same distribution for every subtype
  set.seed(2)
  subs <- rep(c("S1", "S2", "S3"), each = 20)
  ids  <- paste0("case", seq_along(subs))
  purity_df <- data.frame(
    sample_id   = ids,
    ImmuneScore = rnorm(60, 0, 500),
    TumorPurity = rnorm(60, 0.6, 0.05),
    stringsAsFactors = FALSE
  )
  subtypes <- data.frame(sample_id = ids, subtype = subs, stringsAsFactors = FALSE)

  # Act
  res <- fn_subtype_purity_test(purity_df, subtypes)

  # Assert
  expect_false(res$is_purity_proxy)
  expect_equal(res$gate_arm, "none")
})

# ADDED BEYOND THE PLAN. The plan's verdict reads the PURITY arm only, while
# both the gate's own name and every sentence written about it claim it covers
# "tumour-purity / immune-infiltration". Immune infiltration is the more likely
# confound of the two in ccRCC, so a cohort whose subtypes are a pure immune
# proxy is exactly the case the gate exists to catch -- and with an
# unconsumed `immune_p` it was the case the gate cleared.
test_that("the gate FIRES on an immune-infiltration proxy even when purity is flat", {
  # Arrange: the real cohort shape (n = 524; S1 20 / S2 306 / S3 76 / S4 122).
  # ImmuneScore separates hard by subtype, TumorPurity does not. The counts are
  # the measured subtype sizes; the SCORES are synthetic stand-ins and are not
  # a measurement of anything.
  set.seed(11)
  sizes <- c(S1 = 20L, S2 = 306L, S3 = 76L, S4 = 122L)
  subs  <- rep(names(sizes), sizes)
  ids   <- paste0("case", seq_along(subs))
  purity_df <- data.frame(
    sample_id   = ids,
    ImmuneScore = rnorm(length(subs),
                        c(S1 = -1500, S2 = 0, S3 = 1200, S4 = 2400)[subs], 300),
    TumorPurity = rnorm(length(subs), 0.6, 0.05),
    stringsAsFactors = FALSE
  )
  subtypes <- data.frame(sample_id = ids, subtype = subs, stringsAsFactors = FALSE)

  # Act
  res <- fn_subtype_purity_test(purity_df, subtypes)

  # Assert: the purity arm is silent, the immune arm carries the verdict
  expect_gt(res$purity_p, PURITY_PROXY_ALPHA)
  expect_lt(res$purity_eta2, PURITY_PROXY_ETA2)
  expect_lt(res$immune_p, PURITY_PROXY_ALPHA)
  expect_gte(res$immune_eta2, PURITY_PROXY_ETA2)
  expect_true(res$is_purity_proxy)
  expect_equal(res$gate_arm, "immune")
  expect_equal(res$n, 524L)
})

test_that("the gate names both arms when both fire", {
  # Arrange: purity AND immune both track subtype
  set.seed(13)
  subs <- rep(c("S1", "S2", "S3"), each = 20)
  ids  <- paste0("case", seq_along(subs))
  purity_df <- data.frame(
    sample_id   = ids,
    ImmuneScore = c(rnorm(20, -1000), rnorm(20, 0), rnorm(20, 1000)),
    TumorPurity = c(rnorm(20, 0.9, .02), rnorm(20, 0.6, .02), rnorm(20, 0.3, .02)),
    stringsAsFactors = FALSE
  )
  subtypes <- data.frame(sample_id = ids, subtype = subs, stringsAsFactors = FALSE)

  # Act
  res <- fn_subtype_purity_test(purity_df, subtypes)

  # Assert
  expect_true(res$is_purity_proxy)
  expect_equal(res$gate_arm, "purity+immune")
})

# ADDED BEYOND THE PLAN, deliberately. Task 6.4's only test is
# `skip_if_not_installed("estimate")`, and `estimate` is an R-Forge package that
# is absent both in CI and on the dev box -- so as the plan stands, EVERY line
# of fn_read_estimate_gct and fn_estimate_purity ships never having been
# executed anywhere. The GCT parser needs no `estimate` install to exercise: it
# reads a text file. This test pins its contract against a synthetic scored GCT
# in the exact shape estimate::estimateScore writes on platform = "illumina"
# (three score rows, no TumorPurity row -- that one is Affymetrix-only, which is
# why fn_estimate_purity derives purity itself).
test_that("fn_read_estimate_gct tidies a scored GCT into one row per sample", {
  # Arrange: a scored GCT exactly as estimateScore(platform = "illumina") writes
  gct <- tempfile(fileext = ".gct")
  ids <- paste0("TCGA-B0-494", 1:3)   # hyphens: the parser must not mangle them
  writeLines(
    c("#1.2",
      "3\t3",
      paste(c("NAME", "Description", ids), collapse = "\t"),
      # decimals, as estimateScore actually writes them -- whole numbers would
      # parse as integer and hide a real coercion difference from the live file
      paste(c("StromalScore",  "StromalScore",  -1000.25, 12.5, 1500.75),
            collapse = "\t"),
      paste(c("ImmuneScore",   "ImmuneScore",    -500.50, 20.25, 2500.125),
            collapse = "\t"),
      paste(c("ESTIMATEScore", "ESTIMATEScore", -1500.75, 32.75, 4000.875),
            collapse = "\t")),
    gct
  )

  # Act
  out <- fn_read_estimate_gct(gct)

  # Assert
  expect_s3_class(out, "data.frame")
  expect_equal(out$sample_id, ids)
  expect_setequal(names(out),
                  c("sample_id", "StromalScore", "ImmuneScore", "ESTIMATEScore"))
  expect_equal(out$ESTIMATEScore, c(-1500.75, 32.75, 4000.875))
  expect_type(out$ESTIMATEScore, "double")
})

test_that("fn_estimate_purity returns one purity row per sample with the expected columns", {
  skip_if_not_installed("estimate")

  # Arrange: tiny genes x samples matrix with recognisable gene symbols.
  # min_coverage = 0 because 12 genes cannot clear the production floor; the
  # floor itself is asserted separately below, without needing `estimate`.
  set.seed(3)
  genes <- c("VHL", "PBRM1", "BAP1", "CD3D", "CD14", "PECAM1",
             "EPCAM", "CA9", "ACTA2", "VWF", "COL1A1", "MS4A1")
  expr <- matrix(rnorm(length(genes) * 6, 8, 2),
                 nrow = length(genes),
                 dimnames = list(genes, paste0("case", 1:6)))

  # Act
  out <- fn_estimate_purity(expr, min_coverage = 0)

  # Assert
  expect_s3_class(out, "data.frame")
  expect_setequal(names(out),
                  c("sample_id", "StromalScore", "ImmuneScore",
                    "ESTIMATEScore", "TumorPurity"))
  expect_equal(nrow(out), 6L)
  expect_true(all(out$TumorPurity >= -1 & out$TumorPurity <= 1))
  # The achieved common-gene coverage travels with the result.
  expect_true(is.numeric(attr(out, "common_gene_coverage")))
})

test_that("fn_estimate_gene_coverage reads the filtered GCT's row count", {
  # Arrange: the header filterCommonGenes writes -- version tag, then dims
  gct <- tempfile(fileext = ".gct")
  writeLines(c("#1.2", "9371\t6", "NAME\tDescription\tcase1"), gct)

  # Act / Assert
  expect_equal(fn_estimate_gene_coverage(gct),
               9371 / ESTIMATE_COMMON_GENES_N)
})

test_that("fn_estimate_gene_coverage refuses an unparseable dims line", {
  gct <- tempfile(fileext = ".gct")
  writeLines(c("#1.2", "not-a-count\t6"), gct)

  expect_error(fn_estimate_gene_coverage(gct), "not parseable")
})

test_that("the coverage floor catches a variance-filtered input matrix", {
  # A top-5000-variance matrix can supply at most 5000/10412 = 0.48 of
  # ESTIMATE's common genes; the floor exists so purity_bulk being re-wired
  # from rna_full back to rna_mat fails loudly instead of silently corrupting
  # the gate verdict (Phase 6 review, finding 2).
  expect_error(
    fn_assert_estimate_coverage(5000 / ESTIMATE_COMMON_GENES_N),
    "rna_full"
  )
  expect_invisible(fn_assert_estimate_coverage(0.9))
  expect_equal(fn_assert_estimate_coverage(0.9), 0.9)
})
