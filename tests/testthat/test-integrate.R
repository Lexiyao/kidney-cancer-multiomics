# Module 2 (integrate): MOFA2 is the main integration method, SNFtool is a
# cheap sensitivity analysis, and mutation is an external annotation label —
# never a MOFA view. Tests needing the Python side guard on mofapy2.

test_that("fn_run_mofa returns a trained MOFA object with the requested views", {
  # Arrange
  skip_if_no_mofapy2()
  views <- load_view_list()

  # Act
  model <- fn_run_mofa(views, n_factors = 5L, maxiter = 100L)

  # Assert
  expect_s4_class(model, "MOFA")
  expect_setequal(MOFA2::views_names(model), names(views))
})

test_that("fn_extract_factors and fn_variance_explained expose factors and per-omics R2", {
  # Arrange
  skip_if_no_mofapy2()
  views <- load_view_list()
  model <- fn_run_mofa(views, n_factors = 5L, maxiter = 100L)

  # Act
  factors <- fn_extract_factors(model)
  varexp  <- fn_variance_explained(model)

  # Assert
  expect_true(is.matrix(factors))
  expect_equal(nrow(factors), ncol(views$RNA))          # samples
  expect_setequal(rownames(factors), colnames(views$RNA))
  expect_true(is.matrix(varexp))
  expect_true(all(colnames(varexp) %in% names(views)))  # one column per omics
  expect_true(all(varexp >= 0 & varexp <= 100))         # R2 percentages
  # scale_views = TRUE is set so no single omics layer can dominate: every view
  # must end up with a non-trivial share of explained variance.
  expect_true(all(colSums(varexp) > 0))
})

test_that("fn_assign_subtypes returns k seeded, reproducible clusters over sample IDs", {
  # Arrange: two clearly separated blobs in 2-factor space
  set.seed(1)
  a <- matrix(rnorm(40, mean = -5), ncol = 2)
  b <- matrix(rnorm(40, mean =  5), ncol = 2)
  fm <- rbind(a, b)
  rownames(fm) <- paste0("TCGA-", seq_len(nrow(fm)))
  colnames(fm) <- c("Factor1", "Factor2")

  # Act
  s1 <- fn_assign_subtypes(fm, k = 2L)
  s2 <- fn_assign_subtypes(fm, k = 2L)

  # Assert
  expect_s3_class(s1, "factor")
  expect_equal(nlevels(s1), 2L)
  expect_equal(names(s1), rownames(fm))
  expect_identical(s1, s2)                                  # reproducible
  expect_equal(length(unique(s1[1:20])), 1L)               # blob A one cluster
  expect_equal(length(unique(s1[21:40])), 1L)              # blob B one cluster
  # nlevels() reads the DECLARED levels (always seq_len(k) by construction) and
  # the two block-homogeneity checks above are both satisfied by a constant
  # labelling, so without these a stub that ignores factor_matrix and returns
  # "S1" for every sample would pass. The blobs must land in DIFFERENT clusters.
  expect_equal(length(unique(s1)), 2L)
  expect_false(identical(as.character(s1[1]), as.character(s1[40])))
})

test_that("fn_assign_subtypes seeds without mutating the caller's RNG state", {
  # Arrange
  set.seed(6)
  fm <- rbind(matrix(rnorm(40, mean = -5), ncol = 2),
              matrix(rnorm(40, mean =  5), ncol = 2))
  rownames(fm) <- paste0("TCGA-", seq_len(nrow(fm)))
  colnames(fm) <- c("Factor1", "Factor2")
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  saved <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(fn_rng_restore(saved), add = TRUE)

  # Act / Assert: an existing RNG state survives the seeded k-means untouched.
  set.seed(99)
  before <- get(".Random.seed", envir = globalenv())
  invisible(fn_assign_subtypes(fm, k = 2L))
  expect_identical(get(".Random.seed", envir = globalenv()), before)

  # Act / Assert: with no prior RNG state, none is left behind either.
  rm(".Random.seed", envir = globalenv())
  invisible(fn_assign_subtypes(fm, k = 2L))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("fn_run_snf fuses views and returns spectral clusters over sample IDs", {
  # Arrange: two views, two separable sample blocks (features x samples)
  set.seed(2)
  mk <- function() {
    m <- cbind(matrix(rnorm(200, -3), nrow = 20),
               matrix(rnorm(200,  3), nrow = 20))
    colnames(m) <- paste0("TCGA-", seq_len(ncol(m)))
    m
  }
  views <- list(RNA = mk(), CNV = mk())

  # Act
  clusters <- fn_run_snf(views, k = 5L, n_clusters = 2L)

  # Assert
  expect_s3_class(clusters, "factor")
  expect_equal(length(clusters), 20L)
  expect_equal(names(clusters), colnames(views$RNA))
  expect_equal(length(unique(clusters[1:10])),  1L)
  expect_equal(length(unique(clusters[11:20])), 1L)
  # Both clusters must actually be USED and the two blocks must land in
  # DIFFERENT ones: without these, a stub returning a single constant label for
  # all 20 samples satisfies every assertion above.
  expect_equal(length(unique(clusters)), 2L)
  expect_false(identical(as.character(clusters[1]), as.character(clusters[20])))
})

test_that("fn_run_snf fuses information across views (no single view suffices)", {
  # Arrange: the same two sample blocks in both views, but with a signal so weak
  # that NEITHER view alone recovers them — only the fused network does. Without
  # this, both views carry the identical strong signal and an implementation
  # that silently discards every view after the first still passes.
  set.seed(32)
  mk <- function(mu = 0.35) {
    m <- cbind(matrix(rnorm(200, -mu), nrow = 20),
               matrix(rnorm(200,  mu), nrow = 20))
    colnames(m) <- paste0("TCGA-", seq_len(ncol(m)))
    m
  }
  views <- list(RNA = mk(), CNV = mk())

  # Act
  clusters <- fn_run_snf(views, k = 5L, n_clusters = 2L)

  # Assert: exact recovery of the true 10/10 partition. (Verified: clustering
  # either view's affinity matrix alone misassigns samples and fails this.)
  expect_equal(length(unique(clusters[1:10])),  1L)
  expect_equal(length(unique(clusters[11:20])), 1L)
  expect_false(identical(as.character(clusters[1]), as.character(clusters[20])))
})

test_that("fn_run_snf rejects views whose sample columns disagree", {
  # Arrange: SNF fuses affinities POSITIONALLY and labels from view 1, so a
  # permuted view would silently attach labels to the wrong samples.
  set.seed(5)
  mk <- function() {
    m <- cbind(matrix(rnorm(200, -3), nrow = 20),
               matrix(rnorm(200,  3), nrow = 20))
    colnames(m) <- paste0("TCGA-", seq_len(ncol(m)))
    m
  }
  views <- list(RNA = mk(), CNV = mk())
  permuted <- views
  permuted$CNV <- permuted$CNV[, rev(colnames(permuted$CNV)), drop = FALSE]
  unnamed <- views
  colnames(unnamed$RNA) <- NULL

  # Act / Assert
  expect_error(fn_run_snf(permuted, k = 5L, n_clusters = 2L))
  expect_error(fn_run_snf(unnamed, k = 5L, n_clusters = 2L))
})

test_that("fn_cluster_concordance returns ARI = 1 for identical partitions and aligns by name", {
  # Arrange
  a <- factor(setNames(c("S1", "S1", "S2", "S2"), paste0("s", 1:4)))
  # Same partition by NAME (s1,s2 -> one cluster; s3,s4 -> the other), different
  # labels, and an order that is NOT invariant under positional zipping: pairing
  # a and b by position gives four singleton cells (ARI = -0.5), so an
  # implementation that drops the `[common]` name indexing fails here.
  b <- factor(setNames(c("C1", "C2", "C1", "C2"), c("s1", "s3", "s2", "s4")))

  # Act
  res <- fn_cluster_concordance(a, b)

  # Assert
  expect_named(res, c("ari", "contingency", "n"))
  expect_equal(res$n, 4L)
  expect_equal(res$ari, 1)
})

test_that("fn_cluster_concordance returns ARI near 0 for independent random labels", {
  # Arrange
  set.seed(3)
  ids <- paste0("s", 1:200)
  a <- factor(setNames(sample(c("S1", "S2"), 200, TRUE), ids))
  b <- factor(setNames(sample(c("C1", "C2"), 200, TRUE), ids))

  # Act
  res <- fn_cluster_concordance(a, b)

  # Assert
  expect_lt(abs(res$ari), 0.15)
})

test_that("fn_annotate_mutation flags the factor that tracks a gene's mutation status", {
  # Arrange
  set.seed(4)
  ids <- paste0("TCGA-", 1:100)
  status <- rep(c(0L, 1L), each = 50)                      # 50 wt / 50 mutant
  fm <- cbind(
    Factor1 = rnorm(100),                                  # noise
    Factor2 = status * 3 + rnorm(100, sd = 0.5)            # tracks mutation
  )
  rownames(fm) <- ids
  mut <- data.frame(BAP1 = status, PBRM1 = sample(0:1, 100, TRUE),
                    row.names = ids)
  # Break the incidental order/length coincidence between fm and mut. The join
  # is by ROWNAME; with both frames the same length and in the same order, an
  # implementation that drops the `mut_annot[common, ]` re-alignment is
  # indistinguishable from the correct one, so a real misalignment bug (factor
  # values paired with the wrong patient's mutation status) would ship green.
  mut <- mut[sample(nrow(mut)), , drop = FALSE]
  fm  <- rbind(fm, matrix(rnorm(20), nrow = 10,
                          dimnames = list(paste0("TCGA-ZZ-", 1:10),
                                          colnames(fm))))   # factor-only samples

  # Act
  res <- fn_annotate_mutation(fm, mut, genes = c("BAP1", "PBRM1"))

  # Assert
  expect_named(res, c("gene", "factor", "p_value", "median_wt",
                      "median_mut", "q_value"))
  top <- res[1, ]
  expect_equal(top$gene, "BAP1")
  expect_equal(top$factor, "Factor2")
  expect_lt(top$p_value, 1e-6)
  # median_wt / median_mut carry the entire directional biology of the
  # deliverable ("BAP1-mutant tumours sit HIGH on Factor2"); shape-only checks
  # let a swap of the two definitions invert it everywhere, silently.
  expect_gt(top$median_mut, top$median_wt)
  expect_equal(top$median_wt,  stats::median(fm[ids[status == 0L], "Factor2"]))
  expect_equal(top$median_mut, stats::median(fm[ids[status == 1L], "Factor2"]))
  # q_value must be a real BH correction over the 4 gene x factor tests, not a
  # copy of p_value: Modules 3/5 read this column as an FDR.
  expect_equal(res$q_value, stats::p.adjust(res$p_value, method = "BH"))
  expect_gt(res$q_value[1], res$p_value[1])
})

test_that("fn_annotate_mutation returns NA rows for a gene that is constant in the overlap", {
  # Arrange: Module 1 GUARANTEES a driver gene absent from the MAF becomes an
  # all-zero column (test-ingest.R, "absent genes all-zero"), and a ~5% gene can
  # be 0% mutant in a restricted cohort. wilcox.test() would abort the whole
  # mutation_factor_annot target on it.
  set.seed(8)
  ids <- paste0("TCGA-", 1:60)
  status <- rep(c(0L, 1L), each = 30)
  fm <- cbind(Factor1 = rnorm(60), Factor2 = status * 3 + rnorm(60, sd = 0.5))
  rownames(fm) <- ids
  mut <- data.frame(BAP1 = status,
                    MTOR = rep(0L, 60),          # absent from the MAF
                    KDM5C = rep(1L, 60),         # the 100% branch
                    row.names = ids)

  # Act
  res <- fn_annotate_mutation(fm, mut, genes = c("BAP1", "MTOR", "KDM5C"))

  # Assert
  expect_equal(nrow(res), 6L)                                  # 3 genes x 2 factors
  expect_true(all(is.na(res$p_value[res$gene %in% c("MTOR", "KDM5C")])))
  expect_true(all(is.na(res$q_value[res$gene %in% c("MTOR", "KDM5C")])))
  expect_true(all(!is.na(res$p_value[res$gene == "BAP1"])))
  expect_equal(res$gene[1], "BAP1")                            # NAs sort last
})

test_that("fn_annotate_mutation reports a diagnostic error on a degenerate overlap", {
  # Arrange
  set.seed(9)
  ids <- paste0("TCGA-", 1:60)
  fm <- matrix(rnorm(120), nrow = 60,
               dimnames = list(ids, c("Factor1", "Factor2")))
  mut <- data.frame(BAP1 = rep(c(0L, 1L), each = 30), row.names = ids)

  # Act / Assert: the message must name the actual overlap and the likely cause,
  # not just restate the failed predicate.
  expect_error(fn_annotate_mutation(fm[1:20, , drop = FALSE], mut, genes = "BAP1"),
               "overlaps only 20")
  expect_error(fn_annotate_mutation(fm, mut, genes = "NOTAGENE"),
               "none of the requested genes")
})

test_that("fn_run_snf survives non-finite values, which SNFtool cannot handle", {
  # Arrange: real TCGA 450k M-values carry NAs (failed probes). SNFtool's
  # affinityMatrix() does apply(Diff, 2, sort), and sort() DROPS NAs, so the
  # columns come back ragged, apply() returns a list instead of a matrix, and
  # the next subscript loses its dim -> "dim(X) must have a positive length".
  # MOFA2 handles missing data natively; SNF must be given complete features.
  set.seed(11)
  n <- 60L
  mk <- function(p) matrix(rnorm(p * n), p, n,
                           dimnames = list(paste0("f", seq_len(p)),
                                           paste0("s", seq_len(n))))
  views <- list(RNA = mk(200L), Methylation = mk(200L), CNV = mk(200L))
  clean <- fn_run_snf(views)

  dirty <- views
  dirty$Methylation[1L, 3L] <- NA_real_
  dirty$Methylation[2L, 5L] <- Inf

  # Act
  got <- fn_run_snf(dirty)

  # Assert: it runs, and dropping exactly the two offending features gives the
  # same partition as pre-filtering by hand -- so the guard removes features,
  # never samples, and changes nothing else.
  expect_s3_class(got, "factor")
  expect_equal(names(got), colnames(views$RNA))
  by_hand <- views
  by_hand$Methylation <- by_hand$Methylation[-(1:2), , drop = FALSE]
  expect_equal(as.character(got), as.character(fn_run_snf(by_hand)))
  expect_false(identical(as.character(got), rep(as.character(got)[1], n)))
})

test_that("fn_run_snf errors clearly when a view has too few complete features", {
  set.seed(12)
  n <- 60L
  mk <- function(p) matrix(rnorm(p * n), p, n,
                           dimnames = list(paste0("f", seq_len(p)),
                                           paste0("s", seq_len(n))))
  views <- list(RNA = mk(200L), Methylation = mk(200L))
  views$Methylation[, 1L] <- NA_real_   # one bad sample kills every feature

  expect_error(fn_run_snf(views), "Methylation")
  expect_error(fn_run_snf(views), "retains 0 of 200")
})
