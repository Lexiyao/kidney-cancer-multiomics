#' Train a MOFA2 model on Gaussian omics views (mutation is NOT a view).
#'
#' MOFA2 is the MAIN integration method (SNFtool is only a cheap sensitivity
#' analysis; iClusterPlus is cut). Views are RNA / Methylation / CNV, each a
#' features x samples matrix on the shared harmonised cohort. Mutation status
#' stays OUT of the model and enters later as an external annotation label
#' (spec 6a).
#'
#' @param view_list named list of features x samples numeric matrices; the
#'   names become MOFA view names.
#' @param n_factors upper bound on factors; MOFA prunes inactive ones.
#' @param seed training seed (reproducibility).
#' @param maxiter maximum training iterations.
#' @param outfile path for the trained-model HDF5 file.
#' @return a trained `MOFA` S4 object (a new object; inputs are not mutated).
#' @note `scale_views = TRUE` is set explicitly: the three views are on
#'   incomparable measurement scales, so MOFA2's unscaled default would let the
#'   highest-raw-variance view dominate the factors.
fn_run_mofa <- function(view_list,
                        n_factors = MOFA_N_FACTORS,
                        seed = MOFA_SEED,
                        maxiter = MOFA_MAXITER,
                        outfile = tempfile(fileext = ".hdf5")) {
  stopifnot(is.list(view_list), length(view_list) >= 2L,
            !is.null(names(view_list)))
  fn_view_samples(view_list)   # enforce the shared-sample-columns contract

  mofa <- MOFA2::create_mofa(view_list)

  data_opts  <- MOFA2::get_default_data_options(mofa)
  # The three views live on incomparable measurement scales — RNA is log2 RSEM,
  # Methylation is unbounded M-values, CNV is GISTIC integers in [-2, 2] spread
  # over ~5x more features — so the MOFA2 default (scale_views = FALSE) weights
  # each view by its RAW total variance and lets whichever view happens to carry
  # the most drive the factors, and with them the headline per-omics R2 table.
  # Scale so the three views contribute on comparable terms (MOFA2's own
  # guidance for views that differ substantially in scale).
  data_opts$scale_views <- TRUE
  model_opts <- MOFA2::get_default_model_options(mofa)
  model_opts$num_factors <- n_factors
  train_opts <- MOFA2::get_default_training_options(mofa)
  train_opts$convergence_mode <- "fast"
  train_opts$maxiter <- maxiter
  train_opts$seed    <- seed

  mofa <- MOFA2::prepare_mofa(
    object           = mofa,
    data_options     = data_opts,
    model_options    = model_opts,
    training_options = train_opts
  )

  # use_basilisk = FALSE -> use the container's RETICULATE_PYTHON / mofapy2
  # (basilisk external-env resolution is set in the scaffold, Module 0).
  MOFA2::run_mofa(mofa, outfile = outfile, use_basilisk = FALSE, save_data = TRUE)
}

#' Sample x factor matrix from a trained MOFA model (single group).
#'
#' The factor matrix is the low-dimensional representation every downstream
#' module consumes (subtyping, mutation annotation, survival). Rows are
#' samples, columns are `Factor1...FactorN`.
#'
#' @param mofa_model a trained `MOFA` S4 object from `fn_run_mofa`.
#' @return numeric matrix, samples x factors (a new object).
fn_extract_factors <- function(mofa_model) {
  stopifnot(methods::is(mofa_model, "MOFA"))
  f <- MOFA2::get_factors(mofa_model, factors = "all", as.data.frame = FALSE)
  mat <- if (is.list(f)) f[[1L]] else f   # single group -> first element
  as.matrix(mat)
}

#' Per-omics-per-factor variance explained (R2, percent), factors x views.
#'
#' Answers "how much of each omics layer does each factor explain?" — the
#' headline interpretability output of the MOFA integration.
#'
#' @param mofa_model a trained `MOFA` S4 object from `fn_run_mofa`.
#' @return numeric matrix of R2 percentages, factors x views (a new object).
fn_variance_explained <- function(mofa_model) {
  stopifnot(methods::is(mofa_model, "MOFA"))
  ve  <- MOFA2::calculate_variance_explained(mofa_model)
  r2f <- ve$r2_per_factor            # list per group
  mat <- if (is.list(r2f)) r2f[[1L]] else r2f
  as.matrix(mat)                     # rows = factors, cols = views
}

#' Assign MOFA factor subtypes via seeded k-means (returns a new factor vector).
#'
#' Clusters samples in the MOFA latent space into `k` subtypes. `K_SUBTYPES`
#' defaults to 4, matching the four documented KIRC mRNA EXPRESSION subtypes
#' m1-m4 (TCGA Nature 2013) — an expression result, not a methylation one; see
#' METHYL_N_STRATA in R/constants.R. The RNG is snapshotted and restored so
#' seeding never mutates caller state.
#'
#' @param factor_matrix numeric matrix, samples x factors (`fn_extract_factors`).
#' @param k number of subtypes.
#' @param seed k-means seed (reproducibility).
#' @return named `factor` of length nrow(factor_matrix), levels `S1...Sk`,
#'   names = sample IDs (a new object; the input is not mutated).
fn_assign_subtypes <- function(factor_matrix,
                               k = K_SUBTYPES,
                               seed = SUBTYPE_SEED) {
  stopifnot(is.matrix(factor_matrix), k >= 2L, nrow(factor_matrix) >= k)
  old <- fn_rng_snapshot()
  on.exit(fn_rng_restore(old), add = TRUE)
  set.seed(seed)
  km <- stats::kmeans(factor_matrix, centers = k, nstart = 25L)
  labels <- factor(km$cluster, levels = seq_len(k),
                   labels = paste0("S", seq_len(k)))
  stats::setNames(labels, rownames(factor_matrix))
}

# Snapshot/restore the RNG so seeding does not mutate caller state.
# (The plan named these `.Random.seed_snapshot` / `.Random.seed_restore`; that
# name fails this repo's own object_name_linter gate, so they carry the
# project's snake_case `fn_` prefix instead. Behaviour is unchanged.)
fn_rng_snapshot <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv())
  } else {
    NULL
  }
}
fn_rng_restore <- function(snapshot) {
  if (is.null(snapshot)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    # `.Random.seed` is base R's reserved name, not ours to rename.
    assign(".Random.seed", snapshot, envir = globalenv()) # nolint: object_name_linter.
  }
}

#' Keep only features with no missing or infinite value, per view.
#'
#' SNFtool cannot tolerate non-finite input and fails obscurely rather than
#' loudly: `affinityMatrix()` runs `apply(Diff, 2, sort)`, and `sort()` DROPS
#' NAs, so the sorted columns come back ragged, `apply()` degrades from a matrix
#' to a list, and the next subscript loses its dim -- surfacing four calls later
#' as `dim(X) must have a positive length`. Real TCGA 450k M-values carry NAs
#' from failed probes, so this fires on real data (it did: target `snf_clusters`,
#' GitHub Actions run 30716612865).
#'
#' Dropping is deliberate, and confined to SNF. MOFA2 models missing data
#' natively, so the MAIN integration keeps every feature; imputing here would
#' inject invented structure into what is only a sensitivity check. The
#' consequence -- SNF sees fewer features than MOFA -- is a conservative
#' comparison and must be reported, not hidden.
#'
#' The guard is RELATIVE, not an absolute feature floor: a caller is entitled to
#' hand SNF a deliberately small view, and rejecting that would be this
#' function's opinion rather than a defect. What must never pass silently is a
#' view whose features were largely destroyed by missingness -- e.g. one bad
#' sample turning every feature incomplete -- so the check is the surviving
#' FRACTION, plus the mathematical minimum of two features.
#'
#' @param view_list named list of features x samples numeric matrices.
#' @param min_frac minimum fraction of features that must survive per view.
#' @return new list of matrices restricted to complete features.
fn_complete_features <- function(view_list, min_frac = MIN_SNF_COMPLETE_FRAC) {
  lapply(stats::setNames(names(view_list), names(view_list)), function(nm) {
    m <- as.matrix(view_list[[nm]])
    keep <- apply(is.finite(m), 1L, all)
    n_keep <- sum(keep)
    if (n_keep < 2L || n_keep < min_frac * nrow(m)) {
      stop(sprintf(
        paste0("view '%s' retains %d of %d complete features (%.1f%%, floor %.0f%%); ",
               "%d were dropped for non-finite values"),
        nm, n_keep, nrow(m), 100 * n_keep / nrow(m), 100 * min_frac, sum(!keep)
      ))
    }
    if (any(!keep)) {
      message(sprintf("fn_run_snf: dropped %d/%d non-finite features from view '%s'",
                      sum(!keep), nrow(m), nm))
    }
    m[keep, , drop = FALSE]
  })
}

#' SNFtool sensitivity clustering; view_list is features x samples.
#'
#' A CHEAP SECOND OPINION on the MOFA subtypes, not a competing main method:
#' each view is turned into a sample-sample affinity matrix, the affinities are
#' fused by similarity network fusion, and the fused network is spectrally
#' clustered. The resulting partition is only ever compared against
#' `fn_assign_subtypes` via `fn_cluster_concordance`.
#'
#' @param view_list named list of features x samples numeric matrices sharing
#'   the same sample columns (the same views handed to `fn_run_mofa`).
#' @param k SNF K-nearest-neighbours for the affinity matrices.
#' @param alpha SNF affinity hyperparameter (sigma).
#' @param t number of SNF diffusion iterations.
#' @param n_clusters number of spectral clusters to cut.
#' @return named `factor` of length ncol(view_list[[1]]), levels
#'   `C1...Cn_clusters`, names = sample IDs (a new object; inputs unchanged).
fn_run_snf <- function(view_list,
                       k = SNF_K,
                       alpha = SNF_ALPHA,
                       t = SNF_T,
                       n_clusters = K_SUBTYPES) {
  stopifnot(is.list(view_list), length(view_list) >= 2L)
  samples <- fn_view_samples(view_list)
  view_list <- fn_complete_features(view_list)

  affinities <- lapply(view_list, function(m) {
    x <- SNFtool::standardNormalization(t(as.matrix(m)))  # samples x features
    # dist2() returns SQUARED euclidean distance; affinityMatrix() wants a
    # DISTANCE matrix (its own example square-roots dist2 before calling it).
    # Feeding squares is not a rescaling: affinityMatrix builds a Gaussian
    # kernel with a locally-adaptive bandwidth derived from the same input, so
    # the kernel would decay in d^4 instead of d^2 and fuse a different network.
    d <- sqrt(SNFtool::dist2(as.matrix(x), as.matrix(x)))
    SNFtool::affinityMatrix(d, K = k, sigma = alpha)
  })

  fused <- SNFtool::SNF(affinities, K = k, t = t)
  cl    <- SNFtool::spectralClustering(fused, K = n_clusters)

  factor(stats::setNames(paste0("C", cl), samples),
         levels = paste0("C", seq_len(n_clusters)))
}

#' Shared sample IDs of a view list, enforcing identical sample columns.
#'
#' SNF fuses affinity matrices POSITIONALLY and takes its labels from the first
#' view, so views whose columns disagree in content or order silently attach
#' cluster labels to the wrong samples (SNFtool only emits a soft "Dim names not
#' consistent" warning). This turns that silent corruption into a hard stop, and
#' enforces the "sharing the same sample columns" contract both integration
#' entry points document.
#'
#' @param view_list named list of features x samples numeric matrices.
#' @return character vector of the shared sample IDs (the first view's colnames).
fn_view_samples <- function(view_list) {
  samples <- colnames(view_list[[1L]])
  stopifnot(
    !is.null(samples),
    all(vapply(view_list,
               function(m) identical(colnames(m), samples),
               logical(1)))
  )
  samples
}

#' Adjusted Rand index between two named cluster labellings, aligned by name.
#'
#' This is the two-method MOFA-vs-SNF CONCORDANCE check that replaces the old
#' "consensus clustering" framing: it asks how far two independent partitions
#' of the same cohort agree, label-name- and order-independently. Samples are
#' matched on the intersection of the two sample-ID sets.
#'
#' @param labels_a named cluster labels (e.g. `subtypes_mofa`).
#' @param labels_b named cluster labels (e.g. `snf_clusters`).
#' @return `list(ari = <numeric>, contingency = <table>, n = <integer>)`
#'   (a new object; inputs are not mutated).
fn_cluster_concordance <- function(labels_a, labels_b) {
  stopifnot(!is.null(names(labels_a)), !is.null(names(labels_b)))
  common <- intersect(names(labels_a), names(labels_b))
  stopifnot(length(common) >= 2L)
  tab <- table(labels_a[common], labels_b[common])
  list(
    ari         = fn_adjusted_rand_index(tab),
    contingency = tab,
    n           = length(common)
  )
}

#' ARI from a contingency table (Hubert & Arabie 1985).
#'
#' Computed in-house so the concordance check adds no extra dependency.
#'
#' @param contingency a two-way contingency `table` of the two labellings.
#' @return numeric adjusted Rand index (1 = identical partitions, ~0 = chance).
fn_adjusted_rand_index <- function(contingency) {
  choose2 <- function(x) x * (x - 1) / 2
  n     <- sum(contingency)
  index <- sum(choose2(contingency))
  a     <- sum(choose2(rowSums(contingency)))
  b     <- sum(choose2(colSums(contingency)))
  expected  <- a * b / choose2(n)
  max_index <- (a + b) / 2
  if (isTRUE(all.equal(max_index, expected))) return(0)
  (index - expected) / (max_index - expected)
}

#' Test each factor against each driver gene's mutation status (annotation only).
#'
#' Mutation is NEVER a MOFA view (spec 6a): it enters here as an EXTERNAL LABEL
#' used to interpret the latent space, answering "which factor tracks BAP1 /
#' PBRM1 status?". Samples are matched on the intersection of the factor matrix
#' and the mutation subset (n = 417 cases, 413 of them inside the 524-case main
#' cohort), and `MIN_MUT_ANNOT_SAMPLES` guards against a degenerate overlap.
#'
#' @param factor_matrix numeric matrix, samples x factors (`fn_extract_factors`).
#' @param mut_annot data.frame with sample-ID rownames and one 0/1 integer
#'   column per driver gene (`fn_extract_mutation_status`).
#' @param genes character vector of driver genes to test.
#' @return data.frame with columns `gene, factor, p_value, median_wt,
#'   median_mut, q_value`, ordered by `p_value` (a new object). Genes that are
#'   constant across the overlap (0% or 100% mutant) get `NA` p/q values rather
#'   than aborting the pipeline — see `fn_test_factor_gene`.
fn_annotate_mutation <- function(factor_matrix, mut_annot,
                                 genes = DRIVER_GENES) {
  stopifnot(is.matrix(factor_matrix), is.data.frame(mut_annot))
  common <- intersect(rownames(factor_matrix), rownames(mut_annot))
  # A bare stopifnot here reports only "length(common) >= MIN_MUT_ANNOT_SAMPLES
  # is not TRUE", which tells an operator neither the actual overlap nor the
  # likeliest cause (two different ID spaces, e.g. full aliquot barcodes versus
  # harmonised patient IDs).
  if (length(common) < MIN_MUT_ANNOT_SAMPLES) {
    stop("mutation annotation overlaps only ", length(common), " of ",
         nrow(factor_matrix), " factor samples (need >= ",
         MIN_MUT_ANNOT_SAMPLES, "); check that both use harmonised ",
         PATIENT_BARCODE_LEN, "-char patient IDs")
  }

  fm      <- factor_matrix[common, , drop = FALSE]
  mm      <- mut_annot[common, , drop = FALSE]
  missing <- setdiff(genes, colnames(mm))
  genes   <- intersect(genes, colnames(mm))
  if (length(genes) < 1L) {
    stop("none of the requested genes are columns of mut_annot; missing: ",
         paste(missing, collapse = ", "))
  }

  rows <- lapply(genes, function(g) fn_test_factor_gene(fm, as.integer(mm[[g]]), g))
  out  <- do.call(rbind, rows)
  out$q_value <- stats::p.adjust(out$p_value, method = "BH")
  out[order(out$p_value), , drop = FALSE]
}

#' Wilcoxon rank-sum of every factor against one binary mutation vector.
#'
#' @param factor_matrix numeric matrix, samples x factors, already restricted to
#'   the samples shared with the mutation annotation.
#' @param status integer 0/1 mutation status, in the same sample order.
#' @param gene name of the gene `status` refers to.
#' @return data.frame, one row per factor (a new object). If `status` is
#'   constant over the overlap the gene is NOT informative and every column but
#'   `gene`/`factor` is `NA_real_`; no test is run.
fn_test_factor_gene <- function(factor_matrix, status, gene) {
  # A driver gene can legitimately be constant here: Module 1 GUARANTEES that a
  # gene absent from the MAF becomes an all-zero column, and a low-frequency
  # gene (MTOR, KDM5C at ~5%) can be 0% or 100% mutant in a restricted cohort.
  # wilcox.test() would then abort the whole mutation_factor_annot target with
  # "grouping factor must have exactly 2 levels", taking Modules 3-5 with it.
  # NA rows keep the gene visible; p.adjust() and order() both propagate NA.
  if (length(unique(status[!is.na(status)])) < 2L) {
    return(data.frame(
      gene       = gene,
      factor     = colnames(factor_matrix),
      p_value    = NA_real_,
      median_wt  = NA_real_,
      median_mut = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, lapply(colnames(factor_matrix), function(fn) {
    vals <- factor_matrix[, fn]
    wt   <- stats::wilcox.test(vals ~ status)
    data.frame(
      gene       = gene,
      factor     = fn,
      p_value    = wt$p.value,
      median_wt  = stats::median(vals[status == 0L]),
      median_mut = stats::median(vals[status == 1L]),
      stringsAsFactors = FALSE
    )
  }))
}
