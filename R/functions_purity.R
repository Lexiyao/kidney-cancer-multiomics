# Module 6 confound guard: ESTIMATE purity + subtype-purity proxy test.
# Immutable: every function returns a new object; inputs are never mutated.

#' Convert an ESTIMATE score to tumour purity (Yoshihara et al. 2013).
#'
#' Applies the published Affymetrix-calibrated model
#' `cos(intercept + slope * ESTIMATEScore)`. The coefficients live in
#' `R/constants.R`. NOTE: the transform is not clamped, so the returned value
#' is in [-1, 1], not [0, 1] — `cos` crosses zero once the ESTIMATEScore
#' exceeds ~6580, which is INSIDE the range real high-stromal/high-immune ccRCC
#' samples reach. Callers that need a proportion must handle that themselves;
#' silently clamping here would hide an out-of-calibration input.
#'
#' @param estimate_score numeric vector of ESTIMATEScore values.
#' @return numeric vector of tumour-purity estimates.
fn_purity_from_estimate_score <- function(estimate_score) {
  stopifnot(is.numeric(estimate_score))
  cos(ESTIMATE_PURITY_INTERCEPT + ESTIMATE_PURITY_SLOPE * estimate_score)
}

#' Epsilon-corrected eta-squared effect size for a Kruskal-Wallis test.
#'
#' `(H - k + 1) / (n - k)` (Tomczak & Tomczak 2014). DEVIATION FROM THE PLAN,
#' documentation only: the plan's roxygen claims the return is "in [0, 1]".
#' It is not bounded below by zero — under the null H concentrates around
#' `k - 1`, so any `H < k - 1` yields a NEGATIVE value. That is the correct
#' behaviour (it reads as "no effect") and is left unclamped, but the contract
#' is written honestly here so no caller treats the lower bound as guaranteed.
#'
#' @param h_stat numeric Kruskal-Wallis H statistic.
#' @param n_levels integer number of groups.
#' @param n integer total number of observations.
#' @return numeric eta-squared, upper-bounded by 1 and possibly negative.
fn_kruskal_eta_squared <- function(h_stat, n_levels, n) {
  stopifnot(n > n_levels)
  (h_stat - n_levels + 1) / (n - n_levels)
}

#' Test whether bulk MOFA subtypes are a tumour-purity / immune proxy.
#'
#' THE GATE (spec section 10). Runs on BULK data only and BEFORE any
#' bulk->single-cell signature mapping, so its verdict is available to be
#' consumed rather than rationalised after the fact. `is_purity_proxy = TRUE`
#' means the subtypes largely re-encode tumour purity / immune infiltration and
#' the single-cell mapping must be RE-FRAMED (Task 6.7/6.9), not reported as a
#' cellular basis for the subtypes.
#'
#' DEVIATION FROM THE PLAN, and it widens the gate rather than narrowing it.
#' The plan's verdict reads the PURITY arm alone while still returning
#' `immune_p`, so the immune half of the test was computed, returned, and
#' consumed by nothing — no effect size, no verdict contribution, no line on
#' the dashboard. In ccRCC, immune infiltration is the MORE likely of the two
#' confounds, and a cohort whose subtypes are a pure immune proxy (flat purity,
#' hard immune separation) was cleared as "not explained by purity" and then
#' described as cell-intrinsic. Both arms now carry an epsilon-corrected
#' eta-squared and the verdict is their OR: the gate fires if EITHER purity or
#' immune infiltration is both significant and large. `gate_arm` records which.
#' See the immune-only case in tests/testthat/test-purity.R.
#'
#' @param purity_df data.frame with columns sample_id, ImmuneScore, TumorPurity.
#' @param subtypes data.frame with columns sample_id, subtype.
#' @return list(purity_p, immune_p, purity_eta2, immune_eta2, is_purity_proxy,
#'   gate_arm, n). `gate_arm` is one of "none", "purity", "immune",
#'   "purity+immune".
fn_subtype_purity_test <- function(purity_df, subtypes) {
  stopifnot(all(c("sample_id", "ImmuneScore", "TumorPurity") %in% names(purity_df)))
  stopifnot(all(c("sample_id", "subtype") %in% names(subtypes)))

  merged <- merge(purity_df, subtypes, by = "sample_id")
  merged$subtype <- factor(merged$subtype)
  n_levels <- nlevels(merged$subtype)
  n <- nrow(merged)

  kw_purity <- stats::kruskal.test(TumorPurity ~ subtype, data = merged)
  kw_immune <- stats::kruskal.test(ImmuneScore ~ subtype, data = merged)
  purity_eta2 <- fn_kruskal_eta_squared(unname(kw_purity$statistic), n_levels, n)
  immune_eta2 <- fn_kruskal_eta_squared(unname(kw_immune$statistic), n_levels, n)

  fires <- function(p_value, eta2) {
    (p_value < PURITY_PROXY_ALPHA) && (eta2 >= PURITY_PROXY_ETA2)
  }
  purity_fires <- fires(kw_purity$p.value, purity_eta2)
  immune_fires <- fires(kw_immune$p.value, immune_eta2)
  arms <- c("purity", "immune")[c(purity_fires, immune_fires)]

  list(
    purity_p        = unname(kw_purity$p.value),
    immune_p        = unname(kw_immune$p.value),
    purity_eta2     = purity_eta2,
    immune_eta2     = immune_eta2,
    is_purity_proxy = purity_fires || immune_fires,
    gate_arm        = if (length(arms) == 0L) "none" else paste(arms, collapse = "+"),
    n               = n
  )
}

#' Parse an ESTIMATE output GCT into a tidy per-sample data.frame.
#'
#' GCT layout: line 1 is the version tag, line 2 the dimensions, line 3 the
#' header — hence `skip = 2`. `check.names = FALSE` is load-bearing: TCGA
#' barcodes contain hyphens and would otherwise be silently rewritten to dots,
#' breaking every downstream join on `sample_id`.
#'
#' @param gct_path path to the scored GCT written by estimate::estimateScore.
#' @return data.frame with sample_id + one column per ESTIMATE score row.
fn_read_estimate_gct <- function(gct_path) {
  raw <- utils::read.delim(gct_path, skip = 2, header = TRUE,
                           check.names = FALSE, stringsAsFactors = FALSE)
  score_names <- raw$NAME
  values <- raw[, setdiff(names(raw), c("NAME", "Description")), drop = FALSE]
  tidy <- as.data.frame(t(values), stringsAsFactors = FALSE)
  names(tidy) <- score_names
  data.frame(sample_id = rownames(tidy), tidy,
             row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
}

#' Fraction of ESTIMATE's common-gene background surviving in a filtered GCT.
#'
#' `estimate::filterCommonGenes` writes a GCT whose second line is
#' "<n_rows>\t<n_cols>"; the row count is exactly the number of the package's
#' 10412 common genes the input matrix could supply. Split out of
#' fn_estimate_purity so the arithmetic is testable on a fabricated GCT header
#' without the R-Forge package installed.
#'
#' @param filtered_gct_path path to the GCT written by filterCommonGenes.
#' @return length-1 numeric in [0, 1]: n_genes_kept / ESTIMATE_COMMON_GENES_N.
fn_estimate_gene_coverage <- function(filtered_gct_path) {
  header <- readLines(filtered_gct_path, n = 2)
  if (length(header) < 2L) {
    stop("filtered GCT at ", filtered_gct_path, " has no dimensions line")
  }
  n_kept <- suppressWarnings(as.integer(strsplit(header[2], "\t")[[1]][1]))
  if (is.na(n_kept) || n_kept < 0L) {
    stop("filtered GCT dimensions line is not parseable: ", header[2])
  }
  n_kept / ESTIMATE_COMMON_GENES_N
}

#' Refuse an ESTIMATE input that supplies too little of the common-gene space.
#'
#' The guard finding 2 of the Phase 6 review demanded: ESTIMATE scored on a
#' variance-filtered matrix is not the published quantity, and because the
#' purity gate's verdict derives from those scores, a wrong input here corrupts
#' the gate silently. The floor is calibrated so any top-5000-variance matrix
#' (max coverage 5000/10412 = 0.48) fails loudly; see
#' ESTIMATE_MIN_GENE_COVERAGE in R/constants.R.
#'
#' @param coverage length-1 numeric from fn_estimate_gene_coverage.
#' @param min_coverage the floor below which this stops.
#' @return `coverage`, invisibly, when it clears the floor.
fn_assert_estimate_coverage <- function(coverage,
                                        min_coverage =
                                          ESTIMATE_MIN_GENE_COVERAGE) {
  stopifnot(is.numeric(coverage), length(coverage) == 1L)
  if (coverage < min_coverage) {
    stop(sprintf(
      paste0("only %.1f%% of ESTIMATE's %d common genes survived in the input ",
             "matrix (floor: %.0f%%). ESTIMATE must be scored on the FULL ",
             "gene space (the rna_full target), never on a variance-filtered ",
             "matrix such as rna_mat: the filter drops signature genes and ",
             "changes the within-sample ranks ssGSEA is computed from, so the ",
             "scores would silently stop being the published quantity."),
      100 * coverage, ESTIMATE_COMMON_GENES_N, 100 * min_coverage
    ))
  }
  invisible(coverage)
}

#' Run ESTIMATE on a bulk expression matrix and append tumour purity.
#'
#' UNVERIFIED LOCALLY — `estimate` is an R-Forge package present in neither CI
#' nor the dev box (it IS installed in the project Docker image, where the
#' flag-ON Module 6 run executes), so this function has never been executed;
#' its guarded test skips. Two things it must be checked against on the first
#' container run that has `estimate`:
#'   1. the `Description` column written into `input_f`. `filterCommonGenes`
#'      reads that file with `row.names = 1`, so `NAME` becomes the row names
#'      and `Description` survives as an ordinary column — i.e. it may reach the
#'      filtered GCT as a phantom extra "sample" of gene symbols. The plan's own
#'      `expect_equal(nrow(out), 6L)` is the falsifier: a leaked `Description`
#'      column yields 7 rows, not 6.
#'   2. `platform = "illumina"` emits Stromal/Immune/ESTIMATE scores ONLY. The
#'      Affymetrix-only TumorPurity row is deliberately not requested; purity is
#'      derived here via fn_purity_from_estimate_score, which is an
#'      out-of-platform approximation (see R/constants.R).
#'
#' @param expr_mat genes x samples matrix (rownames = gene symbols) of
#'   log-transformed normalised expression over the FULL gene space — the
#'   `rna_full` target (~20500 genes), NEVER the top-variable-filtered
#'   `rna_mat`. ESTIMATE is a published 141+141-gene ssGSEA panel scored
#'   against a 10412 common-gene background: a data-driven variance filter
#'   both drops signature genes and changes every within-sample rank, so a
#'   filtered input yields scores that are not the published quantity (the
#'   same rule that keeps the Brannon ccA/ccB panel on `rna_full`). Enforced
#'   below via `min_coverage`.
#' @param min_coverage floor on the fraction of ESTIMATE's common genes the
#'   input must supply; the achieved fraction is recorded on the result as
#'   attribute `common_gene_coverage`.
#' @return data.frame(sample_id, StromalScore, ImmuneScore, ESTIMATEScore, TumorPurity).
fn_estimate_purity <- function(expr_mat,
                               min_coverage =
                                 ESTIMATE_MIN_GENE_COVERAGE) {
  stopifnot(!is.null(rownames(expr_mat)), !is.null(colnames(expr_mat)))

  input_f    <- tempfile(fileext = ".txt")
  filtered_f <- tempfile(fileext = ".gct")
  scored_f   <- tempfile(fileext = ".gct")

  input_df <- data.frame(NAME = rownames(expr_mat),
                         Description = rownames(expr_mat),
                         as.data.frame(expr_mat, check.names = FALSE),
                         check.names = FALSE, stringsAsFactors = FALSE)
  utils::write.table(input_df, input_f, sep = "\t",
                     quote = FALSE, row.names = FALSE)

  estimate::filterCommonGenes(input.f = input_f, output.f = filtered_f,
                              id = "GeneSymbol")
  # Refuse a filtered/incomplete gene space BEFORE scoring; record the
  # achieved coverage on the result so the first real run documents it.
  coverage <- fn_estimate_gene_coverage(filtered_f)
  fn_assert_estimate_coverage(coverage, min_coverage)
  # platform = "illumina" -> Stromal/Immune/ESTIMATE scores only; purity is
  # computed via fn_purity_from_estimate_score (Affymetrix-calibrated approx).
  estimate::estimateScore(input.ds = filtered_f, output.ds = scored_f,
                          platform = "illumina")

  scores <- fn_read_estimate_gct(scored_f)
  out <- data.frame(
    sample_id     = scores$sample_id,
    StromalScore  = scores$StromalScore,
    ImmuneScore   = scores$ImmuneScore,
    ESTIMATEScore = scores$ESTIMATEScore,
    TumorPurity   = fn_purity_from_estimate_score(scores$ESTIMATEScore),
    stringsAsFactors = FALSE
  )
  attr(out, "common_gene_coverage") <- coverage
  out
}
