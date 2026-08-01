#' Log-transform RSEM upper-quartile normalised expression.
#'
#' Produces "log-transformed normalised expression" (spec 2a). This is NOT
#' DESeq2::vst: KIRC_RNASeq2GeneNorm is already normalised, not raw counts.
#'
#' @param rsem_mat gene x sample matrix of RSEM-normalised values (>= 0).
#' @return new gene x sample matrix of log2(x + 1) values.
fn_log2_normalise_rna <- function(rsem_mat) {
  m <- as.matrix(rsem_mat)
  if (any(m < 0, na.rm = TRUE)) {
    stop("RNA matrix has negative values; expected RSEM-normalised expression")
  }
  log2(m + 1)
}

#' Convert methylation beta values to M-values.
#'
#' @param beta_mat CpG x sample matrix of beta values in [0, 1].
#' @param eps clamp bound keeping beta in [eps, 1 - eps] before the logit.
#' @return new CpG x sample matrix of M-values (log2(beta / (1 - beta))).
fn_beta_to_mvalue <- function(beta_mat, eps = MVALUE_CLAMP_EPS) {
  m <- as.matrix(beta_mat)
  m <- pmin(pmax(m, eps), 1 - eps)
  log2(m / (1 - m))
}

#' Drop SNP-adjacent and sex-chromosome methylation probes.
#'
#' CpGs missing from `annotation` (e.g. the ~1.6k HM27-only probes absent from
#' the 450k manifest) are dropped: an unmatched rowname yields an all-NA
#' annotation row, and an NA logical index would fabricate all-NA rows named
#' `<NA>` rather than filtering.
#'
#' Row subsetting only, so the input may be on any scale (beta or M-value):
#' the pipeline filters beta values first to avoid realising the HDF5-backed
#' HM450 matrix before ~95% of its probes are discarded.
#'
#' @param mval_mat CpG x sample matrix (beta or M-value scale).
#' @param annotation data.frame row-aligned to CpGs with `chr` and `is_snp`.
#' @return new matrix with flagged CpGs removed.
fn_drop_bad_probes <- function(mval_mat, annotation) {
  stopifnot(all(c("chr", "is_snp") %in% colnames(annotation)))
  ann <- annotation[rownames(mval_mat), , drop = FALSE]
  drop <- ann$chr %in% SEX_CHROMOSOMES | ann$is_snp
  drop[is.na(drop)] <- TRUE   # unannotated probes are dropped, never kept as NA rows
  mval_mat[!drop, , drop = FALSE]
}

#' Merge HM27 and HM450 methylation on their common CpGs.
#'
#' @param hm27_mat  CpG x sample M-value matrix (HumanMethylation27).
#' @param hm450_mat CpG x sample M-value matrix (HumanMethylation450).
#' @return new CpG x sample matrix over shared CpGs, all samples column-bound.
fn_merge_methyl_platforms <- function(hm27_mat, hm450_mat) {
  common <- intersect(rownames(hm27_mat), rownames(hm450_mat))
  if (length(common) == 0L) {
    stop("no common CpGs between HM27 and HM450 platforms")
  }
  cbind(hm27_mat[common, , drop = FALSE], hm450_mat[common, , drop = FALSE])
}

#' Prepare GISTIC gene-level thresholded CNV as a Gaussian MOFA view.
#'
#' @param gistic_mat gene x sample matrix of thresholded copy-number (-2..2).
#' @return new numeric matrix with all-complete rows only.
fn_prep_cnv <- function(gistic_mat) {
  m <- matrix(as.numeric(gistic_mat), nrow = nrow(gistic_mat),
              dimnames = dimnames(gistic_mat))
  keep <- rowSums(is.na(m)) == 0L
  m[keep, , drop = FALSE]
}

#' Keep the most variable features (rows) of a matrix.
#'
#' @param mat feature x sample matrix.
#' @param n_top number of highest-variance features to retain.
#' @return new matrix with up to n_top rows, in original row order.
fn_top_variable <- function(mat, n_top) {
  stopifnot(n_top >= 1L)
  m <- as.matrix(mat)
  v <- apply(m, 1L, stats::var, na.rm = TRUE)
  keep <- utils::head(order(v, decreasing = TRUE), min(n_top, nrow(m)))
  m[sort(keep), , drop = FALSE]
}

#' Build a CpG annotation (chr + SNP flag) from the Illumina 450k manifest.
#'
#' The annotation package must be ATTACHED, not merely installed or reached
#' through `::`. `minfi::getAnnotation()` does not read the object it is
#' handed: it re-resolves each annotation table BY NAME off the search list
#' (`get(<table>, envir = as.environment("package:<pkg>"))`), so a plain
#' `pkg::object` call dies with
#'   no item called "package:IlluminaHumanMethylation450kanno.ilmn12.hg19"
#'   on the search list
#' — which is exactly how `tar_make()` failed on the `methyl_anno` target.
#' This is the same attach-vs-namespace trap documented for MultiAssayExperiment
#' in `_targets.R`. Attaching here (rather than trusting the caller) keeps the
#' function correct in EVERY context: the targets pipeline, testthat, and a bare
#' `Rscript -e`. Do NOT "simplify" this back into a single `::` call.
#'
#' @return data.frame with rownames = CpG IDs and columns `chr`, `is_snp`.
fn_load_methyl_annotation <- function() {
  pkg <- METHYL_ANNO_PKG
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("methylation annotation package '", pkg, "' is not installed; ",
         "it is a DESCRIPTION Import required by the HEAVY_PULL run")
  }
  if (!paste0("package:", pkg) %in% search()) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  # getExportedValue(pkg, pkg) is `pkg::pkg`: the package exports its
  # IlluminaMethylationAnnotation object under its own name.
  ann <- minfi::getAnnotation(getExportedValue(pkg, pkg))
  data.frame(
    chr    = as.character(ann$chr),
    is_snp = !is.na(ann$Probe_rs) | !is.na(ann$CpG_rs) | !is.na(ann$SBE_rs),
    row.names = rownames(ann),
    stringsAsFactors = FALSE
  )
}
