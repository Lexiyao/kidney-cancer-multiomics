#' Truncate TCGA barcodes to a harmonised case/sample ID.
#'
#' @param barcodes character vector of full TCGA aliquot barcodes.
#' @param level "patient" (12 chars) or "sample" (15 chars, keeps type code).
#' @return character vector of the same length; returns a new object.
fn_harmonise_ids <- function(barcodes, level = c("patient", "sample")) {
  level <- match.arg(level)
  n <- if (level == "patient") PATIENT_BARCODE_LEN else SAMPLE_BARCODE_LEN
  substr(as.character(barcodes), 1L, n)
}

#' Complete-case intersection of harmonised IDs across modalities.
#'
#' @param id_lists list of character vectors (one per modality), harmonised.
#' @return character vector of IDs present in every element; new object.
fn_intersect_cases <- function(id_lists) {
  if (!is.list(id_lists) || length(id_lists) < 2L) {
    stop("fn_intersect_cases needs a list of at least two ID vectors")
  }
  Reduce(intersect, lapply(id_lists, unique))
}

#' Align a feature x sample matrix onto the common cohort.
#'
#' Keeps one column per patient (first occurrence), relabels columns to
#' patient IDs, and reorders to `common_ids`. Returns a new matrix.
#'
#' @param mat matrix with full TCGA barcodes as column names.
#' @param common_ids character vector of harmonised patient IDs (target order).
#' @return matrix (features x length(common_ids)).
fn_align_samples <- function(mat, common_ids) {
  m <- as.matrix(mat)
  patient <- fn_harmonise_ids(colnames(m))
  keep <- !duplicated(patient) & patient %in% common_ids
  sub <- m[, keep, drop = FALSE]
  colnames(sub) <- patient[keep]
  sub[, common_ids, drop = FALSE]
}

#' Fetch a curatedTCGAData experiment by assay stub.
#'
#' curatedTCGAData names experiments "<DISEASE>_<stub>-<SNAPSHOT_DATE>"
#' (e.g. "KIRC_RNASeq2GeneNorm-20160128"), and `ExperimentList[[<missing>]]`
#' returns NULL silently, so a bare `mae[["RNASeq2GeneNorm"]]` fails much
#' later inside `assay(NULL)`. This resolves the real name and errors loudly.
#'
#' @param mae MultiAssayExperiment from curatedTCGAData.
#' @param stub assay stub as listed in CURATED_ASSAYS.
#' @return the experiment object (matrix / SummarizedExperiment / RaggedExperiment).
fn_experiment <- function(mae, stub) {
  nm <- paste0(CURATED_CANCER, "_", stub, "-", SNAPSHOT_DATE)
  if (!nm %in% names(mae)) {
    stop("experiment not found in MAE: ", nm,
         " (have: ", paste(names(mae), collapse = ", "), ")")
  }
  mae[[nm]]
}
