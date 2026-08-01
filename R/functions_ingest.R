#' Load the versioned KIRC MultiAssayExperiment (heavy pull; local only).
#'
#' @return MultiAssayExperiment from curatedTCGAData 1.34.0 (snapshot 20160128).
fn_load_mae <- function() {
  curatedTCGAData::curatedTCGAData(
    diseaseCode = CURATED_CANCER,
    assays      = CURATED_ASSAYS,
    version     = CURATED_VERSION,
    dry.run     = FALSE
  )
}

#' QC a KIRC MultiAssayExperiment: keep primary-tumour columns only.
#'
#' @param mae MultiAssayExperiment.
#' @return new MultiAssayExperiment restricted to primary-tumour aliquots.
fn_qc_mae <- function(mae) {
  stopifnot(methods::is(mae, "MultiAssayExperiment"))
  keep <- lapply(
    MultiAssayExperiment::colnames(mae),
    function(bc) bc[substr(bc, 14L, 15L) == PRIMARY_TUMOUR_CODE]
  )
  mae[, keep, ]
}

#' Non-silent gene x sample 0/1 status from a synthetic (dense) Mutation matrix.
#'
#' Mirrors fn_mutation_status_ragged's guarantee: one row per requested gene,
#' in `genes` order, genes absent from `mut` filled with 0.
#'
#' @param mut gene x sample character matrix of variant classifications.
#' @param genes driver genes to report.
#' @return integer gene x sample matrix (rows = genes, in `genes` order).
fn_mutation_status_dense <- function(mut, genes) {
  gxs <- if (methods::is(mut, "matrix")) {
    as.matrix(mut)
  } else {
    as.matrix(SummarizedExperiment::assay(mut))
  }
  status <- matrix(0L, nrow = length(genes), ncol = ncol(gxs),
                   dimnames = list(genes, colnames(gxs)))
  present <- intersect(genes, rownames(gxs))
  if (length(present)) {
    sub <- gxs[present, , drop = FALSE]
    nonsilent <- !is.na(sub) & !(sub %in% SILENT_CLASSES)
    status[present, ] <- matrix(as.integer(nonsilent), nrow = length(present),
                                ncol = ncol(gxs))
  }
  status
}

#' Non-silent gene x sample 0/1 status from a real RaggedExperiment Mutation.
#'
#' Aggregates to gene level via the Hugo_Symbol metadata column; a sample is
#' mutant for a gene if any of its ranges hit that gene with a non-silent
#' Variant_Classification. Does NOT rely on compactAssay row names being genes.
#'
#' @param mut RaggedExperiment with Hugo_Symbol + Variant_Classification cols.
#' @param genes driver genes to report.
#' @return integer gene x sample matrix (rows = genes, in `genes` order).
fn_mutation_status_ragged <- function(mut, genes) {
  hugo <- RaggedExperiment::sparseAssay(mut, i = "Hugo_Symbol", withDimnames = TRUE)
  vc   <- RaggedExperiment::sparseAssay(mut, i = "Variant_Classification", withDimnames = TRUE)
  nonsilent <- !is.na(vc) & !(vc %in% SILENT_CLASSES)
  status <- matrix(0L, nrow = length(genes), ncol = ncol(hugo),
                   dimnames = list(genes, colnames(hugo)))
  for (g in genes) {
    hit <- (hugo == g) & nonsilent
    hit[is.na(hit)] <- FALSE
    status[g, colSums(hit) > 0L] <- 1L
  }
  status
}

#' Extract per-patient non-silent driver-mutation status (annotation only).
#'
#' Mutation is never a MOFA view: this is an external label for factor
#' interpretation (which factor tracks BAP1/PBRM1) on the n=417 subset.
#'
#' @param mae MultiAssayExperiment containing a "Mutation" experiment.
#' @param genes driver genes to report.
#' @return data.frame: character `sample_id` (harmonised patient IDs, unique)
#'   plus one integer 0/1 column per driver gene. New object; rownames are set
#'   to `sample_id` so Module 2's `fn_annotate_mutation` can index by rowname
#'   while Modules 3/4 use the `sample_id` column.
fn_extract_mutation_status <- function(mae, genes = DRIVER_GENES) {
  mut <- fn_experiment(mae, "Mutation")
  status_gs <- if (methods::is(mut, "RaggedExperiment")) {
    fn_mutation_status_ragged(mut, genes)
  } else {
    fn_mutation_status_dense(mut, genes)
  }
  ids <- fn_harmonise_ids(colnames(status_gs))
  keep <- !duplicated(ids)
  status_gs <- status_gs[, keep, drop = FALSE]
  ids <- ids[keep]

  out <- data.frame(sample_id = ids, stringsAsFactors = FALSE)
  sample_by_gene <- t(status_gs)  # sample x gene
  for (g in colnames(sample_by_gene)) {
    out[[g]] <- as.integer(sample_by_gene[, g])
  }
  rownames(out) <- ids
  out
}
