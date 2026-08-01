# Regenerate all Module-1 fixtures.
#   HEAVY_PULL=true  -> subsample the real curatedTCGAData MAE (local only)
#   otherwise        -> build a deterministic synthetic MAE (CI-safe default)
suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
})
source("R/constants.R")
source("R/functions_utils.R")   # fn_experiment resolves curatedTCGAData naming

fixture_dir <- "tests/fixtures"

# curatedTCGAData names every experiment "<DISEASE>_<stub>-<SNAPSHOT_DATE>".
# The synthetic MAE must use the same naming or the fixtures exercise a lookup
# path the real HEAVY_PULL never takes.
curated_name <- function(stub) paste0(CURATED_CANCER, "_", stub, "-", SNAPSHOT_DATE)

make_barcodes <- function(patients, type_code) {
  sprintf("TCGA-%s-%s-%sA-01R-1289-07", substr(patients, 1L, 2L),
          substr(patients, 3L, 6L), type_code)
}

fn_build_synthetic_mae <- function(seed = 1L) {
  set.seed(seed)
  patients <- sprintf("B0%04d", 1:24)                  # 24 synthetic patients
  tumour   <- make_barcodes(patients, PRIMARY_TUMOUR_CODE)
  normal   <- make_barcodes(patients[1:2], "11")       # 2 normal-tissue columns
  cols     <- c(tumour, normal)

  rna_genes <- c(paste0("GENE", 1:40), DRIVER_GENES)
  rna <- matrix(round(2^runif(length(rna_genes) * length(cols), 0, 14)),
                nrow = length(rna_genes),
                dimnames = list(rna_genes, cols))

  cpgs27  <- paste0("cg", sprintf("%08d", 1:60))
  cpgs450 <- paste0("cg", sprintf("%08d", 31:120))     # overlap cg31..cg60 with 27k
  m27_cols  <- c(tumour[1:20], normal)               # HM27 covers patients 1-20 + normals
  m450_cols <- tumour[5:24]                          # HM450 overlaps HM27 on patients 5-20
  m27  <- matrix(runif(length(cpgs27)  * length(m27_cols)),  nrow = length(cpgs27),
                 dimnames = list(cpgs27,  m27_cols))
  m450 <- matrix(runif(length(cpgs450) * length(m450_cols)), nrow = length(cpgs450),
                 dimnames = list(cpgs450, m450_cols))

  cnv_genes <- rna_genes
  cnv <- matrix(sample(-2:2, length(cnv_genes) * length(cols), replace = TRUE),
                nrow = length(cnv_genes), dimnames = list(cnv_genes, cols))

  mut <- matrix(sample(c(NA, "Missense_Mutation", "Silent"),
                       length(DRIVER_GENES) * length(tumour), replace = TRUE),
                nrow = length(DRIVER_GENES), dimnames = list(DRIVER_GENES, tumour))

  exps <- list(rna, m27, m450, cnv, mut)
  names(exps) <- vapply(
    c("RNASeq2GeneNorm", "Methylation_methyl27", "Methylation_methyl450",
      "GISTIC_ThresholdedByGene", "Mutation"),
    curated_name, character(1L)
  )
  MultiAssayExperiment(experiments = ExperimentList(exps))
}

heavy <- tolower(Sys.getenv("HEAVY_PULL", "false")) %in% c("true", "1", "yes")
mae <- if (heavy) {
  library(curatedTCGAData)
  full <- curatedTCGAData(CURATED_CANCER, CURATED_ASSAYS,
                          version = CURATED_VERSION, dry.run = FALSE)
  # Subset by SHARED patient IDs, not by position: MultiAssayExperiment has no
  # nrow/ncol method (min(30L, ncol(full)) silently degrades to seq_len(30) and
  # errors), and positional `j` takes columns 1..30 of each experiment
  # independently, leaving the assays with near-zero patient overlap.
  ids <- utils::head(rownames(MultiAssayExperiment::colData(full)), 30L)
  full[seq_len(200L), ids, ]
} else {
  fn_build_synthetic_mae(1L)
}

# The real experiments are SummarizedExperiment / RaggedExperiment objects,
# which as.matrix() cannot coerce directly.
as_mat <- function(x) if (is.matrix(x)) x else as.matrix(SummarizedExperiment::assay(x))

saveRDS(mae, file.path(fixture_dir, "kirc_mae_subset.rds"))
saveRDS(as_mat(fn_experiment(mae, "RNASeq2GeneNorm")),          file.path(fixture_dir, "rna_subset.rds"))
saveRDS(as_mat(fn_experiment(mae, "Methylation_methyl450")),    file.path(fixture_dir, "methyl_subset.rds"))
saveRDS(as_mat(fn_experiment(mae, "GISTIC_ThresholdedByGene")), file.path(fixture_dir, "cnv_subset.rds"))
mut <- fn_experiment(mae, "Mutation")
saveRDS(
  if (methods::is(mut, "RaggedExperiment")) {
    RaggedExperiment::sparseAssay(mut, i = "Variant_Classification", withDimnames = TRUE)
  } else {
    as_mat(mut)
  },
  file.path(fixture_dir, "mutation_subset.rds")
)
cat("fixtures written to", fixture_dir, "\n")
