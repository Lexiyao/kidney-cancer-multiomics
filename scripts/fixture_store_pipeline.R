# A SYNTHETIC `_targets` store, built so the dashboard's STORE-PRESENT render
# path can be exercised somewhere.
#
# READ THIS FIRST: every value below is a made-up stand-in. NOT ONE of them is a
# measurement, none may be quoted anywhere, and nothing rendered from this store
# may ever be published. Its only job is to make the branches that the storeless
# render can never reach — the plots, the tables, the value boxes — actually
# execute, in CI and locally.
#
# WHY IT EXISTS. Every .qmd page has two branches: store-absent (a stated gap)
# and store-present (the figure). Until this script, only the first was ever
# rendered outside a manual `pages.yml` dispatch, which is how a defect that
# aborts the whole render on the store-present path — `facet_wrap()` over an
# empty driver-gene selection, "Faceting variables must have at least one
# value." — sat undetected in `dashboard.qmd`.
#
# The shapes here are the PRODUCER shapes documented in `_targets.R`; they are
# what makes this a real test rather than a mock. Deliberately small: 40 cases,
# 8 genes, 15 factors.
#
# Usage (from the repo root):
#   Rscript -e 'targets::tar_make(script = "scripts/fixture_store_pipeline.R", store = "_targets")'
#   quarto render dashboard
#   rm -rf _targets _site

library(targets)

tar_option_set(packages = c("survival"))

FIXTURE_N        <- 40L
FIXTURE_FACTORS  <- paste0("Factor", 1:15)
FIXTURE_VIEWS    <- c("RNA", "Methylation", "CNV")
FIXTURE_IDS      <- sprintf("TCGA-FX-%04d", seq_len(FIXTURE_N))
FIXTURE_SEED     <- 20260812L

list(
  tar_target(fixture_ids, FIXTURE_IDS),

  # Module 1 shape: gene x sample. Deliberately includes only THREE of the six
  # DRIVER_GENES, so the partial-overlap branch of fn_driver_gene_frame() —
  # plot plus a note naming the dropped genes — is the one that renders here.
  tar_target(rna_mat, {
    set.seed(FIXTURE_SEED)
    genes <- c("VHL", "PBRM1", "BAP1", "AAA1", "AAA2", "AAA3", "AAA4", "AAA5")
    m <- matrix(stats::rnorm(length(genes) * FIXTURE_N, 8, 2),
                nrow = length(genes), dimnames = list(genes, fixture_ids))
    m
  }),

  tar_target(methyl_platform, {
    stats::setNames(rep(c("HM27", "HM450"), length.out = FIXTURE_N), fixture_ids)
  }),

  # Module 2 shapes.
  tar_target(mofa_factors, {
    set.seed(FIXTURE_SEED + 1L)
    matrix(stats::rnorm(FIXTURE_N * length(FIXTURE_FACTORS)), nrow = FIXTURE_N,
           dimnames = list(fixture_ids, FIXTURE_FACTORS))
  }),
  tar_target(mofa_varexp, {
    set.seed(FIXTURE_SEED + 2L)
    matrix(runif(length(FIXTURE_FACTORS) * length(FIXTURE_VIEWS), 0, 5),
           nrow = length(FIXTURE_FACTORS),
           dimnames = list(FIXTURE_FACTORS, FIXTURE_VIEWS))
  }),
  # One degenerate and one unscored row, because fn_factor_verdict() has
  # branches for both and a fixture that only exercises clean/confounded would
  # let either regress silently.
  tar_target(factor_platform, {
    q <- c(0.99, 1e-40, 2e-9, 0.22, 1e-30, 0, rep(0.4, 9))
    data.frame(
      factor  = FIXTURE_FACTORS,
      n       = FIXTURE_N,
      auc     = c(0.50, 0.89, 0.66, 0.53, 0.82, 0.74, rep(0.52, 8), NA_real_),
      p_value = q,
      q_value = q,
      degenerate = c(rep(FALSE, 5), TRUE, rep(FALSE, 9)),
      stringsAsFactors = FALSE
    )
  }),
  tar_target(mutation_factor_annot, {
    set.seed(FIXTURE_SEED + 3L)
    genes <- c("VHL", "PBRM1", "BAP1")
    out <- expand.grid(gene = genes, factor = FIXTURE_FACTORS,
                       stringsAsFactors = FALSE)
    out$p_value    <- runif(nrow(out))
    out$median_wt  <- stats::rnorm(nrow(out))
    out$median_mut <- stats::rnorm(nrow(out))
    out$q_value    <- stats::p.adjust(out$p_value, method = "BH")
    out[order(out$p_value), , drop = FALSE]
  }),
  tar_target(concordance, list(ari = 0.35)),

  # Module 5 consumption shapes.
  tar_target(subtypes_df, {
    data.frame(sample_id = fixture_ids,
               subtype = rep(c("S1", "S2", "S3", "S4"), length.out = FIXTURE_N),
               stringsAsFactors = FALSE)
  }),
  tar_target(survival_df, {
    set.seed(FIXTURE_SEED + 4L)
    data.frame(
      sample_id = fixture_ids,
      time      = stats::rexp(FIXTURE_N, 1 / 900),
      status    = stats::rbinom(FIXTURE_N, 1, 0.4),
      Factor1   = mofa_factors[, "Factor1"],
      Factor4   = mofa_factors[, "Factor4"],
      age_years = stats::rnorm(FIXTURE_N, 60, 9),
      stage_num = sample(1:4, FIXTURE_N, replace = TRUE),
      platform  = unname(methyl_platform),
      stringsAsFactors = FALSE
    )
  }),
  tar_target(km_subtype_df, merge(
    survival_df[, c("sample_id", "time", "status")], subtypes_df,
    by = "sample_id"
  )),

  # Module 4 shapes. A REAL coxph object, because `cox_fit` deserialises one on
  # the survival page and a list of numbers would not prove that path works.
  tar_target(cox_fit, {
    preds <- c("Factor1", "Factor4", "age_years", "stage_num", "platform")
    idx   <- seq_len(floor(nrow(survival_df) * 0.7))
    train <- survival_df[idx, , drop = FALSE]
    test  <- survival_df[-idx, , drop = FALSE]
    form  <- stats::as.formula(paste("survival::Surv(time, status) ~",
                                     paste(preds, collapse = " + ")))
    fit   <- survival::coxph(form, data = train, x = TRUE)
    list(model = fit, predictors = preds, train = train, test = test,
         n_events = sum(train$status == 1L),
         n_events_cohort = sum(survival_df$status == 1L),
         n_terms = length(preds), max_predictors = 12L)
  }),
  tar_target(survival_metrics, list(
    cindex   = list(cox = 0.70, penalised = 0.69, rsf = 0.71),
    optimism = list(cox = 0.02),
    calibration = data.frame(predicted = c(0.2, 0.5, 0.8),
                             observed  = c(0.25, 0.48, 0.75))
  )),
  tar_target(bap1_auroc, list(heldout_auroc = 0.80, cv_auroc = 0.78)),

  # Module 3 shape: sanity_table must carry a FAILING row, because the red-light
  # styling is the thing most worth proving renders.
  tar_target(sanity_table, data.frame(
    check  = c("mutation frequencies", "BAP1 worse OS", "m1-m4"),
    passed = c(TRUE, TRUE, FALSE),
    detail = c("synthetic stand-in, not a measurement",
               "synthetic stand-in; UNDERPOWERED",
               "synthetic stand-in; RED: platform ARI over the veto"),
    stringsAsFactors = FALSE
  )),

  # The live panel, offline. fn_query_gdc() would hit the network; a render
  # smoke test that needs the GDC API is a flaky test.
  tar_target(gdc_live_panel, list(
    project_id  = "TCGA-KIRC",
    total_cases = 999L,
    facets = list(
      demographic.vital_status = data.frame(category = c("Alive", "Dead"),
                                            n = c(600L, 399L),
                                            stringsAsFactors = FALSE),
      demographic.sex_at_birth = data.frame(category = c("male", "female"),
                                            n = c(650L, 349L),
                                            stringsAsFactors = FALSE),
      diagnoses.ajcc_pathologic_stage = data.frame(
        category = c("Stage I", "Stage III", "Stage IV", "Stage II"),
        n = c(500L, 200L, 180L, 119L), stringsAsFactors = FALSE
      )
    ),
    retrieved_at = as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
  )),

  # Module 6 (v1.1) shapes, for the same reason as every other target here:
  # without them the store-present branch of dashboard/singlecell.qmd is
  # rendered NOWHERE. The real `sc_mapping` is declared only behind
  # `run_singlecell`, which is off in CI and in the release, so the page's
  # populated path would otherwise be exercised for the first time on a manual
  # publish. These are stand-ins with no data behind them, exactly like the
  # rest of this file, and the gate is set to the FAILING arm because that is
  # the arm carrying the caveat a reader must not miss.
  # Shape must track fn_subtype_purity_test: fn_purity_gate_line REFUSES a
  # purity_check without the immune-arm evidence (immune_eta2, gate_arm), so a
  # stand-in lacking them aborts the store-present render this script exists to
  # exercise. Values are internally consistent stand-ins: both arms clear
  # p < 0.05 and eta2 >= 0.14, hence gate_arm = "purity+immune".
  tar_target(purity_check, list(
    purity_p = 1.2e-08, immune_p = 3.4e-05,
    purity_eta2 = 0.31, immune_eta2 = 0.22,
    is_purity_proxy = TRUE, gate_arm = "purity+immune", n = 40L
  )),
  # Shape must track map_bulk_signature: fn_sc_score_table REFUSES a mapping
  # without per-signature `coverage`, so the stand-in carries it too.
  tar_target(sc_mapping, list(
    scores_by_celltype = list(
      S1 = list(Tumor_epithelial = 0.42, T_cell = -0.11, Myeloid = 0.05),
      S2 = list(Tumor_epithelial = -0.07, T_cell = 0.38, Myeloid = 0.12)
    ),
    coverage = list(
      S1 = list(n_genes_requested = 20L, n_genes_used = 17L,
                low_coverage = FALSE),
      S2 = list(n_genes_requested = 20L, n_genes_used = 20L,
                low_coverage = FALSE)
    ),
    purity_confounded = TRUE,
    status = "exploratory",
    interpretation = paste(
      "SYNTHETIC STAND-IN, not a result. Bulk subtypes are a tumour-purity /",
      "immune-infiltration proxy (purity_check failed): these single-cell",
      "scores reflect cell-type composition, NOT a tumour-cell-intrinsic",
      "subtype program."
    )
  ))
)
