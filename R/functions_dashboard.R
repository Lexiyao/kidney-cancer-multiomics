# Module 5 -- adapters between the frozen research targets and the dashboard.
#
# Presentation only. Nothing here recomputes, re-fits or re-thresholds anything:
# it reformats what Modules 2-4 already measured so a page can render it. If a
# quantity is not in the store it does not appear here either — a dashboard that
# fills a gap with a plausible number is worse than one with a gap.

#' Flatten one `sanity_results` element into a one-line evidence string.
#'
#' The landing page renders `sanity_table` as check / passed / detail. A blank
#' detail column collapses a nuanced verdict into a bare tick or cross, and two
#' of the five checks cannot survive that: BAP1 survival PASSES its `hr > 1`
#' criterion while being underpowered, and m1-m4 FAILS for a specific,
#' defensible reason that `fn_check_methyl_strata` already put into words.
#'
#' @param el one element of `sanity_results` (a named list).
#' @return length-1 character; never empty.
fn_sanity_detail <- function(el) {
  stopifnot(is.list(el))
  if (!is.null(el[["detail"]])) return(as.character(el[["detail"]])[1])
  caveat <- if (isTRUE(el[["underpowered"]])) SANITY_UNDERPOWERED_NOTE else NULL
  parts <- c(as.character(el[["message"]]), caveat, fn_sanity_scalars(el))
  if (length(parts) == 0L) {
    return("passed/failed flag only; this check recorded no scalar evidence")
  }
  paste(parts, collapse = "; ")
}

#' Format the length-1 numeric fields of a sanity element as "name = value".
#' Only what the check actually measured is printed, so nothing can be quoted
#' here that the pipeline did not produce.
#' @param el one element of `sanity_results`.
#' @return character vector, possibly length 0.
fn_sanity_scalars <- function(el) {
  keep <- vapply(el, function(x) is.numeric(x) && length(x) == 1L, logical(1))
  scalars <- el[keep]
  if (length(scalars) == 0L) return(character(0))
  fmt <- function(nm) {
    paste0(nm, " = ", format(scalars[[nm]], digits = SANITY_DETAIL_DIGITS))
  }
  vapply(names(scalars), fmt, character(1), USE.NAMES = FALSE)
}

# --- Reading the frozen store from a rendering page ---------------------------

#' Read one frozen research target for a dashboard page.
#'
#' Same contract as `read_pipeline_target()` in tests/testthat/helper-fixtures.R,
#' and deliberately so: an ABSENT store (fresh clone, dev machine, CI before the
#' release asset is restored) is the one honest reason to have no number, and
#' the page degrades to a stated gap. A PRESENT store that cannot supply the
#' target is a broken restore, and it throws — a page that swallowed that would
#' show a gap where a real result exists upstream, which is a different lie from
#' the one this project usually worries about but a lie all the same.
#'
#' @param name character(1) target name.
#' @param store path to the `_targets` store.
#' @return the target's value, or NULL if the store holds no pipeline metadata.
fn_dashboard_read <- function(name, store = DASHBOARD_STORE_PATH) {
  stopifnot(is.character(name), length(name) == 1L)
  if (!file.exists(file.path(store, "meta", "meta"))) {
    return(NULL)
  }
  meta <- targets::tar_meta(store = store)
  if (!name %in% meta$name) {
    stop("the _targets store at ", store, " has pipeline metadata but records ",
         "no `", name, "` target: the restored release asset is older than ",
         "this page, or an upstream target failed before it.")
  }
  err <- meta$error[meta$name == name]
  if (!is.na(err[[1]])) {
    stop(name, " is recorded in the _targets store but ERRORED: ", err[[1]])
  }
  targets::tar_read_raw(name, store = store)
}

#' Read one CONDITIONALLY DECLARED target for the v1.1 single-cell page.
#'
#' Same guard as `fn_dashboard_read()` with one deliberate difference: a
#' populated store that does not record the target returns NULL instead of
#' throwing. For a Module 0-5 target that situation is a broken restore and must
#' throw. For Module 6 it is the NORMAL state — `sc_object`/`sc_mapping` are
#' declared only when `config$singlecell$run_singlecell` is TRUE (Task 6.8), and
#' CI and the frozen release build with it FALSE, so the release store legitimately
#' has no such target and the page must degrade rather than take the site down.
#' `purity_check` is read through here for the same reason: it is a v1.1 target
#' and a store frozen before Module 6 will not carry it.
#'
#' A target that IS recorded and ERRORED still throws. Absent and failed are not
#' the same thing, and rendering a failure as "not built yet" hides it.
#'
#' @param name character(1) target name.
#' @param store path to the `_targets` store.
#' @return the target's value, or NULL when the store cannot supply it.
fn_dashboard_read_optional <- function(name, store = DASHBOARD_STORE_PATH) {
  stopifnot(is.character(name), length(name) == 1L)
  if (!file.exists(file.path(store, "meta", "meta"))) {
    return(NULL)
  }
  meta <- targets::tar_meta(store = store)
  if (!name %in% meta$name) {
    return(NULL)
  }
  err <- meta$error[meta$name == name]
  if (!is.na(err[[1]])) {
    stop(name, " is recorded in the _targets store but ERRORED: ", err[[1]])
  }
  targets::tar_read_raw(name, store = store)
}

# --- Module 6: the purity confound gate, as a reader sees it ------------------

#' One-line verdict for the subtype-vs-purity confound gate (spec section 10).
#'
#' The gate decides whether bulk subtypes may be mapped onto single cells at
#' face value, so its two verdicts must be distinguishable at a glance: FAILED
#' is rendered in the same red used for a failing positive control, passed in
#' the same green. The evidence for BOTH tested arms (TumorPurity and
#' ImmuneScore: p, eta-squared, n), which arm fired, and the two thresholds
#' that produced the verdict all travel with it, so the sentence can be quoted
#' without the reader having to trust the word alone.
#'
#' The printed decision rule is the one the code applies, character for
#' character: `fn_subtype_purity_test` compares p with a STRICT `<` and
#' eta-squared with an inclusive `>=`. This line previously printed "p <=",
#' which is a different rule from the one that produced the verdict beside it.
#'
#' An NA or non-logical verdict STOPS rather than defaulting. `isTRUE(NA)` is
#' FALSE, so a defaulting formatter would print "passed" for a gate that never
#' returned one — the same unsafe direction `_validate_gate_verdict` refuses on
#' the Python side.
#'
#' @param purity_check the `purity_check` target (fn_subtype_purity_test).
#' @return character(1) of inline HTML.
fn_purity_gate_line <- function(purity_check) {
  stopifnot(is.list(purity_check))
  required <- c("purity_p", "purity_eta2", "immune_p", "immune_eta2",
                "n", "is_purity_proxy", "gate_arm")
  missing_fields <- setdiff(required, names(purity_check))
  if (length(missing_fields) > 0L) {
    stop("purity_check lacks: ", paste(missing_fields, collapse = ", "),
         ". The gate verdict may not be shown without the evidence that ",
         "produced it.")
  }
  proxy <- purity_check$is_purity_proxy
  if (!is.logical(proxy) || length(proxy) != 1L || is.na(proxy)) {
    stop("purity_check$is_purity_proxy must be a single TRUE or FALSE; a ",
         "missing verdict must not be rendered as a passing one.")
  }
  sprintf(
    paste0("Kruskal-Wallis, n = %d. TumorPurity ~ subtype: p = %.3g, ",
           "epsilon-corrected eta&sup2; = %.3f. ",
           "ImmuneScore ~ subtype: p = %.3g, eta&sup2; = %.3f. ",
           "(The gate fires when EITHER arm has p &lt; %g AND eta&sup2; &ge; ",
           "%g; arm fired: %s.) %s"),
    purity_check$n,
    purity_check$purity_p, purity_check$purity_eta2,
    purity_check$immune_p, purity_check$immune_eta2,
    PURITY_PROXY_ALPHA, PURITY_PROXY_ETA2, purity_check$gate_arm,
    fn_purity_gate_span(proxy)
  )
}

#' The verdict half of `fn_purity_gate_line()`, coloured like every other
#' pass/fail verdict on the site.
#' @param proxy length-1 logical, already validated.
#' @return character(1) of inline HTML.
fn_purity_gate_span <- function(proxy) {
  if (proxy) {
    return(paste0("<span class='verdict-red'>GATE FAILED &mdash; the bulk ",
                  "subtypes ARE a tumour-purity / immune-infiltration ",
                  "proxy.</span>"))
  }
  # Deliberately narrow. A gate that did not fire is not evidence that the
  # subtypes are tumour-cell-intrinsic; it is evidence that neither tested
  # confound reached the stated effect size. At n = 524 with S1 at 20 samples,
  # a real but moderate gradient can leave the eta-squared arm silent, so the
  # sentence must not be readable as a positive licence.
  paste0("<span class='verdict-green'>gate did not fire &mdash; neither ",
         "tumour purity nor immune infiltration separates the subtypes at ",
         "this effect size. That is not evidence the subtypes are ",
         "tumour-cell-intrinsic.</span>")
}

#' Reshape `sc_mapping$scores_by_celltype` into a cell-type x subtype table.
#'
#' The mapping arrives as subtype -> cell type -> mean score (the shape
#' `map_bulk_signature` returns). Cell types are the rows because that is the
#' comparison the gate is about: if the subtype signatures separate on cell type
#' rather than within tumour cells, the reader should be able to see it in the
#' row pattern.
#'
#' Every column header carries the signature's GENE COVERAGE
#' ("S1 (18/20 genes)"): the signatures are derived on bulk TCGA expression and
#' intersected with the post-QC scRNA gene space, so partial dropout is the
#' expected case, and a score computed from 2 of 20 genes must not print with
#' the authority of one computed from all 20. A mapping without per-signature
#' coverage is REFUSED for the same reason the empty mapping is: a number
#' without the count behind it reads as more than was measured. An
#' unmeasurable signature arrives as NaN (map_bulk_signature never
#' zero-fills) and is rendered as NaN, not as 0.000.
#'
#' @param sc_mapping the `sc_mapping` target.
#' @return data.frame(cell_type, <one numeric column per subtype, coverage in
#'   the column name>).
fn_sc_score_table <- function(sc_mapping) {
  stopifnot(is.list(sc_mapping))
  scores <- sc_mapping[["scores_by_celltype"]]
  if (!is.list(scores) || length(scores) == 0L) {
    stop("sc_mapping carries no `scores_by_celltype`; there is nothing to ",
         "tabulate and an empty table would read as 'no signal'.")
  }
  coverage <- sc_mapping[["coverage"]]
  if (!is.list(coverage) || !all(names(scores) %in% names(coverage))) {
    stop("sc_mapping carries no per-signature `coverage`; a score may not be ",
         "shown without the gene count that produced it (a value computed ",
         "from 2 of 20 genes must not print with full-panel authority).")
  }
  cell_types <- unique(unlist(lapply(scores, names), use.names = FALSE))
  columns <- lapply(scores, function(s) as.numeric(unlist(s)[cell_types]))
  names(columns) <- vapply(names(scores), function(nm) {
    cov <- coverage[[nm]]
    sprintf("%s (%s/%s genes)", nm,
            format(cov$n_genes_used), format(cov$n_genes_requested))
  }, character(1))
  do.call(data.frame, c(list(cell_type = cell_types), columns,
                        list(check.names = FALSE, stringsAsFactors = FALSE)))
}

#' Why the single-cell panels are empty, stated in their own terms.
#'
#' Deliberately NOT `fn_dashboard_pending()`: that sentence tells the reader to
#' restore the `targets-store` release asset, and here the release asset is
#' complete and simply contains no single-cell targets, because it was built
#' with `run_singlecell: false`. A true gap explained by a false cause sends the
#' reader to fix something that is not broken.
#'
#' @param what character(1) description of the missing quantity.
#' @return character(1).
fn_singlecell_gap <- function(what) {
  stopifnot(is.character(what), length(what) == 1L)
  paste0(
    "NOT IN THIS BUILD &mdash; ", what, " is not shown because the Module 6 ",
    "targets were not built here. They are declared only when ",
    "<code>run_singlecell</code> is TRUE in <code>config/params.yml</code>, ",
    "and CI and the frozen release both build with it FALSE so the GSE159115 ",
    "download never runs on a runner. Nothing is displayed in its place."
  )
}

# --- Per-factor platform verdicts ---------------------------------------------

#' Annotate MOFA factors with their measured association to the assay split.
#'
#' `methyl_mat` is `cbind(HM27, HM450)` with no batch correction, and several
#' factors are substantially that split rather than biology (Factor2 AUC 0.888).
#' Nothing on the site may show a factor without this number beside it. The
#' verdict uses the SAME `q <= SANITY_MAX_P` rule Phase 4 used to select
#' `PLATFORM_CLEAN_MOFA_FACTORS`, so the page cannot call a factor clean that
#' the survival model treated as dirty.
#'
#' @param factor_names character vector of factors, in display order.
#' @param factor_platform the `factor_platform` target
#'   (`data.frame(factor, n, auc, p_value, q_value, degenerate)`).
#' @return data.frame(factor, auc, q_value, verdict, label), one row per
#'   requested factor, in the requested order.
fn_factor_platform_labels <- function(factor_names, factor_platform) {
  stopifnot(is.character(factor_names), is.data.frame(factor_platform))
  required <- c("factor", "auc", "q_value", "degenerate")
  missing_cols <- setdiff(required, names(factor_platform))
  if (length(missing_cols) > 0L) {
    stop("factor_platform lacks column(s): ",
         paste(missing_cols, collapse = ", "))
  }
  idx <- match(factor_names, factor_platform$factor)
  auc <- factor_platform$auc[idx]
  qv  <- factor_platform$q_value[idx]
  verdict <- fn_factor_verdict(auc, qv, factor_platform$degenerate[idx])
  data.frame(
    factor  = factor_names,
    auc     = auc,
    q_value = qv,
    verdict = verdict,
    label   = fn_factor_label(factor_names, auc, verdict),
    stringsAsFactors = FALSE
  )
}

#' Verdict for one factor's platform association. Vectorised.
#' An unscored factor is an UNKNOWN and is never defaulted to clean.
#' @return character vector of FACTOR_VERDICT_* values.
fn_factor_verdict <- function(auc, q_value, degenerate) {
  ifelse(
    is.na(auc) | is.na(q_value), FACTOR_VERDICT_UNSCORED,
    ifelse(
      fn_is_true(degenerate), FACTOR_VERDICT_DEGENERATE,
      ifelse(q_value <= SANITY_MAX_P, FACTOR_VERDICT_CONFOUNDED,
             FACTOR_VERDICT_CLEAN)
    )
  )
}

#' NA-safe elementwise isTRUE, so an unscored row cannot be read as degenerate.
#' @return logical vector, no NAs.
fn_is_true <- function(x) !is.na(x) & x

#' Axis/table label carrying the AUC beside the factor name.
#' @return character vector, same length as `factor_names`.
fn_factor_label <- function(factor_names, auc, verdict) {
  ifelse(
    is.na(auc),
    paste0(factor_names, "  [", verdict, "]"),
    sprintf("%s  [AUC %.3f, %s]", factor_names, auc, verdict)
  )
}

# --- Verdict cells: a failing result must LOOK failing ------------------------

#' HTML span for a positive control's PASS/FAIL verdict.
#'
#' dashboard/styles.css states the rule this implements: styling a failing
#' anchor like a passing one is how a negative result gets laundered into a
#' positive one. It lives here, and not in a page, because index.qmd and
#' dashboard.qmd render the SAME `sanity_table` and had already drifted apart —
#' index emitted the red span, dashboard printed "FAIL (RED)" in body text.
#' Every caller must pass `escape = FALSE` to `knitr::kable()`.
#'
#' @param passed logical vector from `sanity_table$passed`.
#' @return character vector of HTML spans, same length.
fn_pass_span <- function(passed) {
  stopifnot(is.logical(passed))
  ifelse(passed,
         "<span class='verdict-green'>PASS</span>",
         "<span class='verdict-red'>FAIL (RED)</span>")
}

#' HTML span for one MOFA factor's platform verdict.
#'
#' A degenerate AUC is styled RED, not neutral: "AUC UNUSABLE" is a refusal to
#' certify the factor, and a reader scanning the colour column must not read it
#' as a clean bill of health. An unscored factor is amber — unknown, which is
#' also not clean.
#'
#' @param verdict character vector of FACTOR_VERDICT_* values.
#' @return character vector of HTML spans, same length.
fn_verdict_span <- function(verdict) {
  stopifnot(is.character(verdict))
  classes <- c(
    "verdict-green",   # FACTOR_VERDICT_CLEAN
    "verdict-red",     # FACTOR_VERDICT_CONFOUNDED
    "verdict-red",     # FACTOR_VERDICT_DEGENERATE
    "verdict-pending"  # FACTOR_VERDICT_UNSCORED
  )
  names(classes) <- c(FACTOR_VERDICT_CLEAN, FACTOR_VERDICT_CONFOUNDED,
                      FACTOR_VERDICT_DEGENERATE, FACTOR_VERDICT_UNSCORED)
  unknown <- setdiff(verdict, names(classes))
  if (length(unknown) > 0L) {
    stop("no styling rule for verdict(s): ", paste(unknown, collapse = ", "),
         ". A verdict with no colour renders as ordinary body text, which is ",
         "exactly how a failing result stops looking failing.")
  }
  sprintf("<span class='%s'>%s</span>", classes[verdict], verdict)
}

#' Standing sentence for a quantity the pipeline has not produced in THIS render.
#'
#' The gap is shown as a gap. Filling it with a plausible value is the failure
#' mode this whole file is written against.
#' @param what character(1) description of the missing quantity.
#' @return character(1).
fn_dashboard_pending <- function(what) {
  stopifnot(is.character(what), length(what) == 1L)
  paste0(
    "PENDING — ", what, " is not shown because the frozen `_targets` store is ",
    "not present in this render. It is a measured pipeline output, not a ",
    "number this page may estimate, so nothing is displayed in its place. ",
    "Restore the `targets-store` release asset and re-render."
  )
}

# --- Survival page: the number, and the two things that qualify it ------------

#' Held-out C-index per model arm, with the Cox optimism REPORTED beside it.
#'
#' `survival_metrics$optimism$cox` is apparent(train) - validated(held-out). It
#' is displayed, never subtracted away: subtypes and factors were derived from
#' the same omics, so the apparent figure is optimistic by construction and the
#' gap is the finding. Only the Cox arm has an optimism estimate; the other two
#' rows stay NA rather than borrowing it.
#'
#' @param survival_metrics the `survival_metrics` target.
#' @return data.frame(arm, heldout, optimism, apparent); NA where not measured.
fn_survival_cindex_table <- function(survival_metrics) {
  stopifnot(is.list(survival_metrics))
  ci <- survival_metrics[["cindex"]]
  if (!is.list(ci)) {
    stop("survival_metrics carries no `cindex` list; there is no held-out ",
         "C-index to report and none may be substituted")
  }
  pick <- function(x, nm) {
    v <- x[[nm]]
    if (is.null(v) || length(v) != 1L) NA_real_ else as.numeric(v)
  }
  heldout <- c(pick(ci, "cox"), pick(ci, "penalised"), pick(ci, "rsf"))
  opt_cox <- pick(survival_metrics[["optimism"]], "cox")
  optimism <- c(opt_cox, NA_real_, NA_real_)
  data.frame(
    arm = c("Cox (primary)", "Penalised Cox (glmnet)",
            "Random survival forest"),
    heldout  = heldout,
    optimism = optimism,
    apparent = heldout + optimism,
    stringsAsFactors = FALSE
  )
}

#' Model-size panel: split, event counts and EPV budget, read off `cox_fit`.
#'
#' Every figure here was recorded by fn_fit_cox at fit time. The page reads
#' them; it does not recompute them, because a second derivation of a number
#' the pipeline owns is how the two silently drift apart.
#'
#' @param cox_fit the `cox_fit` target.
#' @return data.frame(quantity, value), both character-friendly for kable.
fn_cox_fit_summary <- function(cox_fit) {
  stopifnot(is.list(cox_fit))
  required <- c("train", "test", "predictors", "n_events", "n_events_cohort",
                "n_terms", "max_predictors")
  missing_fields <- setdiff(required, names(cox_fit))
  if (length(missing_fields) > 0L) {
    stop("cox_fit lacks: ", paste(missing_fields, collapse = ", "))
  }
  data.frame(
    quantity = c("Training rows", "Held-out rows", "OS events in the cohort",
                 "OS events in TRAINING (the EPV base)",
                 "Model terms used vs EPV-10 cap", "Predictors fitted"),
    value = c(
      format(nrow(cox_fit$train)), format(nrow(cox_fit$test)),
      format(cox_fit$n_events_cohort), format(cox_fit$n_events),
      paste(cox_fit$n_terms, "of", cox_fit$max_predictors),
      paste(cox_fit$predictors, collapse = ", ")
    ),
    stringsAsFactors = FALSE
  )
}

# --- Main dashboard: value boxes and the KM frame -----------------------------

#' Build one Quarto dashboard value box, or a stated gap.
#'
#' The plan formats each box as `sprintf("%.3f", metrics$cindex$cox)`. On a
#' target that is absent from this render that expression returns
#' `character(0)`, and Quarto renders an EMPTY box — a blank standing where a
#' reader expects the headline metric, with nothing to say why. Three of the
#' four boxes are container-produced numbers that simply are not there without
#' the restored store, so the absent case is the normal case, not an edge case.
#'
#' @param value length-1 numeric, or NULL/NA when the target is unavailable.
#' @param icon bootstrap icon name.
#' @param colour value-box colour used when the value IS available.
#' @param fmt sprintf format for the available case.
#' @return list(value, icon, color) in Quarto's value-box shape.
fn_dashboard_valuebox <- function(value, icon, colour, fmt = "%.3f") {
  stopifnot(is.character(icon), is.character(colour))
  if (!fn_metric_usable(value)) {
    return(list(value = DASHBOARD_PENDING_VALUE, icon = icon,
                color = DASHBOARD_PENDING_COLOUR))
  }
  list(value = sprintf(fmt, value), icon = icon, color = colour)
}

#' Is this a single measured value the page may print?
#'
#' NULL (target absent from this render), NA, `character(0)`/`numeric(0)` (the
#' shape `x$y$z` returns off a NULL list) and length > 1 all mean "the pipeline
#' did not hand me one number", and they are all treated identically.
#' @param value the candidate.
#' @return TRUE/FALSE, length 1.
fn_metric_usable <- function(value) {
  !is.null(value) && length(value) == 1L && !is.na(value)
}

#' Format one headline metric for a table cell, or state the gap.
#'
#' The landing page's summary table is built with `data.frame()`, and a cell
#' that evaluates to `character(0)` does not render blank there -- it aborts the
#' whole render with "arguments imply differing number of rows". This always
#' returns exactly one string, and when the value is missing that string SAYS
#' the value is missing rather than substituting a plausible one.
#'
#' @param value length-1 numeric, or NULL/NA when the target is unavailable.
#' @param fmt sprintf format applied when the value IS available.
#' @return character(1).
fn_metric_text <- function(value, fmt = "%.3f") {
  stopifnot(is.character(fmt), length(fmt) == 1L)
  if (!fn_metric_usable(value)) return(DASHBOARD_PENDING_VALUE)
  sprintf(fmt, value)
}

# --- Main dashboard: the driver-gene panel ------------------------------------

#' Long frame of driver-gene expression by subtype, plus what is NOT in it.
#'
#' `rna_mat` is `fn_top_variable(rna_full, N_TOP_GENES)` — the 5000 most-variable
#' genes of roughly 20500 — so there is NO guarantee that VHL/PBRM1/SETD2/BAP1/
#' MTOR/KDM5C survive the filter. Two failures follow from ignoring that, and
#' both are why the selection lives here rather than inline in a chunk:
#' `facet_wrap()` on a zero-row frame aborts the whole site render, and a
#' partial panel that silently shows 2 of 6 genes under-reports without saying
#' so. The dropped genes are returned so the page can name them.
#'
#' @param rna_mat gene x sample expression matrix (the `rna_mat` target).
#' @param subtypes_df data.frame(sample_id, subtype) (the `subtypes_df` target).
#' @param genes character vector of genes to show; defaults to `DRIVER_GENES`.
#' @return list(present, dropped, frame); `frame` is data.frame(gene, expr,
#'   subtype) and has zero rows when there is nothing to plot.
fn_driver_gene_frame <- function(rna_mat, subtypes_df, genes = DRIVER_GENES) {
  stopifnot(!is.null(rownames(rna_mat)), is.data.frame(subtypes_df),
            is.character(genes))
  present <- genes[genes %in% rownames(rna_mat)]
  common  <- intersect(colnames(rna_mat), subtypes_df$sample_id)
  empty <- data.frame(gene = character(0), expr = numeric(0),
                      subtype = character(0), stringsAsFactors = FALSE)
  frame <- if (length(present) == 0L || length(common) == 0L) {
    empty
  } else {
    data.frame(
      gene = rep(present, each = length(common)),
      expr = as.numeric(t(rna_mat[present, common, drop = FALSE])),
      subtype = rep(
        as.character(subtypes_df$subtype[match(common, subtypes_df$sample_id)]),
        times = length(present)
      ),
      stringsAsFactors = FALSE
    )
  }
  list(present = present, dropped = setdiff(genes, present), frame = frame)
}

#' Why the driver-gene panel is empty, stated in the panel's own terms.
#'
#' Deliberately NOT routed through `fn_dashboard_pending()`: that sentence blames
#' an unrestored release asset, and this gap occurs with the store PRESENT. A
#' true gap explained by a false cause is its own kind of dishonesty.
#'
#' @param gd the value of `fn_driver_gene_frame()`.
#' @return character(1).
fn_driver_gene_gap <- function(gd) {
  stopifnot(is.list(gd))
  cause <- if (length(gd$present) == 0L) {
    paste0("none of the driver genes (", paste(gd$dropped, collapse = ", "),
           ") survived the top-", N_TOP_GENES, " variance filter that produces ",
           "`rna_mat`")
  } else {
    "the expression matrix and the subtype table share no sample"
  }
  paste0("NOT SHOWN — driver-gene expression by subtype is not plotted because ",
         cause, ". The frozen store is present in this render; the gap is in ",
         "the feature set, not in the store, and no substitute genes are ",
         "plotted in their place.")
}

#' One-line under-reporting note for a partially populated driver-gene panel.
#' @param dropped character vector of genes absent from `rna_mat`.
#' @return character(1), or NULL when nothing was dropped.
fn_driver_gene_note <- function(dropped) {
  stopifnot(is.character(dropped))
  if (length(dropped) == 0L) return(NULL)
  paste0("Showing ", length(DRIVER_GENES) - length(dropped), " of ",
         length(DRIVER_GENES), " driver genes. ", paste(dropped, collapse = ", "),
         " did not survive the top-", N_TOP_GENES, " variance filter applied to ",
         "`rna_mat`, so ", if (length(dropped) == 1L) "it is" else "they are",
         " absent from this panel — not absent from the tumours.")
}

#' Expand a `survfit` object into a tidy step-plot frame for Plotly.
#'
#' Shared by survival.qmd and dashboard.qmd so the two pages cannot draw the
#' same curves from two slightly different expansions.
#'
#' @param fit a `survival::survfit` object stratified by subtype.
#' @return data.frame(time, surv, strata); `strata` carries the bare level.
fn_km_curve_df <- function(fit) {
  stopifnot(inherits(fit, "survfit"))
  data.frame(
    time   = fit$time,
    surv   = fit$surv,
    strata = rep(sub("^[^=]*=", "", names(fit$strata)), fit$strata),
    stringsAsFactors = FALSE
  )
}
