# Module 5 -- live GDC panel. Independent of the frozen research core.
#
# EVERYTHING IN THIS FILE READS THE LIVE PORTAL. Nothing here touches the
# frozen 20160128 curatedTCGAData snapshot the research targets are built from,
# and no value produced here is comparable to a cohort count: this is the whole
# TCGA-KIRC project as GDC holds it TODAY. Any page that renders these numbers
# must label them live, with `retrieved_at`, next to the snapshot figures they
# will otherwise be mistaken for.
#
# The file is split deliberately: `fn_parse_gdc_counts`, `fn_extract_gdc_buckets`
# and `fn_parse_gdc_total` are PURE and are unit-tested offline against recorded
# responses (tests/testthat/test-gdc-live.R); `fn_query_gdc` is the only part
# that needs the network and is verified as pipeline wiring, not in the suite.

#' Convert a GDC facet-aggregation bucket list into a tidy count table.
#' @param buckets list of lists, each with `key` and `doc_count`.
#' @return data.frame(category, n) sorted descending by n.
fn_parse_gdc_counts <- function(buckets) {
  stopifnot(is.list(buckets))
  if (length(buckets) == 0L) {
    return(data.frame(category = character(0), n = integer(0),
                      stringsAsFactors = FALSE))
  }
  category <- vapply(buckets, function(b) as.character(b[["key"]]), character(1))
  n        <- vapply(buckets, function(b) as.integer(b[["doc_count"]]), integer(1))
  out <- data.frame(category = category, n = n, stringsAsFactors = FALSE)
  out[order(-out$n), , drop = FALSE]
}

#' Pull one facet's bucket list out of a parsed GDC /cases response.
#' Pure: takes the already-parsed body, so it is testable without the network.
#' @param parsed parsed JSON body of a GDC /cases POST.
#' @param field facet field name, e.g. "demographic.vital_status".
#' @return list of buckets, ready for fn_parse_gdc_counts.
fn_extract_gdc_buckets <- function(parsed, field) {
  stopifnot(is.list(parsed), is.character(field), length(field) == 1L)
  agg <- parsed[["data"]][["aggregations"]][[field]]
  # A field GDC no longer serves comes back with NO aggregation and a
  # `warnings$facets` note. Degrading to an empty table would publish a blank
  # panel and call it a result; stop instead, so the weekly cron reports the
  # drift it exists to detect.
  if (is.null(agg)) {
    stop("GDC returned no aggregation for facet `", field,
         "`; the field name is stale or unsupported. Response warnings: ",
         paste(unlist(parsed[["warnings"]]), collapse = "; "))
  }
  buckets <- agg[["buckets"]]
  if (is.null(buckets)) list() else buckets
}

#' Read the total case count out of a parsed GDC /cases response.
#' Pure counterpart of the `size = 0` total query.
#' @param parsed parsed JSON body of a GDC /cases POST.
#' @return integer total number of cases matching the filter.
fn_parse_gdc_total <- function(parsed) {
  stopifnot(is.list(parsed))
  total <- parsed[["data"]][["pagination"]][["total"]]
  if (is.null(total)) {
    stop("GDC response carries no data$pagination$total; the endpoint or the ",
         "response format changed.")
  }
  as.integer(total)
}

#' Query current TCGA-KIRC case counts + clinical distribution from GDC.
#' @return list(project_id, total_cases, facets = named list of data.frames,
#'   retrieved_at). Never touches the frozen research targets.
fn_query_gdc <- function(project_id  = GDC_PROJECT_ID,
                         facet_fields = GDC_FACET_FIELDS,
                         base_url     = GDC_API_BASE) {
  endpoint <- paste0(base_url, "/cases")
  project_filter <- list(
    op = "in",
    content = list(field = "project.project_id", value = list(project_id))
  )
  # httr::accept_json() is LOAD-BEARING. Without it GDC content-negotiates to
  # `text/xml` (verified 2026-08-10), and forcing that body through the JSON
  # parser below errors. The plan's version omitted the header.
  post_json <- function(body) {
    resp <- httr::POST(url = endpoint, httr::accept_json(),
                       body = body, encode = "json")
    httr::stop_for_status(resp)
    httr::content(resp, as = "parsed", type = "application/json")
  }
  post_facet <- function(field) {
    parsed <- post_json(list(filters = project_filter, facets = field, size = 0L))
    fn_parse_gdc_counts(fn_extract_gdc_buckets(parsed, field))
  }
  total_parsed <- post_json(list(filters = project_filter, size = 0L))
  list(
    project_id   = project_id,
    total_cases  = fn_parse_gdc_total(total_parsed),
    facets       = stats::setNames(lapply(facet_fields, post_facet), facet_fields),
    retrieved_at = Sys.time()
  )
}

# --- Presentation of the live facets ------------------------------------------

#' Human-readable heading for a GDC facet field path.
#'
#' GDC facet keys are field paths ("demographic.vital_status"). The live panel
#' shows one chart per facet, and a raw field path as a heading is both ugly and
#' easy to misread. Kept here, next to the query that produces those keys, and
#' unit-tested, so live-gdc.qmd carries no untested string logic.
#'
#' @param field character vector of GDC field paths.
#' @return character vector of headings, same length.
fn_gdc_facet_title <- function(field) {
  stopifnot(is.character(field))
  leaf <- sub("^.*\\.", "", field)
  words <- gsub("_", " ", leaf, fixed = TRUE)
  paste0(toupper(substring(words, 1, 1)), substring(words, 2))
}

#' Format a live-panel retrieval timestamp for display, always in UTC.
#'
#' The live charts render from a `_targets` store restored from a release asset,
#' so `retrieved_at` can be weeks older than the render that shows it. A figure
#' with no date on it is indistinguishable from a current one once it is lifted
#' out of the page onto a slide, which is the FROZEN/LIVE confusion the panel's
#' provenance banner exists to prevent. An absent timestamp ERRORS rather than
#' rendering an undated chart.
#'
#' @param retrieved_at POSIXct(1), `gdc_live_panel$retrieved_at`.
#' @return character(1), "YYYY-MM-DD HH:MM UTC".
fn_gdc_stamp <- function(retrieved_at) {
  if (is.null(retrieved_at) || length(retrieved_at) != 1L ||
        is.na(retrieved_at)) {
    stop("the live panel carries no usable retrieval timestamp; a live figure ",
         "must not be rendered undated")
  }
  paste(format(retrieved_at, "%Y-%m-%d %H:%M", tz = "UTC"), "UTC")
}

#' Caption for one live facet chart: which GDC field, and when it was queried.
#'
#' @param field character(1) GDC field path.
#' @param retrieved_at POSIXct(1), `gdc_live_panel$retrieved_at`.
#' @return character(1).
fn_gdc_facet_note <- function(field, retrieved_at) {
  stopifnot(is.character(field), length(field) == 1L)
  sprintf("%s — live GDC field `%s`, retrieved %s",
          fn_gdc_facet_title(field), field, fn_gdc_stamp(retrieved_at))
}
