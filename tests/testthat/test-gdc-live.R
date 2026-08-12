# Module 5 -- live GDC panel. EVERY TEST IN THIS FILE IS OFFLINE.
#
# The plan (Task 5.1) verifies the parser with an inline `Rscript -e` snippet on
# a hand-written bucket list and leaves the response-shaped half of the parse
# untested. That snippet cannot catch the two things that actually broke here:
# GDC content-negotiates to XML unless you ask for JSON, and it renamed
# `demographic.gender`. So the pure parse is unit-tested against RECORDED REAL
# responses committed under tests/fixtures/, captured 2026-08-10 from
# https://api.gdc.cancer.gov/cases. The network wrapper `fn_query_gdc` is NOT
# exercised here -- a test that needs the internet is a flaky test.
#
# The recorded counts below are a 2026-08-10 reading of the LIVE portal for the
# whole TCGA-KIRC project. They are NOT the frozen-snapshot cohort figures and
# must never be quoted as such; they are asserted only to prove the parser
# reproduces the bytes it was given.

gdc_fixture <- function(name) {
  jsonlite::fromJSON(
    testthat::test_path("..", "fixtures", name),
    simplifyVector = FALSE
  )
}

test_that("fn_parse_gdc_counts tidies buckets and sorts descending by n", {
  buckets <- list(list(key = "alive", doc_count = 372L),
                  list(key = "dead",  doc_count = 165L))

  res <- fn_parse_gdc_counts(buckets)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2L)
  expect_equal(res$category[1], "alive")
  expect_equal(res$n[1], 372L)
  expect_type(res$n, "integer")
})

test_that("fn_parse_gdc_counts re-sorts a response that is not already ordered", {
  buckets <- list(list(key = "small", doc_count = 1L),
                  list(key = "big",   doc_count = 99L))

  res <- fn_parse_gdc_counts(buckets)

  expect_equal(res$category, c("big", "small"))
})

test_that("fn_parse_gdc_counts returns an empty frame for no buckets", {
  res <- fn_parse_gdc_counts(list())

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0L)
  expect_named(res, c("category", "n"))
})

test_that("fn_extract_gdc_buckets pulls the facet out of a recorded response", {
  parsed <- gdc_fixture("gdc_cases_facet_response.json")

  buckets <- fn_extract_gdc_buckets(parsed, "demographic.vital_status")
  counts  <- fn_parse_gdc_counts(buckets)

  # Recorded 2026-08-10 from the LIVE GDC API -- not the 2016 snapshot.
  expect_equal(counts$category, c("alive", "dead"))
  expect_equal(counts$n, c(360L, 177L))
})

test_that("fn_extract_gdc_buckets fails loudly on a facet GDC does not serve", {
  parsed <- gdc_fixture("gdc_cases_facet_response.json")

  # `demographic.gender` is exactly the field the plan asked for and the field
  # the current GDC cases index dropped. Returning an empty table here would
  # put a silently blank panel on the site; the drift must surface.
  expect_error(
    fn_extract_gdc_buckets(parsed, "demographic.gender"),
    "demographic.gender"
  )
})

test_that("fn_parse_gdc_total reads the case total off a recorded response", {
  parsed <- gdc_fixture("gdc_cases_total_response.json")

  total <- fn_parse_gdc_total(parsed)

  expect_type(total, "integer")
  # Recorded 2026-08-10 from the LIVE GDC API. The frozen snapshot's colData
  # carries 536 cases; these are different quantities from different sources.
  expect_equal(total, 537L)
})

test_that("fn_parse_gdc_total rejects a response with no pagination block", {
  expect_error(fn_parse_gdc_total(list(data = list())), "pagination")
})

# --- Task 5.7: presentation of the live facets --------------------------------
# The facet KEYS are GDC field paths ("demographic.vital_status"). Printing them
# raw as headings on the page is unreadable, and munging them inline in the .qmd
# would put untested string logic on the one page whose whole job is to be read
# correctly.

test_that("fn_gdc_facet_title humanises a GDC field path", {
  # Arrange
  field <- "demographic.vital_status"

  # Act
  title <- fn_gdc_facet_title(field)

  # Assert
  expect_equal(title, "Vital status")
})

test_that("fn_gdc_facet_title keeps only the leaf of a nested path", {
  expect_equal(fn_gdc_facet_title("diagnoses.ajcc_pathologic_stage"),
               "Ajcc pathologic stage")
  expect_equal(fn_gdc_facet_title("demographic.sex_at_birth"),
               "Sex at birth")
})

test_that("fn_gdc_facet_title is vectorised and leaves a bare field alone", {
  out <- fn_gdc_facet_title(c("gender", "demographic.vital_status"))

  expect_equal(out, c("Gender", "Vital status"))
})

# --- Task 5.7: the live panel's retrieval stamp -------------------------------
#
# live-gdc.qmd promises "the timestamp is printed with every figure", and the
# per-facet charts carried no date at all. These charts render from a store
# restored from a release asset, so `retrieved_at` can be weeks older than the
# render, and a facet chart lifted onto a slide would carry no indication of
# when the portal was queried -- the FROZEN/LIVE confusion the page's own banner
# exists to prevent.

test_that("fn_gdc_stamp prints the retrieval time in UTC, not local time", {
  # Built in a non-UTC zone on purpose: the stamp must convert, not relabel.
  ts <- as.POSIXct("2026-01-02 03:04:05", tz = "America/Los_Angeles")

  expect_equal(fn_gdc_stamp(ts), "2026-01-02 11:04 UTC")
})

test_that("fn_gdc_stamp refuses a missing timestamp rather than printing none", {
  # A chart with a silently absent date is the failure this whole helper is for.
  expect_error(fn_gdc_stamp(NULL), "retrieval timestamp")
  expect_error(fn_gdc_stamp(as.POSIXct(NA)), "retrieval timestamp")
})

test_that("fn_gdc_facet_note carries both the field and the retrieval stamp", {
  note <- fn_gdc_facet_note("demographic.vital_status",
                            as.POSIXct("2026-01-02 03:04:05", tz = "UTC"))

  expect_match(note, "demographic.vital_status", fixed = TRUE)
  expect_match(note, "2026-01-02 03:04 UTC", fixed = TRUE)
  expect_match(note, fn_gdc_facet_title("demographic.vital_status"),
               fixed = TRUE)
})
