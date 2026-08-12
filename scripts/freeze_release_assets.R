#!/usr/bin/env Rscript
# Package the frozen `_targets/` store into a tarball and publish it as a
# GitHub release asset, so CI / Pages / cron can restore the cache WITHOUT
# re-running the research core. The store itself is git-ignored (Task 0.1).
#
# WHAT THIS FREEZES IS THE 2016 SNAPSHOT, NOT A LIVE PULL. Every research
# target in the store descends from the curatedTCGAData 20160128 snapshot and
# is byte-stable across re-runs. The one target that is NOT frozen by this
# tarball in any meaningful sense is `gdc_live_panel`: it carries
# `tar_cue(mode = "always")`, so whatever value is captured here is a stale
# reading of a live API and is re-queried on the next run. Never quote a
# number restored from this asset as a current GDC portal figure.
store   <- "_targets"
tag     <- Sys.getenv("RELEASE_TAG", unset = "targets-store")
tarball <- "targets-store.tar.gz"

stopifnot(dir.exists(store))
# `compression=`, NOT the plan's `compress=`. utils::tar has BOTH `compression`
# and `compression_level`, so the abbreviation is ambiguous and R stops with
# "argument 3 matches multiple formal arguments" — verified on R 4.6.0. The
# plan's verbatim line never produces a tarball; this one does.
utils::tar(tarball, files = store, compression = "gzip")

# Create the release if absent (no-op if it already exists), then upload the
# asset with --clobber so the newest store always replaces the previous one.
system2("gh", c("release", "create", tag,
                "--title", shQuote("Frozen _targets store"),
                "--notes", shQuote("Cached pipeline outputs for CI/Pages render.")),
        stdout = FALSE, stderr = FALSE)
res <- system2("gh", c("release", "upload", tag, tarball, "--clobber"))
if (!identical(res, 0L)) stop("gh release upload failed")
cat("Uploaded", tarball, "to release", tag, "\n")
