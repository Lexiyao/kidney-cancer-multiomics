#!/usr/bin/env Rscript
# Local full pipeline run. Executes the entire `targets` DAG with heavy pulls
# ON (HM450 HDF5 download, MOFA2 training, etc.). Run ONCE locally; CI and the
# Pages site render from the frozen store this produces — never re-running here.
Sys.setenv(HEAVY_PULL = "true")
targets::tar_make()
cat("Full pipeline complete. Freeze the store with scripts/freeze_release_assets.R\n")
