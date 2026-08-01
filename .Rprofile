source("renv/activate.R")

# --- MOFA2 basilisk/reticulate resolution (design spec section 5) ---
# In-container, RETICULATE_PYTHON + BASILISK_EXTERNAL_DIR are set by the
# Dockerfile; downstream MOFA2 runs call run_mofa(use_basilisk = FALSE) so no
# conda env is ever downloaded. Setting the option here keeps basilisk out of
# the per-user cache if it is ever invoked.
options(basilisk.useSystemDir = TRUE)

if (nzchar(Sys.getenv("BASILISK_EXTERNAL_DIR"))) {
  # honoured by basilisk for external env placement
  invisible(Sys.getenv("BASILISK_EXTERNAL_DIR"))
}
