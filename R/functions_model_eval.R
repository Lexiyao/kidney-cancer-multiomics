# Module 4 evaluation — reuse of model-evaluation-from-scratch.
# C-index + calibration for survival models. No external metric packages.

fn_cindex <- function(time, status, risk) {
  # Harrell's concordance from scratch. A pair (i, j) is permissible when
  # subject i has an event and j OUTLIVES i — a strictly later recorded time,
  # or the SAME recorded time with j censored (j was still event-free when i
  # died, so j outlived i). That second clause is the survival::concordance
  # convention and it is NOT cosmetic: dropping it (requiring a strictly later
  # time) agrees with the reference only when times are continuous. MEASURED
  # (R 4.6.0 / survival 3.8.6): identical to 1.1e-16 on 10 continuous cohorts,
  # but off by up to 4.7e-04 on 25 integer-day cohorts and 1.1e-02 on 25
  # coarse-grid cohorts. TCGA follow-up is recorded in whole DAYS, so the real
  # fit lives in the tied regime. A pair is concordant when i (the one who
  # died first) carries the higher risk; ties in risk score 0.5.
  stopifnot(length(time) == length(status), length(time) == length(risk))
  n <- length(time)
  concordant <- 0
  permissible <- 0
  tied_risk <- 0
  for (i in seq_len(n)) {
    if (status[i] != 1L) next
    for (j in seq_len(n)) {
      if (i == j) next
      later <- time[j] > time[i] || (time[j] == time[i] && status[j] != 1L)
      if (!later) next
      permissible <- permissible + 1
      if (risk[i] > risk[j]) {
        concordant <- concordant + 1
      } else if (risk[i] == risk[j]) {
        tied_risk <- tied_risk + 1
      }
    }
  }
  if (permissible == 0) return(NA_real_)
  (concordant + 0.5 * tied_risk) / permissible
}

fn_km_at <- function(km, horizon) {
  # Observed survival probability at `horizon` from a fitted survfit;
  # extend = TRUE carries the last step forward past the final event time.
  s <- summary(km, times = horizon, extend = TRUE)$surv
  if (length(s) == 0) return(NA_real_)
  s[length(s)]
}

fn_calibration <- function(time, status, pred_surv, horizon,
                           n_bins = CALIBRATION_BINS) {
  # Grouped calibration: bin subjects by predicted survival probability,
  # compare mean predicted vs observed Kaplan-Meier survival at `horizon`.
  stopifnot(length(time) == length(status), length(time) == length(pred_surv))
  brks <- stats::quantile(pred_surv, probs = seq(0, 1, length.out = n_bins + 1),
                          na.rm = TRUE)
  bins <- cut(pred_surv, breaks = brks, include.lowest = TRUE, labels = FALSE)
  rows <- lapply(sort(unique(bins)), function(b) {
    idx <- which(bins == b)
    km <- survival::survfit(survival::Surv(time[idx], status[idx]) ~ 1)
    data.frame(bin = b,
               n = length(idx),
               predicted = mean(pred_surv[idx]),
               observed = fn_km_at(km, horizon))
  })
  do.call(rbind, rows)
}
