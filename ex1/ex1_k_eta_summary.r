## Reads K_ETA_ALL.rds and reports the pre-registered contrasts from ex1_k_eta_screen.r.
##
## The quantity throughout is the WITHIN-DATASET contrast between the arms,
##     gap = stat - nonstat,
## because that -- not either arm's level -- is what the example claims. Configs are then compared
## on that gap, again within dataset. Pairing efficiency is printed first, before any comparison,
## so the reader can see whether the design could resolve the effects it reports.
a <- readRDS("K_ETA_ALL.rds")
stopifnot(nrow(a) > 0)
cfgs <- c("C0", "C1", "C2")
lab <- c(C0 = "K=4  eta=0.10  (committed)", C1 = "K=6  eta=0.10", C2 = "K=6  eta=0.05")

cat(sprintf("=== rows: %d   datasets: %d   configs: %s ===\n", nrow(a), length(unique(a$ds)),
  paste(sort(unique(a$config)), collapse = ", ")))
cat(sprintf("  rhat_M    max %.4f (nonstat) %.4f (stat)   fits over 1.01: %d of %d\n",
  max(a$nonstat_rhat_M), max(a$stat_rhat_M),
  sum(a$nonstat_rhat_M > 1.01) + sum(a$stat_rhat_M > 1.01), 2 * nrow(a)))
cat(sprintf("  divergences  total %d   fits with any: %d of %d\n",
  sum(a$nonstat_n_div) + sum(a$stat_n_div),
  sum(a$nonstat_n_div > 0) + sum(a$stat_n_div > 0), 2 * nrow(a)))

## ---- levels -------------------------------------------------------------------------------------
show <- function(v, digits = 3) {
  cat(sprintf("\n-- %s --\n%-28s %10s %10s %10s\n", v, "config", "nonstat", "stat", "gap"))
  for (cn in cfgs) {
    r <- a[a$config == cn, ]
    ns <- mean(r[[paste0("nonstat_", v)]]); st <- mean(r[[paste0("stat_", v)]])
    cat(sprintf("%-28s %10.*f %10.*f %10.*f\n", lab[cn], digits, ns, digits, st, digits, st - ns))
  }
}
## |mean_1| is the absolute error of the treated unit's period-1 effect: the true effect is exactly
## zero, so the posterior mean IS the error.
a$nonstat_abserr <- abs(a$nonstat_mean_1); a$stat_abserr <- abs(a$stat_mean_1)
for (v in c("noise_abs_tr", "abserr", "absz_1", "pred_perc", "pred_width", "eta_med", "pred_mad")) show(v)
cat("\n-- eta_prior_tail (how far the eta posterior sits into its own prior's upper tail) --\n")
for (cn in cfgs) {
  r <- a[a$config == cn, ]
  cat(sprintf("%-28s nonstat median %.3f     stat median %.3g\n", lab[cn],
    median(r$nonstat_eta_prior_tail), median(r$stat_eta_prior_tail)))
}

## ---- the paired contrasts -----------------------------------------------------------------------
gap_of <- function(v) {
  g <- tapply(a[[paste0("stat_", v)]] - a[[paste0("nonstat_", v)]], list(a$ds, a$config), mean)
  g[, cfgs, drop = FALSE]
}
pair <- function(g, from, to, name) {
  d <- g[, to] - g[, from]; d <- d[!is.na(d)]
  se <- sd(d) / sqrt(length(d))
  cat(sprintf("  %-22s %+.4f   se %.4f   95%% CI [%+.4f, %+.4f]   p %.3f   (n = %d)\n",
    name, mean(d), se, mean(d) - 1.96 * se, mean(d) + 1.96 * se,
    2 * pt(-abs(mean(d) / se), length(d) - 1), length(d)))
}
for (v in c("noise_abs_tr", "abserr", "absz_1")) {
  g <- gap_of(v)
  cat(sprintf("\n=== paired contrasts on the %s gap (stat - nonstat) ===\n", v))
  cat(sprintf("  between-dataset sd of the gap: %s\n",
    paste(sprintf("%s %.4f", cfgs, apply(g, 2, sd, na.rm = TRUE)), collapse = "   ")))
  cat(sprintf("  within-pair sd (C2 - C0):      %.4f   -> pairing buys %.1fx\n",
    sd(g[, "C2"] - g[, "C0"], na.rm = TRUE),
    sqrt(2) * mean(apply(g, 2, sd, na.rm = TRUE)) / sd(g[, "C2"] - g[, "C0"], na.rm = TRUE)))
  pair(g, "C0", "C2", "PRIMARY 1  C2 - C0")
  pair(g, "C0", "C1", "PRIMARY 2  C1 - C0 (K)")
  pair(g, "C1", "C2", "PRIMARY 3  C2 - C1 (eta)")
}

## ---- secondary: does either arm move on its own? -------------------------------------------------
cat("\n=== secondary: each arm's own level, paired across configs ===\n")
for (v in c("noise_abs_tr", "abserr")) for (arm in c("nonstat", "stat")) {
  m <- tapply(a[[paste0(arm, "_", v)]], list(a$ds, a$config), mean)[, cfgs, drop = FALSE]
  d <- m[, "C2"] - m[, "C0"]; d <- d[!is.na(d)]
  cat(sprintf("  %-8s %-14s C0 %.3f -> C2 %.3f   paired diff %+.4f se %.4f\n",
    arm, v, mean(m[, "C0"], na.rm = TRUE), mean(m[, "C2"], na.rm = TRUE),
    mean(d), sd(d) / sqrt(length(d))))
}
