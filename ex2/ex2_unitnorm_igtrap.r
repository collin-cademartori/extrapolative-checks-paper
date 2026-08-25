## Decisive (NOT shipped): the unit-norm ints rhat spikes look like STOCHASTIC trapping in shallow minor
## modes (chains ~10 lp below best; whole-vector rhat inflation; same seed/config irreproducible across
## runs). Test by REPEATING each fit across many Stan seeds and reporting the DISTRIBUTION of rhat/offmode.
## Arms: clean model 4ch, INV-GAMMA model 4ch (does zero-avoidance kill the minor modes?), clean 8ch (does
## more chains help?). 2 datasets that trapped (i1, i4). If clean_4ch has a nonzero spike rate that ig_4ch
## or clean_8ch drives to ~0, we have the cause AND the fix.

library(foreach); library(doParallel); library(posterior)
source("../sample_model.r"); source("../pathfinder_init.r")

ruv <- function(d) { v <- rnorm(d); v / sqrt(sum(v * v)) }
anchor_order <- function(y, K) {
  N <- ncol(y); yc <- scale(y, center = TRUE, scale = FALSE); sel <- 1L
  remaining <- setdiff(seq_len(N), sel)
  while (length(sel) < K && length(remaining) > 0) {
    Q <- qr.Q(qr(yc[, sel, drop = FALSE]))
    resid <- yc[, remaining, drop = FALSE] - Q %*% crossprod(Q, yc[, remaining, drop = FALSE])
    pick <- remaining[which.max(colSums(resid^2))]; sel <- c(sel, pick); remaining <- setdiff(remaining, pick)
  }
  c(sel, remaining)
}
sim_model_intercepts <- function(N_unc = 2, N_comp_true = 2, N_comp_spur = 2, T_times = 30,
                                 T_treated = 5, K_unc = 1, sim = 0.9) {
  N_units <- 1 + N_comp_true + N_comp_spur + N_unc; K_gen <- 2 + K_unc
  f_treat <- 6 + arima.sim(model = list(ar = 0.9), n = T_times); f_treat_sd <- 1.9
  f_alt <- (f_treat - 6) + c(rep(0, T_times - T_treated), rep(-f_treat_sd, T_treated))
  f_unc <- matrix(nrow = K_unc, ncol = T_times); cor_unc <- Inf
  while (cor_unc > 0.01) {
    for (k in 1:K_unc) f_unc[k, ] <- rnorm(1, 0, 2) + arima.sim(model = list(ar = 0.9), n = T_times)
    cor_unc <- max(abs(cor(t(f_unc), t(t(f_treat)))))
  }
  facs <- rbind(f_treat, f_alt, f_unc); loads <- matrix(nrow = N_units, ncol = K_gen)
  loads[1, ] <- c(1, rep(0, K_gen - 1))
  for (n in seq_len(N_comp_true)) loads[1 + n, ] <- c(sqrt(sim), 0, sqrt(1 - sim) * ruv(K_gen - 2))
  for (n in seq_len(N_comp_spur)) loads[1 + N_comp_true + n, ] <- c(0, sqrt(sim), sqrt(1 - sim) * ruv(K_gen - 2))
  for (n in seq_len(N_unc)) loads[1 + N_comp_true + N_comp_spur + n, ] <- c(0, 0, ruv(K_gen - 2))
  lat <- loads %*% facs
  Y <- t(lat + rnorm(nrow(lat) * ncol(lat), sd = 0.1 * max(apply(lat, 1, sd))))
  list(Y = Y)
}

K_LAT <- 3L; NPT <- 5L; run_seed <- 71641
CELLS <- expand.grid(rep = 1:3, N_comp = c(2, 3), sim = 0.9)
SEEDS <- c(42, 7, 100, 500, 1234, 2024)
ARMS <- c("clean_4", "ig_4")  # known-bad seed 42 + others: does inv-gamma PREVENT the trap clean shows?
GRID <- expand.grid(i = c(1, 4), arm = ARMS, s = SEEDS, stringsAsFactors = FALSE)
SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
ROWDIR <- file.path(SP, "ex2_igtrap_rows"); dir.create(ROWDIR, showWarnings = FALSE)  # crash-safe: one file per fit

fit_one <- function(row) {
  i <- GRID$i[row]; arm <- GRID$arm[row]; s <- GRID$s[row]
  nch <- if (arm == "clean_8") 8L else 4L
  ig <- (arm == "ig_4")
  tryCatch({
    set.seed(run_seed + 1000 * CELLS$rep[i] + 10 * CELLS$N_comp[i])
    N_comp <- CELLS$N_comp[i]
    gen <- sim_model_intercepts(N_unc = 5 - N_comp, N_comp_true = 2, N_comp_spur = N_comp, K_unc = 1, sim = 0.9, T_times = 30)
    fit_ys <- gen$Y[, anchor_order(gen$Y, K_LAT)]; N_ <- ncol(fit_ys)
    assign("ife_mod", if (ig) ife_ig else ife_clean, envir = globalenv())
    f <- sample_model(N_units = N_, T_times = nrow(fit_ys), K_latent = K_LAT,
      overall_scales = apply(fit_ys, 2, sd), err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.05,
      data = fit_ys, autocor_a = 90, autocor_b = 10, nonstationary = FALSE, num_treated = NPT,
      include_ints = TRUE, int_scale = 3, int_loc = 4, fit_scales = FALSE, alpha_diag = if (ig) 10 else 0,
      pathfinder_init = TRUE, type = "posterior", quiet = TRUE, ad = 0.8, iter = 500, iter_warm = 500,
      max_treedepth = 12, n_chains = nch, seed = s)
    sd_ <- f$sampler_diag
    df <- data.frame(i = i, N_comp = N_comp, arm = arm, nch = nch, seed = s, rhat = round(sd_$rhat_max, 2),
      offmode = sd_$n_offmode, lp_gap = round(sd_$lp_gap_max, 1), div = sd_$n_div, row.names = NULL)
    saveRDS(df, file.path(ROWDIR, sprintf("row_%02d.rds", row)))
    df
  }, error = function(e) {
    df <- data.frame(i = i, N_comp = NA, arm = arm, nch = nch, seed = s, rhat = NA, offmode = NA,
      lp_gap = NA, div = NA, row.names = NULL)
    saveRDS(df, file.path(ROWDIR, sprintf("row_%02d.rds", row))); df
  })
}

cl <- makeCluster(5, outfile = "")
registerDoParallel(cl); invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, {
  source("../sample_model.r"); source("../pathfinder_init.r")
  ife_clean <- cmdstan_model(stan_file = "../ife_named_unitnorm.stan")
  ife_ig    <- cmdstan_model(stan_file = "../ife_named_unitnorm_ig.stan")
}))
clusterExport(cl, c("fit_one", "ruv", "anchor_order", "sim_model_intercepts", "CELLS", "K_LAT", "NPT", "run_seed", "GRID", "ROWDIR"))
res <- tryCatch(foreach(row = seq_len(nrow(GRID)), .combine = "rbind",
  .packages = c("cmdstanr", "posterior")) %dopar% fit_one(row), error = function(e) NULL)
tryCatch(stopCluster(cl), error = function(e) NULL)
# Fallback: reconstruct from the per-fit row files (survives a combine/cluster stall)
if (is.null(res) || nrow(res) < nrow(GRID)) {
  rfs <- list.files(ROWDIR, pattern = "^row_.*\\.rds$", full.names = TRUE)
  res <- do.call(rbind, lapply(rfs, readRDS))
}

save(res, file = file.path(SP, "ex2_unitnorm_igtrap.RData"))
res <- res[order(res$i, match(res$arm, ARMS), res$seed), ]
cat("\nex2 unit-norm | REPLICATION across seeds | clean_4 vs ig_4 vs clean_8 | 2 trapping datasets (i1,i4) x 8 seeds\n")
cat("If clean_4 spikes stochastically and ig_4 / clean_8 drive spike-rate to ~0 -> shallow minor modes, and the fix.\n\n")
print(res, row.names = FALSE)
cat("\nsummary (spike = rhat>1.1):\n")
for (ii in c(1, 4)) for (a in ARMS) { r <- res[res$i == ii & res$arm == a, ]
  if (nrow(r) == 0) next
  cat(sprintf("  i%d %-8s: rhat med=%.2f max=%.2f | spike-rate=%d/%d | mean offmode=%.1f | mean div=%.0f\n",
    ii, a, median(r$rhat, na.rm = TRUE), max(r$rhat, na.rm = TRUE), sum(r$rhat > 1.1, na.rm = TRUE),
    sum(!is.na(r$rhat)), mean(r$offmode, na.rm = TRUE), mean(r$div, na.rm = TRUE))) }
cat("\nRead: clean_4 spike-rate>0 (stochastic trapping)? does ig_4 kill it (zero-avoidance removes minor modes)? does clean_8 kill it?\n")
