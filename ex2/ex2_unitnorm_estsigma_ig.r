## Test (NOT shipped): estimate sigma on the INV-GAMMA unit-norm ints model (trap-free per the igtrap
## test), sweeping adapt_delta in {0.8, 0.95, 0.99}. Questions: (1) where does the sigma posterior sit vs
## the fixed value sigma_data=sd (multiplier sigma_raw; 1 = matches fixed)? (2) does it place mass NEAR
## ZERO -- the funnel that drives divergences (prior sigma_raw~half-normal(0,5) has its MODE at 0)?
## (3) does tree depth / leapfrog cost shoot up (the runtime risk)? (4) can higher adapt_delta kill the
## divergences at acceptable cost? Crash-safe: one row file per fit. 3 datasets (i1,i4,i6).

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
GRID <- expand.grid(i = c(1, 4, 6), ad = c(0.8, 0.95, 0.99))
SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
ROWDIR <- file.path(SP, "ex2_estig_rows"); dir.create(ROWDIR, showWarnings = FALSE)

fit_one <- function(row) {
  i <- GRID$i[row]; ad <- GRID$ad[row]
  wr <- function(df) { saveRDS(df, file.path(ROWDIR, sprintf("row_%02d.rds", row))); df }
  tryCatch({
    set.seed(run_seed + 1000 * CELLS$rep[i] + 10 * CELLS$N_comp[i])
    N_comp <- CELLS$N_comp[i]
    gen <- sim_model_intercepts(N_unc = 5 - N_comp, N_comp_true = 2, N_comp_spur = N_comp, K_unc = 1, sim = 0.9, T_times = 30)
    fit_ys <- gen$Y[, anchor_order(gen$Y, K_LAT)]; N_ <- ncol(fit_ys)
    sig_data <- apply(fit_ys, 2, sd)                      # the fixed reference (sigma_raw==1)
    assign("ife_mod", ife_ig, envir = globalenv())
    secs <- system.time(f <- sample_model(N_units = N_, T_times = nrow(fit_ys), K_latent = K_LAT,
      overall_scales = sig_data, err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.05,
      data = fit_ys, autocor_a = 90, autocor_b = 10, nonstationary = FALSE, num_treated = NPT,
      include_ints = TRUE, int_scale = 3, int_loc = 4, fit_scales = TRUE, alpha_diag = 10,
      pathfinder_init = TRUE, type = "posterior", quiet = TRUE, ad = ad, iter = 500, iter_warm = 500,
      max_treedepth = 12, n_chains = 4, seed = 42,
      return_draws = c("tau", "sigma", "sigma_raw")))["elapsed"]
    d <- f$draws
    sigm <- sapply(1:N_, function(n) as.numeric(extract_variable_matrix(d, sprintf("sigma[%d]", n))))
    sraw <- sapply(1:N_, function(n) as.numeric(extract_variable_matrix(d, sprintf("sigma_raw[%d]", n))))
    # posterior vs fixed: multiplier median per unit (1 == fixed); mass near zero per unit
    sraw_med_unit <- apply(sraw, 2, median)
    p_near0 <- apply(sraw, 2, function(x) mean(x < 0.1))          # P(multiplier < 0.1) per unit
    sig_q05 <- apply(sigm, 2, quantile, 0.05)                     # 5th pct of sigma per unit (near-0 tail)
    tau <- median(as.numeric(extract_variable_matrix(d, "tau")))
    # tree depth / leapfrog cost from sampler diagnostics
    td <- as.numeric(extract_variable_matrix(f$sampler_draws, "treedepth__"))
    leap <- as.numeric(extract_variable_matrix(f$sampler_draws, "n_leapfrog__"))
    sd_ <- f$sampler_diag
    wr(data.frame(i = i, N_comp = N_comp, ad = ad,
      fix_sig = round(mean(sig_data), 2), sraw_med = round(median(sraw), 2),
      sraw_med_min = round(min(sraw_med_unit), 2),                # smallest per-unit multiplier
      sig_mean = round(mean(colMeans(sigm)), 2), sig_min_q05 = round(min(sig_q05), 2),
      pnear0_max = round(max(p_near0), 3), tau = round(tau, 3),
      div = sd_$n_div, ntree = sd_$n_tree, td_mean = round(mean(td), 1), td_max = max(td),
      leap_mean = round(mean(leap)), rhat = round(sd_$rhat_max, 2), offmode = sd_$n_offmode,
      secs = round(secs), row.names = NULL))
  }, error = function(e) wr(data.frame(i = i, N_comp = NA, ad = ad, fix_sig = NA, sraw_med = NA,
    sraw_med_min = NA, sig_mean = NA, sig_min_q05 = NA, pnear0_max = NA, tau = NA, div = NA, ntree = NA,
    td_mean = NA, td_max = NA, leap_mean = NA, rhat = NA, offmode = NA, secs = NA, row.names = NULL)))
}

cl <- makeCluster(5, outfile = "")
registerDoParallel(cl); invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, {
  source("../sample_model.r"); source("../pathfinder_init.r")
  ife_ig <- cmdstan_model(stan_file = "../ife_named_unitnorm_ig.stan")
}))
clusterExport(cl, c("fit_one", "ruv", "anchor_order", "sim_model_intercepts", "CELLS", "K_LAT", "NPT", "run_seed", "GRID", "ROWDIR"))
res <- tryCatch(foreach(row = seq_len(nrow(GRID)), .combine = "rbind",
  .packages = c("cmdstanr", "posterior")) %dopar% fit_one(row), error = function(e) NULL)
tryCatch(stopCluster(cl), error = function(e) NULL)
if (is.null(res) || nrow(res) < nrow(GRID)) {
  rfs <- list.files(ROWDIR, pattern = "^row_.*\\.rds$", full.names = TRUE)
  res <- do.call(rbind, lapply(rfs, readRDS))
}
save(res, file = file.path(SP, "ex2_unitnorm_estsigma_ig.RData"))
res <- res[order(res$i, res$ad), ]
cat("\nex2 unit-norm INV-GAMMA | ESTIMATE sigma | adapt_delta sweep {0.8,0.95,0.99} | 3 datasets\n")
cat("fix_sig=fixed sd | sraw_med=median multiplier (1==fixed) | sraw_med_min=smallest per-unit mult | pnear0_max=max P(mult<0.1) (funnel)\n")
cat("sig_min_q05=smallest per-unit 5th-pct sigma (near-0 tail) | td_*=tree depth | leap_mean=leapfrog/iter (cost) | secs=wall time\n\n")
print(res, row.names = FALSE)
cat("\nby adapt_delta:\n")
for (a in c(0.8, 0.95, 0.99)) { r <- res[res$ad == a, ]
  cat(sprintf("  ad=%.2f: div(sum)=%d ntree(sum)=%d | td_max=%d leap_mean=%.0f | pnear0_max=%.3f sraw_med=%.2f | rhat=%.2f secs(mean/max)=%.0f/%.0f\n",
    a, sum(r$div, na.rm = TRUE), sum(r$ntree, na.rm = TRUE), max(r$td_max, na.rm = TRUE), mean(r$leap_mean, na.rm = TRUE),
    max(r$pnear0_max, na.rm = TRUE), mean(r$sraw_med, na.rm = TRUE), max(r$rhat, na.rm = TRUE),
    mean(r$secs, na.rm = TRUE), max(r$secs, na.rm = TRUE))) }
cat("\nRead: does estimated sigma sit near fixed (sraw~1) or drift? is pnear0_max large (funnel->divergences)? does higher ad kill div, and at what td/leap/time cost?\n")
