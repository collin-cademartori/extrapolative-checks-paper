## Test the global-scale reparam (ife_named_unitnorm_scalereparam.stan) vs the ig baseline on the SAME
## ints dataset as the geometry probe. Did the per-unit scale ridge shrink (max |cor(sigma_i,sigma_j)|)
## and ESS/leapfrog improve? Baseline (ig, ad=0.9, 800/800): secs~356 div~68 leapfrog~824 ESS_min~555
## (gamma_raw) with sigma_raw pairwise r up to 0.98. Same iter/ad/seed; 4 parallel chains.
library(posterior)
source("../sample_model.r"); source("../pathfinder_init.r")
ife_mod <- cmdstan_model(stan_file = "../ife_named_unitnorm_scalereparam.stan")
# Teach pathfinder the reparam's scale parameters (sigma_raw is gone).
PF_PARAM_BASES <- c(setdiff(PF_PARAM_BASES, "sigma_raw"), "log_sigma_g", "log_sigma_dev_raw")

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
sim_model_intercepts <- function(N_unc = 2, N_comp_true = 2, N_comp_spur = 2, T_times = 30, T_treated = 5, K_unc = 1, sim = 0.9) {
  N_units <- 1 + N_comp_true + N_comp_spur + N_unc; K_gen <- 2 + K_unc
  f_treat <- 6 + arima.sim(model = list(ar = 0.9), n = T_times); f_treat_sd <- 1.9
  f_alt <- (f_treat - 6) + c(rep(0, T_times - T_treated), rep(-f_treat_sd, T_treated))
  f_unc <- matrix(nrow = K_unc, ncol = T_times); cor_unc <- Inf
  while (cor_unc > 0.01) { for (k in 1:K_unc) f_unc[k, ] <- rnorm(1, 0, 2) + arima.sim(model = list(ar = 0.9), n = T_times); cor_unc <- max(abs(cor(t(f_unc), t(t(f_treat))))) }
  facs <- rbind(f_treat, f_alt, f_unc); loads <- matrix(nrow = N_units, ncol = K_gen); loads[1, ] <- c(1, rep(0, K_gen - 1))
  for (n in seq_len(N_comp_true)) loads[1 + n, ] <- c(sqrt(sim), 0, sqrt(1 - sim) * ruv(K_gen - 2))
  for (n in seq_len(N_comp_spur)) loads[1 + N_comp_true + n, ] <- c(0, sqrt(sim), sqrt(1 - sim) * ruv(K_gen - 2))
  for (n in seq_len(N_unc)) loads[1 + N_comp_true + N_comp_spur + n, ] <- c(0, 0, ruv(K_gen - 2))
  lat <- loads %*% facs
  list(Y = t(lat + rnorm(nrow(lat) * ncol(lat), sd = 0.1 * max(apply(lat, 1, sd)))))
}

K_LAT <- 3L; NPT <- 5L
set.seed(88213 + 1000 * 1 + 10 * 3 + 90)   # SAME dataset as the geometry probe
gen <- sim_model_intercepts(N_unc = 2, N_comp_true = 2, N_comp_spur = 3, K_unc = 1, sim = 0.9, T_times = 30)
fit_ys <- gen$Y[, anchor_order(gen$Y, K_LAT)]; N_ <- ncol(fit_ys)

secs <- system.time(f <- sample_model(N_units = N_, T_times = nrow(fit_ys), K_latent = K_LAT,
  overall_scales = apply(fit_ys, 2, sd), err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.05,
  data = fit_ys, autocor_a = 90, autocor_b = 10, nonstationary = FALSE, num_treated = NPT,
  include_ints = TRUE, int_scale = 3, int_loc = 4, fit_scales = TRUE, alpha_diag = 10,
  pathfinder_init = TRUE, type = "posterior", quiet = TRUE, ad = 0.9, iter = 800, iter_warm = 800,
  max_treedepth = 12, n_chains = 4, parallel_chains = 4, seed = 42,
  return_draws = c("lp__", "sigma", "log_sigma_g", "log_sigma_dev_raw", "tau", "gamma_raw", "delta_raw", "rho")))["elapsed"]
d <- f$draws; sd_ <- f$sampler_diag
leap <- as.numeric(extract_variable_matrix(f$sampler_draws, "n_leapfrog__"))
td <- as.numeric(extract_variable_matrix(f$sampler_draws, "treedepth__"))
cat(sprintf("\nSCALEREPARAM: secs=%.0f div=%d rhat=%.2f | treedepth mean=%.1f max=%d | leapfrog mean=%.0f\n",
  secs, sd_$n_div, sd_$rhat_max, mean(td), max(td), mean(leap)))
cat("(ig baseline SAME dataset/ad/iter: secs~356 div~68 rhat~1.01 treedepth~9.2/11 leapfrog~824 ESS_min~555)\n\n")

vars <- setdiff(dimnames(d)$variable, "lp__")
ss <- summarise_draws(subset_draws(d, variable = vars), "rhat", "ess_bulk"); ss$base <- sub("\\[.*", "", ss$variable)
cat("=== per-group ESS_bulk ===\n")
agg <- do.call(rbind, lapply(split(ss, ss$base), function(g) data.frame(group = g$base[1], n = nrow(g),
  ess_min = round(min(g$ess_bulk)), ess_med = round(median(g$ess_bulk)), rhat_max = round(max(g$rhat), 3))))
print(agg[order(agg$ess_min), ], row.names = FALSE)

# did the per-unit SIGMA ridge shrink? (compare to baseline sigma_raw pairwise up to 0.98)
Msig <- sapply(1:N_, function(n) as.numeric(extract_variable_matrix(d, sprintf("sigma[%d]", n))))
Csig <- cor(Msig); Csig[upper.tri(Csig, diag = TRUE)] <- NA
cat(sprintf("\nper-unit sigma[i]~sigma[j]: max|r|=%.2f  median|r|=%.2f   (baseline sigma_raw max|r|~0.98)\n",
  max(abs(Csig), na.rm = TRUE), median(abs(Csig), na.rm = TRUE)))
lg <- as.numeric(extract_variable_matrix(d, "log_sigma_g[1]"))
cat(sprintf("log_sigma_g: mean=%.2f sd=%.2f ess=%.0f  (the shared amplitude axis; wide+low-ESS is fine if per-unit sigma decorrelate)\n",
  mean(lg), sd(lg), ess_bulk(matrix(lg, ncol = 4))))
cat("\nRead: did max|r| among sigma drop well below 0.98 and ESS_min rise above ~555 / leapfrog fall below ~824? -> reparam helps.\n")
