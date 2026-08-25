## Efficiency/geometry probe (NOT shipped): WHY is the unit-norm ints est-sigma fit slow? Fit one ints
## dataset capturing full parameter draws + sampler diagnostics, then locate the compressed directions:
## (i) per-parameter ESS/leapfrog (slowest-mixing = compressed), (ii) strongest pairwise posterior
## correlations (ridges), (iii) the prime suspect -- the non-identified row magnitude ||Lambda_raw[n]||
## (pinned only by prior): its ESS and correlation with sigma_raw[n]/sigma[n]. Also treedepth/leapfrog.
## Runs its 4 chains in PARALLEL (parallel_chains=4) on the spare cores. ad=0.9 (geometry, not precision).

library(posterior)
source("../sample_model.r"); source("../pathfinder_init.r")
ife_mod <- cmdstan_model(stan_file = "../ife_named_unitnorm_ig.stan")

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
  while (cor_unc > 0.01) { for (k in 1:K_unc) f_unc[k, ] <- rnorm(1, 0, 2) + arima.sim(model = list(ar = 0.9), n = T_times); cor_unc <- max(abs(cor(t(f_unc), t(t(f_treat))))) }
  facs <- rbind(f_treat, f_alt, f_unc); loads <- matrix(nrow = N_units, ncol = K_gen)
  loads[1, ] <- c(1, rep(0, K_gen - 1))
  for (n in seq_len(N_comp_true)) loads[1 + n, ] <- c(sqrt(sim), 0, sqrt(1 - sim) * ruv(K_gen - 2))
  for (n in seq_len(N_comp_spur)) loads[1 + N_comp_true + n, ] <- c(0, sqrt(sim), sqrt(1 - sim) * ruv(K_gen - 2))
  for (n in seq_len(N_unc)) loads[1 + N_comp_true + N_comp_spur + n, ] <- c(0, 0, ruv(K_gen - 2))
  lat <- loads %*% facs
  list(Y = t(lat + rnorm(nrow(lat) * ncol(lat), sd = 0.1 * max(apply(lat, 1, sd)))))
}

K_LAT <- 3L; NPT <- 5L
set.seed(88213 + 1000 * 1 + 10 * 3 + 90)          # a sim=0.9, N_comp=3 dataset
gen <- sim_model_intercepts(N_unc = 2, N_comp_true = 2, N_comp_spur = 3, K_unc = 1, sim = 0.9, T_times = 30)
fit_ys <- gen$Y[, anchor_order(gen$Y, K_LAT)]; N_ <- ncol(fit_ys)

secs <- system.time(f <- sample_model(N_units = N_, T_times = nrow(fit_ys), K_latent = K_LAT,
  overall_scales = apply(fit_ys, 2, sd), err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.05,
  data = fit_ys, autocor_a = 90, autocor_b = 10, nonstationary = FALSE, num_treated = NPT,
  include_ints = TRUE, int_scale = 3, int_loc = 4, fit_scales = TRUE, alpha_diag = 10,
  pathfinder_init = TRUE, type = "posterior", quiet = TRUE, ad = 0.9, iter = 800, iter_warm = 800,
  max_treedepth = 12, n_chains = 4, parallel_chains = 4, seed = 42,
  return_draws = c("lp__", "sigma_raw", "sigma", "tau", "Lambda_raw", "gamma_raw", "delta_raw", "rho")))["elapsed"]
d <- f$draws
sd_ <- f$sampler_diag
leap <- as.numeric(extract_variable_matrix(f$sampler_draws, "n_leapfrog__"))
td <- as.numeric(extract_variable_matrix(f$sampler_draws, "treedepth__"))
cat(sprintf("\nFIT: secs=%.0f div=%d rhat=%.2f | treedepth mean=%.1f max=%d | leapfrog mean=%.0f  (ESS/leapfrog = efficiency)\n\n",
  secs, sd_$n_div, sd_$rhat_max, mean(td), max(td), mean(leap)))

# (i) per-parameter ESS: slowest-mixing = most compressed
vars <- setdiff(dimnames(d)$variable, "lp__")
ss <- summarise_draws(subset_draws(d, variable = vars), "rhat", "ess_bulk", "ess_tail")
ss$base <- sub("\\[.*", "", ss$variable)
cat("=== per-parameter-group ESS_bulk (min mixing = compressed direction) ===\n")
agg <- do.call(rbind, lapply(split(ss, ss$base), function(g) data.frame(group = g$base[1], n = nrow(g),
  ess_min = round(min(g$ess_bulk)), ess_med = round(median(g$ess_bulk)), rhat_max = round(max(g$rhat), 3))))
print(agg[order(agg$ess_min), ], row.names = FALSE)
cat(sprintf("\nlowest-ESS single parameters:\n")); print(head(ss[order(ss$ess_bulk), c("variable", "ess_bulk", "rhat")], 8), row.names = FALSE)

# derived: per-row ||Lambda_raw[n,:]|| (the non-identified magnitude) as [iter x chain] matrices
lr_norm <- lapply(1:N_, function(n) {
  ks <- 1:min(K_LAT, n)  # lower-triangular: row n has min(K,n) free entries
  m <- sapply(ks, function(k) as.numeric(extract_variable_matrix(d, sprintf("Lambda_raw[%d,%d]", n, k))))
  matrix(sqrt(rowSums(m^2)), ncol = 1)
})
cat("\n=== non-identified row magnitude ||Lambda_raw[n]||: ESS + corr with sigma_raw[n] (nuisance-direction check) ===\n")
for (n in 1:N_) {
  nm <- as.numeric(lr_norm[[n]]); srn <- as.numeric(extract_variable_matrix(d, sprintf("sigma_raw[%d]", n)))
  sgn <- as.numeric(extract_variable_matrix(d, sprintf("sigma[%d]", n)))
  ess_nm <- tryCatch(round(ess_bulk(matrix(nm, ncol = 4))), error = function(e) NA)
  cat(sprintf("  unit %d: ||Lraw|| mean=%.2f sd=%.2f ess=%s | cor(||Lraw||, sigma_raw)=%.2f cor(||Lraw||, sigma)=%.2f\n",
    n, mean(nm), sd(nm), ess_nm, cor(nm, srn), cor(nm, sgn)))
}

# (ii) strongest pairwise posterior correlations among a curated scalar set
curated <- c(grep("^sigma_raw\\[", vars, value = TRUE), grep("^tau", vars, value = TRUE),
  grep("^gamma_raw\\[", vars, value = TRUE), grep("^delta_raw\\[", vars, value = TRUE),
  grep("^rho\\[", vars, value = TRUE))
M <- sapply(curated, function(v) as.numeric(extract_variable_matrix(d, v)))
for (n in 1:N_) M <- cbind(M, as.numeric(lr_norm[[n]])); colnames(M)[(ncol(M) - N_ + 1):ncol(M)] <- sprintf("normLraw[%d]", 1:N_)
C <- cor(M); C[upper.tri(C, diag = TRUE)] <- NA
idx <- which(abs(C) > 0.6, arr.ind = TRUE)
cat("\n=== strongest posterior correlations |r|>0.6 (ridges = compressed geometry) ===\n")
if (nrow(idx) > 0) { ord <- order(-abs(C[idx])); for (j in head(ord, 20)) cat(sprintf("  %-16s ~ %-16s : r=%.2f\n", rownames(C)[idx[j,1]], colnames(C)[idx[j,2]], C[idx[j,1], idx[j,2]])) } else cat("  (none > 0.6)\n")

SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
saveRDS(list(agg = agg, ss = ss, secs = secs, td = td, leap = leap), file.path(SP, "ex2_geometry.rds"))
cat("\nRead: is ||Lraw|| the low-ESS/compressed nuisance (and correlated with sigma_raw)? which ridges dominate? -> reparam target.\n")
