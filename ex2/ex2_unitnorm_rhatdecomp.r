## Decisive (NOT shipped): separate two explanations for the unit-norm ints rhat~1.7 -- (A) too-few
## chains, vs (B) the rhat METRIC was over different params. The broken runs used sample_model's
## sampler_diag$rhat_max = max over RAW params (Lambda, Phi_innovations, gamma_raw, delta_raw, rho,
## tau_param, Phi_means_param); the "clean" 8-chain diag measured only IDENTIFIED outputs (tau, delta,
## normalized Lambda, sigma, gamma). Here, for the SAME 3 datasets at nch in {4,8}, report BOTH rhats
## AND a per-group max-rhat breakdown, so we see which params carry the 1.7 and whether chains change it.

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
GRID <- expand.grid(i = c(1, 4, 6), nch = c(4L, 8L))
GROUPS <- c("Lambda", "Phi_innovations", "gamma_raw", "delta_raw", "rho", "tau_param")

grp_rhat <- function(d, base) {
  vs <- grep(sprintf("^%s(\\[|$)", base), dimnames(d)$variable, value = TRUE)
  if (length(vs) == 0) return(NA_real_)
  suppressWarnings(max(summarise_draws(subset_draws(d, variable = vs), "rhat")$rhat, na.rm = TRUE))
}

fit_cell <- function(row) {
  i <- GRID$i[row]; nch <- GRID$nch[row]
  tryCatch({
    set.seed(run_seed + 1000 * CELLS$rep[i] + 10 * CELLS$N_comp[i])
    N_comp <- CELLS$N_comp[i]
    gen <- sim_model_intercepts(N_unc = 5 - N_comp, N_comp_true = 2, N_comp_spur = N_comp, K_unc = 1, sim = 0.9, T_times = 30)
    fit_ys <- gen$Y[, anchor_order(gen$Y, K_LAT)]; N_ <- ncol(fit_ys)
    f <- sample_model(N_units = N_, T_times = nrow(fit_ys), K_latent = K_LAT,
      overall_scales = apply(fit_ys, 2, sd), err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.05,
      data = fit_ys, autocor_a = 90, autocor_b = 10, nonstationary = FALSE, num_treated = NPT,
      include_ints = TRUE, int_scale = 3, int_loc = 4, fit_scales = FALSE, alpha_diag = 0,
      pathfinder_init = TRUE, type = "posterior", quiet = TRUE, ad = 0.8, iter = 500, iter_warm = 500,
      max_treedepth = 12, n_chains = nch, seed = 42,
      return_draws = c("lp__", "tau", "delta", "Lambda", "sigma", "gamma",
        "Phi_innovations", "gamma_raw", "delta_raw", "rho", "tau_param"))
    d <- f$draws
    ident <- c(grep("^tau(\\[|$)", dimnames(d)$variable, value = TRUE),
      grep("^delta\\[", dimnames(d)$variable, value = TRUE), grep("^Lambda\\[", dimnames(d)$variable, value = TRUE),
      grep("^sigma\\[", dimnames(d)$variable, value = TRUE), grep("^gamma\\[", dimnames(d)$variable, value = TRUE))
    rhat_ident <- suppressWarnings(max(summarise_draws(subset_draws(d, variable = ident), "rhat")$rhat, na.rm = TRUE))
    gr <- sapply(GROUPS, function(b) grp_rhat(d, b))
    worst_grp <- GROUPS[which.max(gr)]
    lpc <- colMeans(extract_variable_matrix(d, "lp__"))
    out <- data.frame(i = i, N_comp = N_comp, nch = nch, rhat_ident = round(rhat_ident, 2),
      rhat_diag = round(f$sampler_diag$rhat_max, 2), worst_grp = worst_grp, worst_rhat = round(max(gr, na.rm = TRUE), 2),
      offmode = f$sampler_diag$n_offmode, lp_gap = round(f$sampler_diag$lp_gap_max, 1), div = f$sampler_diag$n_div,
      lpc = paste(sprintf("%.0f", sort(lpc, decreasing = TRUE)), collapse = ","), row.names = NULL)
    for (b in GROUPS) out[[paste0("rh_", b)]] <- round(gr[[b]], 2)
    out
  }, error = function(e) data.frame(i = i, N_comp = NA, nch = nch, rhat_ident = NA, rhat_diag = NA,
    worst_grp = conditionMessage(e), worst_rhat = NA, offmode = NA, lp_gap = NA, div = NA, lpc = NA, row.names = NULL))
}

cl <- makeCluster(6, outfile = "")
registerDoParallel(cl); invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, {
  source("../sample_model.r"); source("../pathfinder_init.r")
  ife_mod <- cmdstan_model(stan_file = "../ife_named_unitnorm.stan")
}))
clusterExport(cl, c("fit_cell", "ruv", "anchor_order", "sim_model_intercepts", "CELLS", "K_LAT", "NPT", "run_seed", "GRID", "GROUPS", "grp_rhat"))
res <- foreach(row = seq_len(nrow(GRID)), .combine = "rbind",
  .packages = c("cmdstanr", "posterior")) %dopar% fit_cell(row)
stopCluster(cl)

SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
save(res, file = file.path(SP, "ex2_unitnorm_rhatdecomp.RData"))
res <- res[order(res$i, res$nch), ]
cat("\nex2 unit-norm | rhat DECOMPOSITION: identified-subset vs sampler_diag(raw) vs per-group | nch {4,8} | same 3 datasets, seed 42\n\n")
print(res[, c("i","N_comp","nch","rhat_ident","rhat_diag","worst_grp","worst_rhat","offmode","lp_gap","div")], row.names = FALSE)
cat("\nper-group max rhat:\n")
print(res[, c("i","nch", paste0("rh_", GROUPS))], row.names = FALSE)
cat("\nper-chain lp (sorted):\n")
for (r in seq_len(nrow(res))) cat(sprintf("  i%d nch=%d: [%s]\n", res$i[r], res$nch[r], res$lpc[r]))
cat("\nRead: does rhat_diag stay ~1.7 at BOTH 4 and 8 chains while rhat_ident stays ~1.0? -> it was the RAW-param metric (esp. worst_grp), not chain count.\n")
