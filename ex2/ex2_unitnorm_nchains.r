## Confirm (NOT shipped): is the unit-norm intercepts model's rhat~1.7 a 4-chain / rhat-estimation
## artifact rather than genuine multimodality? Same 3 datasets, same clean model + config + seed as the
## broken 4-chain runs, sweep n_chains in {4, 8} (and a 2nd seed for 4). Report rhat over all chains,
## rhat over on-mode chains, offmode, and per-chain lp so any adrift chain is visible. If 4-chain shows
## high rhat with one lp-adrift chain while 8-chain is clean, the earlier 1.7 was a small-#chains artifact.

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
GRID <- rbind(
  data.frame(i = c(1, 4, 6), nch = 4, fseed = 42),
  data.frame(i = c(1, 4, 6), nch = 8, fseed = 42),
  data.frame(i = c(1, 4, 6), nch = 4, fseed = 7))

fit_cell <- function(row) {
  i <- GRID$i[row]; nch <- GRID$nch[row]; fseed <- GRID$fseed[row]
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
      max_treedepth = 12, n_chains = nch, seed = fseed,
      return_draws = c("lp__", "tau", "delta", "Lambda", "sigma", "gamma"))
    d <- f$draws; vars <- setdiff(dimnames(d)$variable, "lp__")
    lpc <- colMeans(extract_variable_matrix(d, "lp__")); gap <- max(lpc) - lpc; onmode <- which(gap <= 5)
    rhat_all <- max(summarise_draws(subset_draws(d, variable = vars), "rhat")$rhat, na.rm = TRUE)
    rhat_on <- if (length(onmode) >= 2) max(summarise_draws(subset_draws(d, variable = vars, chain = onmode), "rhat")$rhat, na.rm = TRUE) else NA
    data.frame(i = i, N_comp = N_comp, nch = nch, fseed = fseed, rhat_all = round(rhat_all, 2),
      rhat_on = round(rhat_on, 2), n_on = length(onmode), offmode = sum(gap > 5), lp_gap = round(max(gap), 1),
      div = f$sampler_diag$n_div, lpc = paste(sprintf("%.0f", sort(lpc, decreasing = TRUE)), collapse = ","),
      row.names = NULL)
  }, error = function(e) data.frame(i = i, N_comp = NA, nch = nch, fseed = fseed, rhat_all = NA, rhat_on = NA,
    n_on = NA, offmode = NA, lp_gap = NA, div = NA, lpc = conditionMessage(e), row.names = NULL))
}

cl <- makeCluster(6, outfile = "")
registerDoParallel(cl); invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, {
  source("../sample_model.r"); source("../pathfinder_init.r")
  ife_mod <- cmdstan_model(stan_file = "../ife_named_unitnorm.stan")
}))
clusterExport(cl, c("fit_cell", "ruv", "anchor_order", "sim_model_intercepts", "CELLS", "K_LAT", "NPT", "run_seed", "GRID"))
res <- foreach(row = seq_len(nrow(GRID)), .combine = "rbind",
  .packages = c("cmdstanr", "posterior")) %dopar% fit_cell(row)
stopCluster(cl)

SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
save(res, file = file.path(SP, "ex2_unitnorm_nchains.RData"))
res <- res[order(res$i, res$nch, res$fseed), ]
cat("\nex2 unit-norm | n_chains {4,8} x seed | clean model, ints, ad=0.8, same 3 datasets | is rhat~1.7 a 4-chain artifact?\n\n")
print(res[, c("i","N_comp","nch","fseed","rhat_all","rhat_on","n_on","offmode","lp_gap","div")], row.names = FALSE)
cat("\nper-chain lp (sorted):\n")
for (r in seq_len(nrow(res))) cat(sprintf("  i%d nch=%d seed=%d: [%s]\n", res$i[r], res$nch[r], res$fseed[r], res$lpc[r]))
cat("\nRead: does 4-chain show high rhat_all w/ one adrift chain while rhat_on & 8-chain stay ~1.01? -> small-#chains rhat artifact, not multimodality.\n")
