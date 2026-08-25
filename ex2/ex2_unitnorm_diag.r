## Diagnostic (NOT shipped): WHY does the unit-norm intercepts model get rhat ~1.7 despite pathfinder
## init? Distinguish (i) minor modes (some chains stuck in lower-lp basins), (ii) co-equal MAJOR modes
## (top-lp chains disagree on identified params = genuine non-identification), (iii) poor mixing (low
## ESS/autocorrelation within a single mode). Per fit, 8 chains, compare rhat over ALL chains vs only
## ON-MODE chains (mean lp within 5 of best), plus ESS incl. the single best chain. Also test whether
## restoring the zero-avoiding inverse-gamma on the raw diagonal (ife_named_unitnorm_ig.stan, alpha=10)
## kills the multimodality. 3 worst datasets from the pathfinder run (i1,i4,i6). Fixed sigma, ints, ad=0.8.

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
DS_I <- c(1, 4, 6)   # worst rhat datasets from the pathfinder fix run

fit_cell <- function(i, model) {
  tryCatch({
    set.seed(run_seed + 1000 * CELLS$rep[i] + 10 * CELLS$N_comp[i])   # identical dataset to prior runs
    N_comp <- CELLS$N_comp[i]
    gen <- sim_model_intercepts(N_unc = 5 - N_comp, N_comp_true = 2, N_comp_spur = N_comp, K_unc = 1, sim = 0.9, T_times = 30)
    test_ys <- gen$Y; perm <- anchor_order(test_ys, K_LAT); fit_ys <- test_ys[, perm]
    N_ <- ncol(fit_ys)
    assign("ife_mod", if (model == "ig") ife_ig else ife_clean, envir = globalenv())
    ad_arg <- if (model == "ig") 10 else 0
    args <- list(N_units = N_, T_times = nrow(fit_ys), K_latent = K_LAT,
      overall_scales = apply(fit_ys, 2, sd), err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.05,
      data = fit_ys, autocor_a = 90, autocor_b = 10, nonstationary = FALSE, num_treated = NPT,
      include_ints = TRUE, int_scale = 3, int_loc = 4,
      fit_scales = FALSE, alpha_diag = ad_arg, pathfinder_init = TRUE, type = "posterior", quiet = TRUE,
      ad = 0.8, iter = 500, iter_warm = 500, max_treedepth = 12, n_chains = 8, seed = 42,
      return_draws = c("lp__", "tau", "delta", "Lambda", "sigma", "gamma"))
    f <- do.call(sample_model, args); d <- f$draws
    vars <- setdiff(dimnames(d)$variable, "lp__")

    lpm <- extract_variable_matrix(d, "lp__"); lpc <- colMeans(lpm)      # per-chain mean lp
    best <- max(lpc); gap <- best - lpc; onmode <- which(gap <= 5); n_on <- length(onmode)
    ss_all <- summarise_draws(subset_draws(d, variable = vars), "rhat", "ess_bulk", "ess_tail")
    rhat_all <- max(ss_all$rhat, na.rm = TRUE)
    essb_all <- min(ss_all$ess_bulk, na.rm = TRUE); esst_all <- min(ss_all$ess_tail, na.rm = TRUE)
    # on-mode-only rhat/ess: if this is clean, the all-chain rhat is minor-mode contamination
    if (n_on >= 2) {
      ss_on <- summarise_draws(subset_draws(d, variable = vars, chain = onmode), "rhat", "ess_bulk", "ess_tail")
      rhat_on <- max(ss_on$rhat, na.rm = TRUE); essb_on <- min(ss_on$ess_bulk, na.rm = TRUE)
    } else { rhat_on <- NA; essb_on <- NA }
    # single best chain: ESS here isolates within-mode mixing (autocorrelation), no between-chain effect
    best_ch <- which.max(lpc)
    ss_best <- summarise_draws(subset_draws(d, variable = c("tau", "delta[1]"), chain = best_ch), "ess_bulk", "ess_tail")
    essb_best <- min(ss_best$ess_bulk, na.rm = TRUE)
    # which identified var drives the all-chain rhat, and does it disagree among on-mode chains?
    drv <- ss_all$variable[which.max(ss_all$rhat)]
    rhat_drv_on <- if (n_on >= 2) as.numeric(summarise_draws(subset_draws(d, variable = drv, chain = onmode), "rhat")$rhat) else NA
    sd_ <- f$sampler_diag
    data.frame(i = i, N_comp = N_comp, model = model,
      rhat_all = round(rhat_all, 2), rhat_on = round(rhat_on, 2),
      essb_all = round(essb_all), essb_on = round(essb_on), essb_best = round(essb_best), esst_all = round(esst_all),
      n_on = n_on, offmode = sum(gap > 5), lp_gap = round(max(gap), 1),
      drv = drv, rhat_drv_on = round(rhat_drv_on, 2),
      lpc_str = paste(sprintf("%.0f", sort(lpc, decreasing = TRUE)), collapse = ","),
      div = sd_$n_div, rhat_max_diag = round(sd_$rhat_max, 2), row.names = NULL)
  }, error = function(e) data.frame(i = i, N_comp = CELLS$N_comp[i], model = model, rhat_all = NA,
    rhat_on = NA, essb_all = NA, essb_on = NA, essb_best = NA, esst_all = NA, n_on = NA, offmode = NA,
    lp_gap = NA, drv = NA, rhat_drv_on = NA, lpc_str = conditionMessage(e), div = NA, rhat_max_diag = NA, row.names = NULL))
}

combos <- expand.grid(i = DS_I, model = c("clean", "ig"), stringsAsFactors = FALSE)
cl <- makeCluster(6, outfile = "")
registerDoParallel(cl); invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, {
  source("../sample_model.r"); source("../pathfinder_init.r")
  ife_clean <- cmdstan_model(stan_file = "../ife_named_unitnorm.stan")
  ife_ig    <- cmdstan_model(stan_file = "../ife_named_unitnorm_ig.stan")
}))
clusterExport(cl, c("fit_cell", "ruv", "anchor_order", "sim_model_intercepts", "CELLS", "K_LAT", "NPT", "run_seed"))
res <- foreach(k = seq_len(nrow(combos)), .combine = "rbind",
  .packages = c("cmdstanr", "posterior")) %dopar% fit_cell(combos$i[k], combos$model[k])
stopCluster(cl)

SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
save(res, file = file.path(SP, "ex2_unitnorm_diag.RData"))
res <- res[order(res$i, res$model), ]
cat("\nex2 unit-norm MULTIMODALITY DIAGNOSTIC | 8 chains | clean vs inv-gamma(raw diag, a=10) | 3 worst datasets\n")
cat("rhat_all vs rhat_on: if rhat_on clean -> minor-mode contamination; if rhat_on high -> co-equal MAJOR modes.\n")
cat("essb_best low (one good chain) -> slow mixing/autocorrelation, not modes.  lpc_str = per-chain mean lp sorted.\n\n")
print(res[, c("i","N_comp","model","rhat_all","rhat_on","essb_all","essb_on","essb_best","n_on","offmode","lp_gap","drv","rhat_drv_on","div")], row.names = FALSE)
cat("\nper-chain lp means (sorted):\n")
for (r in seq_len(nrow(res))) cat(sprintf("  i%d %-5s: [%s]\n", res$i[r], res$model[r], res$lpc_str[r]))
cat("\nRead: (clean) is rhat_on << rhat_all (minor modes) or still high (major modes)? is essb_best small (mixing)?  (ig) does inv-gamma pull rhat down / raise n_on?\n")
