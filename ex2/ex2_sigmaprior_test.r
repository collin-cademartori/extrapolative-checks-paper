## Test (NOT shipped): does TIGHTENING the sigma_raw prior recover the ex2 diagnostic while keeping
## estimate-sigma? Hypothesis (user): normal(0,5) lets sigma run to several x the data scale, so no_ints
## can hit unit levels via large near-cancelling factor-mean terms, freeing spurious loadings onto
## factor 1 (high cs_spur). Tightening (normal(0,1), normal(0,0.5)) should pull the multiplier toward 1
## and cs_spur down toward the fixed-sigma value (~0.45). no_ints only, est sigma, ad=0.8. Same datasets.
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
unpermute_untreated <- function(v, perm) { out <- numeric(length(v)); out[perm[-1] - 1] <- v; out }
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
  list(Y = t(lat + rnorm(nrow(lat) * ncol(lat), sd = 0.1 * max(apply(lat, 1, sd)))), groups = c("treated", rep("true", N_comp_true), rep("spurious", N_comp_spur), rep("uncorrelated", N_unc)))
}

K_LAT <- 3L; NPT <- 5L; run_seed <- 88213
GRID <- expand.grid(rep = 1:2, N_comp = c(2, 3), sim = 0.9, s = c("5", "1", "05"), stringsAsFactors = FALSE)
SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
ROWDIR <- file.path(SP, "ex2_sprior_rows"); dir.create(ROWDIR, showWarnings = FALSE)
mods <- list("5" = "../ife_named_unitnorm_ig.stan", "1" = "../ife_named_unitnorm_ig_s1.stan", "05" = "../ife_named_unitnorm_ig_s05.stan")

fit_one <- function(row) {
  rp <- GRID$rep[row]; N_comp <- GRID$N_comp[row]; sim <- GRID$sim[row]; s <- GRID$s[row]
  wr <- function(df) { saveRDS(df, file.path(ROWDIR, sprintf("row_%02d.rds", row))); df }
  tryCatch({
    set.seed(run_seed + 1000 * rp + 10 * N_comp + round(100 * sim))
    gen <- sim_model_intercepts(N_unc = 5 - N_comp, N_comp_true = 2, N_comp_spur = N_comp, K_unc = 1, sim = sim, T_times = 30)
    test_ys <- gen$Y; groups <- gen$groups; N_ <- ncol(test_ys)
    perm <- anchor_order(test_ys, K_LAT); fit_ys <- test_ys[, perm]
    os <- apply(fit_ys, 2, function(x) sqrt(mean(x^2)))     # no_ints: RMS
    assign("ife_mod", get(paste0("ife_", s), envir = globalenv()), envir = globalenv())
    secs <- system.time(f <- sample_model(N_units = N_, T_times = nrow(fit_ys), K_latent = K_LAT, overall_scales = os,
      err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.05, data = fit_ys, autocor_a = 90, autocor_b = 10,
      nonstationary = FALSE, num_treated = NPT, include_factor_means = TRUE,
      fit_scales = TRUE, alpha_diag = 10, pathfinder_init = TRUE, type = "posterior", quiet = TRUE,
      ad = 0.8, iter = 500, iter_warm = 500, max_treedepth = 12, n_chains = 4, seed = 42,
      return_draws = c("tau", "sigma", "sigma_raw")))["elapsed"]
    absz <- mean(abs(f$effect_means / f$effect_sds))
    cs <- unpermute_untreated(f$cor_sq, perm); ug <- groups[2:N_]
    cs_by <- sapply(c("true", "spurious", "uncorrelated"), function(g) if (any(ug == g)) mean(cs[ug == g]) else NA)
    d <- f$draws
    mult <- mean(sapply(1:N_, function(n) mean(as.numeric(extract_variable_matrix(d, sprintf("sigma_raw[%d]", n))))))
    maxmult <- max(sapply(1:N_, function(n) mean(as.numeric(extract_variable_matrix(d, sprintf("sigma_raw[%d]", n))))))
    sd_ <- f$sampler_diag
    wr(data.frame(rep = rp, N_comp = N_comp, sim = sim, s = s, absz = round(absz, 3),
      cs_true = round(cs_by[1], 3), cs_spur = round(cs_by[2], 3), cs_unc = round(cs_by[3], 3),
      mult = round(mult, 2), maxmult = round(maxmult, 2), div = sd_$n_div, rhat = round(sd_$rhat_max, 2), secs = round(secs), row.names = NULL))
  }, error = function(e) wr(data.frame(rep = rp, N_comp = N_comp, sim = sim, s = s, absz = NA, cs_true = NA, cs_spur = NA, cs_unc = NA, mult = NA, maxmult = NA, div = NA, rhat = NA, secs = NA, row.names = NULL)))
}

cl <- makeCluster(4, outfile = "")
registerDoParallel(cl); invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, {
  source("../sample_model.r"); source("../pathfinder_init.r")
  ife_5 <- cmdstan_model(stan_file = "../ife_named_unitnorm_ig.stan")
  ife_1 <- cmdstan_model(stan_file = "../ife_named_unitnorm_ig_s1.stan")
  ife_05 <- cmdstan_model(stan_file = "../ife_named_unitnorm_ig_s05.stan")
}))
clusterExport(cl, c("fit_one", "ruv", "anchor_order", "unpermute_untreated", "sim_model_intercepts", "GRID", "K_LAT", "NPT", "run_seed", "ROWDIR"))
res <- tryCatch(foreach(row = seq_len(nrow(GRID)), .combine = "rbind", .packages = c("cmdstanr", "posterior")) %dopar% fit_one(row), error = function(e) NULL)
tryCatch(stopCluster(cl), error = function(e) NULL)
if (is.null(res) || nrow(res) < nrow(GRID)) { rfs <- list.files(ROWDIR, pattern = "^row_.*\\.rds$", full.names = TRUE); res <- do.call(rbind, lapply(rfs, readRDS)) }
save(res, file = file.path(SP, "ex2_sigmaprior_test.RData"))
res <- res[order(match(res$s, c("5", "1", "05")), res$N_comp, res$rep), ]
cat("\nex2 no_ints | sigma_raw prior scale sweep {5, 1, 0.5} | est sigma | does tightening drop cs_spur toward fixed(~0.45)?\n\n")
print(res, row.names = FALSE)
cat("\n=== by prior scale (no_ints means) ===\n")
for (sc in c("5", "1", "05")) { r <- res[res$s == sc & is.finite(res$cs_spur), ]
  if (nrow(r)) cat(sprintf("  normal(0,%-3s): cs_spur=%.2f cs_true=%.2f cs_unc=%.2f | mult(mean/max)=%.2f/%.2f | absz=%.2f div=%d rhat=%.2f secs=%.0f\n",
    ifelse(sc == "05", "0.5", sc), mean(r$cs_spur), mean(r$cs_true), mean(r$cs_unc), mean(r$mult), mean(r$maxmult), mean(r$absz), sum(r$div), max(r$rhat), mean(r$secs))) }
cat("\n(baselines: est-sigma normal(0,5) no_ints cs_spur~0.55; FIXED sigma no_ints cs_spur~0.45)\n")
cat("Read: does cs_spur fall 0.55 -> ~0.45 and mult fall toward 1 as the prior tightens? -> hypothesis holds, keep est-sigma with a tighter prior.\n")
