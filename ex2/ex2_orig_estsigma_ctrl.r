## Probe (NOT shipped): does the CLEAN unit-norm variant make ESTIMATING sigma computationally viable
## in ex2's INTERCEPTS model, and does it change the science? Prior blocker: original model est-sigma
## in ex2 diverged at ad=0.8, needing ad=0.99 -> ~45-min fits. Hypothesis: unit-norm removes the
## sigma-Lambda confound (which fixed ex1 sampling) so est-sigma samples cleanly at ad=0.8. Two arms,
## both the unit-norm model, both ad=0.8, NO pathfinder (matches ex1 clean test; pathfinder may choke
## on the non-identified magnitude): fixed sigma (study baseline) vs estimated sigma. Judge (b) viability
## by div/tree/rhat + WALL-TIME on the est arm; (a) difference by absz / coverage / loc_cor_pval / sigma.
## DGP + config copied from ex2_sim_study.r (not sourceable -- it auto-runs 200 reps).

library(foreach); library(doParallel); library(purrr)
source("../sample_model.r"); source("../pathfinder_init.r")
ife_mod <- cmdstan_model(stan_file = "../ife_named.stan")  # clean unit-norm variant

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
  list(Y = Y, groups = c("treated", rep("true", N_comp_true), rep("spurious", N_comp_spur), rep("uncorrelated", N_unc)))
}

K_LAT <- 3L; NPT <- 5L
CELLS <- expand.grid(rep = 1:3, N_comp = c(2, 3), sim = 0.9)   # 6 datasets, hardest sim
run_seed <- 71641

fit_cell <- function(i, arm) {
  tryCatch({
    set.seed(run_seed + 1000 * CELLS$rep[i] + 10 * CELLS$N_comp[i])   # dataset fixed across arms
    N_comp <- CELLS$N_comp[i]; sim <- CELLS$sim[i]
    gen <- sim_model_intercepts(N_unc = 5 - N_comp, N_comp_true = 2, N_comp_spur = N_comp, K_unc = 1, sim = sim, T_times = 30)
    test_ys <- gen$Y; N_ <- ncol(test_ys); T_ <- nrow(test_ys)
    perm <- anchor_order(test_ys, K_LAT); fit_ys <- test_ys[, perm]
    fit_seed <- 42
    est <- identical(arm, "est")
    args <- list(N_units = ncol(fit_ys), T_times = nrow(fit_ys), K_latent = K_LAT,
      overall_scales = apply(fit_ys, 2, sd), err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.05,
      data = fit_ys, autocor_a = 90, autocor_b = 10, nonstationary = FALSE, num_treated = NPT,
      include_ints = TRUE, int_scale = 3, int_loc = 4,
      fit_scales = est, alpha_diag = 0, pathfinder_init = TRUE, type = "posterior", quiet = TRUE,
      ad = 0.8, iter = 500, iter_warm = 500, max_treedepth = 12, n_chains = 4, seed = fit_seed,
      return_draws = c("tau", "sigma", "Lambda"))
    tm <- system.time(f <- do.call(sample_model, args))["elapsed"]
    d <- f$draws
    tau <- median(as.numeric(posterior::extract_variable_matrix(d, "tau")))
    sig <- mean(sapply(seq_len(N_), function(n) mean(as.numeric(posterior::extract_variable_matrix(d, sprintf("sigma[%d]", n))))))
    min_diag <- min(sapply(1:K_LAT, function(k) mean(as.numeric(posterior::extract_variable_matrix(d, sprintf("Lambda[%d,%d]", k, k))))))
    absz <- mean(abs(f$effect_means / f$effect_sds)); eabs <- mean(abs(f$effect_means))
    yp <- f$y_pred; inc <- 0
    for (n in 1:N_) for (tt in 1:T_) { b <- quantile(yp[, tt, n], c(.005, .995))
      inc <- inc + (fit_ys[tt, n] >= b[1] && fit_ys[tt, n] <= b[2]) }
    sd_ <- f$sampler_diag
    data.frame(i = i, rep = CELLS$rep[i], N_comp = N_comp, sim = sim, arm = arm,
      sig = round(sig, 2), tau = round(tau, 3), min_diag = round(min_diag, 3),
      absz = round(absz, 3), eabs = round(eabs, 3), loc_p = round(f$loc_cor_pval, 3),
      cover = round(inc / (T_ * N_), 3), secs = round(tm, 1),
      div = sd_$n_div, tree = sd_$n_tree, rhat = round(sd_$rhat_max, 3), offmode = sd_$n_offmode, row.names = NULL)
  }, error = function(e) data.frame(i = i, rep = CELLS$rep[i], N_comp = CELLS$N_comp[i], sim = CELLS$sim[i],
    arm = arm, sig = NA, tau = NA, min_diag = NA, absz = NA, eabs = NA, loc_p = NA, cover = NA, secs = NA,
    div = NA, tree = NA, rhat = NA, offmode = NA, row.names = NULL))
}

combos <- expand.grid(i = seq_len(nrow(CELLS)), arm = c("fix", "est"), stringsAsFactors = FALSE)
cl <- makeCluster(7, outfile = "")
registerDoParallel(cl); invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, {
  source("../sample_model.r"); source("../pathfinder_init.r")
  ife_mod <- cmdstan_model(stan_file = "../ife_named.stan")
}))
clusterExport(cl, c("fit_cell", "ruv", "anchor_order", "sim_model_intercepts", "CELLS", "K_LAT", "NPT", "run_seed"))
res <- foreach(k = seq_len(nrow(combos)), .combine = "rbind",
  .packages = c("cmdstanr", "posterior")) %dopar% fit_cell(combos$i[k], combos$arm[k])
stopCluster(cl)

SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
save(res, file = file.path(SP, "ex2_orig_estsigma_ctrl.RData"))
res <- res[order(match(res$arm, c("fix", "est")), res$i), ]
cat("\nex2 INTERCEPTS model | ORIGINAL model (control) | fixed vs ESTIMATED sigma | ad=0.8, no pathfinder | 6 datasets (sim .9, N_comp 2/3)\n")
cat("Viability(b): est arm div/tree/rhat/secs at ad=0.8 (orig model needed ad=.99 ~45min). Difference(a): absz/cover/loc_p/sigma.\n\n")
print(res, row.names = FALSE)
for (a in c("fix", "est")) { r <- res[res$arm == a, ]
  cat(sprintf("%-4s: sig=%.2f tau=%.2f min_diag=%.3f | absz=%.2f eabs=%.2f loc_p=%.3f cover=%.3f | secs(mean/max)=%.0f/%.0f div=%d tree=%d rhat=%.3f offmode=%d\n",
    a, mean(r$sig,na.rm=T), mean(r$tau,na.rm=T), mean(r$min_diag,na.rm=T), mean(r$absz,na.rm=T), mean(r$eabs,na.rm=T),
    mean(r$loc_p,na.rm=T), mean(r$cover,na.rm=T), mean(r$secs,na.rm=T), max(r$secs,na.rm=T),
    sum(r$div,na.rm=T), sum(r$tree,na.rm=T), max(r$rhat,na.rm=T), sum(r$offmode,na.rm=T))) }
cat("\nRead: (viable?) est div~0, tree~0, rhat<1.05, secs comparable to fix at ad=0.8? (difference?) does est shift absz/loc_p/cover, and what sigma does it pick vs fixed sd?\n")
