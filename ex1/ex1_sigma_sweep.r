## Diagnostic (NOT shipped): fixed-sigma multiplier sweep for stat_weak. sigma = c*sd(y) for
## c in {1,2,4} (fit_scales=FALSE) + an estimate-sigma arm; nonstat as a rough SHORT reference.
## Same datasets/config as ex1_tau_ablation (sampling-fixes config; weak arms 4 chains x 1000/500;
## nonstat trimmed to 2 x 500/300). Captures sigma this time (target input + posterior mean).
## Reads: which sigma pulls tau back into its N(0.1,0.1) prior bulk (prior_tail -> ~0.5) and
## restores the eabs / coverage gap vs nonstat, while sampling stays clean.

library(foreach); library(doParallel)
source("../sample_model.r"); source("../pathfinder_init.r")

GEN_K <- 4L; REPS_REF <- 40L
set.seed(40318)
study_units <- sample.int(2 * REPS_REF, size = REPS_REF, replace = FALSE)
pp_seed <- sample.int(.Machine$integer.max, 1)
test_data <- sample_model(overall_scales = rep(1, 8), err_scale = 3, autocor_a = 8, autocor_b = 2,
  nonstationary = TRUE, num_treated = 0, type = "prior_pred", K_latent = GEN_K, iter = 2 * REPS_REF, seed = pp_seed)
DS <- study_units[1:3]

nonstat_cfg <- function(ys) list(
  overall_scales = apply(ys, 2, function(y) sd(diff(y))),
  err_scale = 0, err_scale_mean = 2, err_scale_sd = 2,
  autocor_a = 8, autocor_b = 2, nonstationary = TRUE,
  fit_scales = FALSE, alpha_diag = 10, pathfinder_init = TRUE, m_tau = 2, s_tau = 2,
  nch = 2L, iter = 500L, warm = 300L)
weak_cfg <- function(ys, c_sigma, est) list(
  overall_scales = c_sigma * apply(ys, 2, sd),
  err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.1,
  autocor_a = 97, autocor_b = 3, nonstationary = FALSE,
  fit_scales = est, alpha_diag = 10, pathfinder_init = TRUE, m_tau = 0.1, s_tau = 0.1,
  nch = 4L, iter = 1000L, warm = 500L)

arm_cfg <- function(ys, arm) switch(arm,
  nonstat  = nonstat_cfg(ys),
  weak_c1  = weak_cfg(ys, 1, FALSE),
  weak_c2  = weak_cfg(ys, 2, FALSE),
  weak_c4  = weak_cfg(ys, 4, FALSE),
  weak_est = weak_cfg(ys, 1, TRUE))

fit_cell <- function(ds, arm) {
  tryCatch({
    ys <- test_data$ys[ds, , ]
    cfg <- arm_cfg(ys, arm)
    m_tau <- cfg$m_tau; s_tau <- cfg$s_tau; nch <- cfg$nch; it <- cfg$iter; wm <- cfg$warm
    cfg$m_tau <- cfg$s_tau <- cfg$nch <- cfg$iter <- cfg$warm <- NULL
    sig_tgt <- mean(cfg$overall_scales)
    args <- c(list(N_units = ncol(ys), T_times = nrow(ys), K_latent = GEN_K,
      data = ys, num_treated = 5, type = "posterior", quiet = TRUE,
      iter = it, iter_warm = wm, n_chains = nch, seed = 42, parallel_chains = 1,
      return_draws = c("lp__", "tau", "delta", "sigma")), cfg)
    t <- system.time(f <- do.call(sample_model, args))
    d <- f$draws
    tau <- as.numeric(posterior::extract_variable_matrix(d, "tau")); tau_med <- median(tau)
    ptail <- (1 - pnorm(tau_med, m_tau, s_tau)) / (1 - pnorm(0, m_tau, s_tau))
    sig_post <- tryCatch(mean(sapply(seq_len(ncol(ys)), function(n)
      mean(as.numeric(posterior::extract_variable_matrix(d, sprintf("sigma[%d]", n)))))), error = function(e) NA_real_)
    eff <- rowMeans(sapply(1:5, function(j) as.numeric(posterior::extract_variable_matrix(d, sprintf("delta[%d]", j)))))
    q0 <- mean(eff < 0)
    yp <- f$y_pred; T_ <- nrow(ys); N_ <- ncol(ys); inc <- 0; wid <- 0
    for (n in 1:N_) { rng <- max(ys[, n]) - min(ys[, n])
      for (tt in 1:T_) { b <- quantile(yp[, tt, n], c(.025, .975))
        inc <- inc + (ys[tt, n] >= b[1] && ys[tt, n] <= b[2]); wid <- wid + (b[2] - b[1]) / rng } }
    sd_ <- f$sampler_diag
    data.frame(ds = ds, arm = arm, secs = round(t[["elapsed"]]),
      sig_tgt = round(sig_tgt, 2), sig_post = round(sig_post, 2),
      tau_med = round(tau_med, 3), prior_tail = signif(ptail, 2),
      eabs = round(mean(abs(eff)), 3), tail0 = round(min(q0, 1 - q0), 3),
      coverage = round(inc / (T_ * N_), 3), pred_width = round(wid / (T_ * N_), 2),
      rhat = round(sd_$rhat_max, 3), offmode = sd_$n_offmode, tree = sd_$n_tree, row.names = NULL)
  }, error = function(e) data.frame(ds = ds, arm = arm, secs = NA, sig_tgt = NA, sig_post = NA,
    tau_med = NA, prior_tail = NA, eabs = NA, tail0 = NA, coverage = NA, pred_width = NA,
    rhat = NA, offmode = NA, tree = NA, row.names = NULL))
}

arms <- c("nonstat", "weak_c1", "weak_c2", "weak_c4", "weak_est")
combos <- expand.grid(ds = DS, arm = arms, stringsAsFactors = FALSE)

cl <- makeCluster(7, outfile = "")
registerDoParallel(cl); invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, { source("../sample_model.r"); source("../pathfinder_init.r") }))
clusterExport(cl, c("fit_cell", "arm_cfg", "nonstat_cfg", "weak_cfg", "test_data", "GEN_K"))
res <- foreach(k = seq_len(nrow(combos)), .combine = "rbind",
  .packages = c("cmdstanr", "posterior")) %dopar% fit_cell(combos$ds[k], combos$arm[k])
stopCluster(cl)

SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
save(res, file = file.path(SP, "ex1_sigma_sweep.RData"))
res <- res[order(res$ds, match(res$arm, arms)), ]
cat(sprintf("\nex1 fixed-sigma sweep | weak: sigma=c*sd(y), 4ch x 1000/500 | nonstat: 2ch x 500/300 (rough ref) | datasets %s\n",
  paste(DS, collapse = ", ")))
cat("sig_tgt = mean input scale (meaningless for weak_est); sig_post = posterior mean sigma over units\n")
cat("prior_tail = P_prior(tau > posterior median) under stat_weak's N_{>0}(0.1,0.1) (~0.5 = prior binding, ~0 = escaped)\n\n")
print(res, row.names = FALSE)
cat("\nRead: as c grows, does tau_med fall and prior_tail rise toward ~0.5, with eabs climbing above nonstat\n")
cat("(intended stat_weak overfitting)? Which c is the sweet spot? Does weak_est land near one of them? Sampling clean (rhat, offmode)?\n")
