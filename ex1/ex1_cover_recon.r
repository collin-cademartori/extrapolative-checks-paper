## Reconcile coverage (NOT shipped): run stat_weak with the EXACT ex1_sim_study (sampling-fixes)
## config across 12 datasets, capturing tau AND coverage per dataset, to see whether they're
## bimodal/correlated and whether the mean matches the ~0.6-0.75 undercoverage in the full run.
## Coverage computed identically to run_sim_stat. nonstat short (2ch x 500/300) as reference.

library(foreach); library(doParallel)
source("../sample_model.r"); source("../pathfinder_init.r")

GEN_K <- 4L; REPS_REF <- 40L
set.seed(40318)
study_units <- sample.int(2 * REPS_REF, size = REPS_REF, replace = FALSE)
pp_seed <- sample.int(.Machine$integer.max, 1)
test_data <- sample_model(overall_scales = rep(1, 8), err_scale = 3, autocor_a = 8, autocor_b = 2,
  nonstationary = TRUE, num_treated = 0, type = "prior_pred", K_latent = GEN_K, iter = 2 * REPS_REF, seed = pp_seed)
DS <- study_units[1:12]

nonstat_cfg <- function(ys) list(
  overall_scales = apply(ys, 2, function(y) sd(diff(y))),
  err_scale = 0, err_scale_mean = 2, err_scale_sd = 2,
  autocor_a = 8, autocor_b = 2, nonstationary = TRUE,
  fit_scales = FALSE, alpha_diag = 10, pathfinder_init = TRUE, m_tau = 2, s_tau = 2,
  nch = 2L, iter = 500L, warm = 300L)
weak_cfg <- function(ys) list(  # EXACT ex1_sim_study stat_weak (c = 1), full 1500/500 x 4
  overall_scales = apply(ys, 2, sd),
  err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.1,
  autocor_a = 97, autocor_b = 3, nonstationary = FALSE,
  fit_scales = FALSE, alpha_diag = 10, pathfinder_init = TRUE, m_tau = 0.1, s_tau = 0.1,
  nch = 4L, iter = 1500L, warm = 500L)

arm_cfg <- function(ys, arm) if (arm == "nonstat") nonstat_cfg(ys) else weak_cfg(ys)

fit_cell <- function(ds, arm) {
  tryCatch({
    ys <- test_data$ys[ds, , ]
    cfg <- arm_cfg(ys, arm)
    m_tau <- cfg$m_tau; s_tau <- cfg$s_tau; nch <- cfg$nch; it <- cfg$iter; wm <- cfg$warm
    cfg$m_tau <- cfg$s_tau <- cfg$nch <- cfg$iter <- cfg$warm <- NULL
    args <- c(list(N_units = ncol(ys), T_times = nrow(ys), K_latent = GEN_K,
      data = ys, num_treated = 5, type = "posterior", quiet = TRUE,
      iter = it, iter_warm = wm, n_chains = nch, seed = 42, parallel_chains = 1,
      return_draws = c("lp__", "tau", "delta", "sigma")), cfg)
    f <- do.call(sample_model, args)
    d <- f$draws
    tau <- as.numeric(posterior::extract_variable_matrix(d, "tau")); tau_med <- median(tau)
    ptail <- (1 - pnorm(tau_med, m_tau, s_tau)) / (1 - pnorm(0, m_tau, s_tau))
    sig_post <- tryCatch(mean(sapply(seq_len(ncol(ys)), function(n)
      mean(as.numeric(posterior::extract_variable_matrix(d, sprintf("sigma[%d]", n)))))), error = function(e) NA_real_)
    # Coverage: EXACTLY as run_sim_stat -- per (t,n) 95% predictive interval, mean over all t,n.
    yp <- f$y_pred; T_ <- nrow(ys); N_ <- ncol(ys); inc <- 0; wid <- 0
    for (n in 1:N_) { rng <- max(ys[, n]) - min(ys[, n])
      for (tt in 1:T_) { b <- quantile(yp[, tt, n], c(.025, .975))
        inc <- inc + (ys[tt, n] >= b[1] && ys[tt, n] <= b[2]); wid <- wid + (b[2] - b[1]) / rng } }
    sd_ <- f$sampler_diag
    data.frame(ds = ds, arm = arm, sig = round(sig_post, 1), tau_med = round(tau_med, 3),
      prior_tail = signif(ptail, 2), coverage = round(inc / (T_ * N_), 3),
      pred_width = round(wid / (T_ * N_), 2), rhat = round(sd_$rhat_max, 3),
      offmode = sd_$n_offmode, row.names = NULL)
  }, error = function(e) data.frame(ds = ds, arm = arm, sig = NA, tau_med = NA, prior_tail = NA,
    coverage = NA, pred_width = NA, rhat = NA, offmode = NA, row.names = NULL))
}

combos <- expand.grid(ds = DS, arm = c("nonstat", "weak"), stringsAsFactors = FALSE)
cl <- makeCluster(7, outfile = "")
registerDoParallel(cl); invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, { source("../sample_model.r"); source("../pathfinder_init.r") }))
clusterExport(cl, c("fit_cell", "arm_cfg", "nonstat_cfg", "weak_cfg", "test_data", "GEN_K"))
res <- foreach(k = seq_len(nrow(combos)), .combine = "rbind",
  .packages = c("cmdstanr", "posterior")) %dopar% fit_cell(combos$ds[k], combos$arm[k])
stopCluster(cl)

SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
save(res, file = file.path(SP, "ex1_cover_recon.RData"))
w <- res[res$arm == "weak", ]; ns <- res[res$arm == "nonstat", ]
cat(sprintf("\nex1 coverage reconciliation | stat_weak EXACT study config (c=1, 4ch x 1500/500) | %d datasets\n\n", length(DS)))
cat("stat_weak per dataset:\n"); print(w[order(w$coverage), ], row.names = FALSE)
cat(sprintf("\nstat_weak  mean coverage=%.3f  [min %.3f, max %.3f] | mean tau=%.2f | mean pred_width=%.2f\n",
  mean(w$coverage, na.rm = TRUE), min(w$coverage, na.rm = TRUE), max(w$coverage, na.rm = TRUE),
  mean(w$tau_med, na.rm = TRUE), mean(w$pred_width, na.rm = TRUE)))
cat(sprintf("nonstat    mean coverage=%.3f  [min %.3f, max %.3f] | mean tau=%.2f\n",
  mean(ns$coverage, na.rm = TRUE), min(ns$coverage, na.rm = TRUE), max(ns$coverage, na.rm = TRUE),
  mean(ns$tau_med, na.rm = TRUE)))
cat("\nRead: is stat_weak coverage bimodal/low-on-average (=> my earlier 3 were lucky)? Does low coverage track small tau / narrow width?\n")
