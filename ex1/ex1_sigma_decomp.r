## Probe (NOT shipped): decompose WHY sigma changes stat_weak's fit. sigma in {RMS, 2*RMS, estimated}
## + nonstat short ref, on ex1 DGP datasets. Channels:
##  (a) genuine overfit -> standardized effect absz = |dbar|/sd(delta) climbs, Lambda magnitude shifts
##  (b) effect-prior widening (delta prior sd ∝ sigma[1]) -> eabs climbs but absz ~ flat
## Invariants: tau*sigma (absolute error) and sig*|Lambda| (effective signal) ~ data-pinned, ~constant.
## Config = sampling-fixes ex1 stat_weak (alpha=10, K=4, ad=0.98 default, pathfinder). Also finally
## records div for ex1 estimate-sigma.

library(foreach); library(doParallel)
source("../sample_model.r"); source("../pathfinder_init.r")

GEN_K <- 4L; REPS_REF <- 40L
set.seed(40318)
study_units <- sample.int(2 * REPS_REF, size = REPS_REF, replace = FALSE)
pp_seed <- sample.int(.Machine$integer.max, 1)
test_data <- sample_model(overall_scales = rep(1, 8), err_scale = 3, autocor_a = 8, autocor_b = 2,
  nonstationary = TRUE, num_treated = 0, type = "prior_pred", K_latent = GEN_K, iter = 2 * REPS_REF, seed = pp_seed)
DS <- study_units[1:6]

rms <- function(y) sqrt(mean(y^2))
nonstat_cfg <- function(ys) list(
  overall_scales = apply(ys, 2, function(y) sd(diff(y))),
  err_scale = 0, err_scale_mean = 2, err_scale_sd = 2, autocor_a = 8, autocor_b = 2,
  nonstationary = TRUE, est = FALSE, m_tau = 2, s_tau = 2, nch = 2L, iter = 500L, warm = 300L)
weak_cfg <- function(ys, scale_fn, est) list(
  overall_scales = apply(ys, 2, scale_fn),
  err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.1, autocor_a = 97, autocor_b = 3,
  nonstationary = FALSE, est = est, m_tau = 0.1, s_tau = 0.1, nch = 4L, iter = 1000L, warm = 500L)

arm_cfg <- function(ys, arm) switch(arm,
  nonstat   = nonstat_cfg(ys),
  weak_rms  = weak_cfg(ys, rms, FALSE),
  weak_2rms = weak_cfg(ys, function(y) 2 * rms(y), FALSE),
  weak_est  = weak_cfg(ys, rms, TRUE))

fit_cell <- function(ds, arm) {
  tryCatch({
    ys <- test_data$ys[ds, , ]
    cfg <- arm_cfg(ys, arm)
    m_tau <- cfg$m_tau; s_tau <- cfg$s_tau; nch <- cfg$nch; it <- cfg$iter; wm <- cfg$warm; est <- cfg$est
    cfg$m_tau <- cfg$s_tau <- cfg$nch <- cfg$iter <- cfg$warm <- cfg$est <- NULL
    args <- c(list(N_units = ncol(ys), T_times = nrow(ys), K_latent = GEN_K, data = ys, num_treated = 5,
      type = "posterior", quiet = TRUE, fit_scales = if (est) TRUE else FALSE, alpha_diag = 10,
      pathfinder_init = TRUE, iter = it, iter_warm = wm, n_chains = nch, seed = 42, parallel_chains = 1,
      return_draws = c("lp__", "tau", "delta", "Lambda", "sigma")), cfg)
    f <- do.call(sample_model, args); d <- f$draws
    tau <- as.numeric(posterior::extract_variable_matrix(d, "tau")); tau_med <- median(tau)
    ptail <- (1 - pnorm(tau_med, m_tau, s_tau)) / (1 - pnorm(0, m_tau, s_tau))
    sig <- tryCatch(mean(sapply(seq_len(ncol(ys)), function(n)
      mean(as.numeric(posterior::extract_variable_matrix(d, sprintf("sigma[%d]", n)))))),
      error = function(e) mean(cfg$overall_scales))
    # delta: per-period posterior mean/sd -> eabs (raw) and absz (standardized)
    dm <- sapply(1:5, function(j) as.numeric(posterior::extract_variable_matrix(d, sprintf("delta[%d]", j))))
    per_mean <- colMeans(dm); per_sd <- apply(dm, 2, sd)
    eff <- rowMeans(dm); q0 <- mean(eff < 0)
    eabs <- mean(abs(per_mean)); absz <- mean(abs(per_mean) / per_sd); dsd <- mean(per_sd)
    # Lambda magnitude: sqrt(sum_ij E[Lambda_ij^2]) (raw); effective signal ~ sig * lam_norm
    lv <- grep("^Lambda\\[", posterior::variables(d), value = TRUE)
    lam_norm <- sqrt(sum(vapply(lv, function(v) mean(as.numeric(posterior::extract_variable_matrix(d, v))^2), 0)))
    yp <- f$y_pred; T_ <- nrow(ys); N_ <- ncol(ys); inc <- 0
    for (n in 1:N_) for (tt in 1:T_) { b <- quantile(yp[, tt, n], c(.025, .975))
      inc <- inc + (ys[tt, n] >= b[1] && ys[tt, n] <= b[2]) }
    sd_ <- f$sampler_diag
    data.frame(ds = ds, arm = arm, sig = round(sig, 1), tau_med = round(tau_med, 3),
      tau_x_sig = round(tau_med * sig, 2), prior_tail = signif(ptail, 2),
      eabs = round(eabs, 3), absz = round(absz, 3), dsd = round(dsd, 3),
      lam_norm = round(lam_norm, 2), sig_x_lam = round(sig * lam_norm, 1),
      coverage = round(inc / (T_ * N_), 3), rhat = round(sd_$rhat_max, 3),
      div = sd_$n_div, offmode = sd_$n_offmode, tail0 = round(min(q0, 1 - q0), 3), row.names = NULL)
  }, error = function(e) data.frame(ds = ds, arm = arm, sig = NA, tau_med = NA, tau_x_sig = NA,
    prior_tail = NA, eabs = NA, absz = NA, dsd = NA, lam_norm = NA, sig_x_lam = NA, coverage = NA,
    rhat = NA, div = NA, offmode = NA, tail0 = NA, row.names = NULL))
}

combos <- expand.grid(ds = DS, arm = c("nonstat", "weak_rms", "weak_2rms", "weak_est"), stringsAsFactors = FALSE)
cl <- makeCluster(7, outfile = "")
registerDoParallel(cl); invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, { source("../sample_model.r"); source("../pathfinder_init.r") }))
clusterExport(cl, c("fit_cell", "arm_cfg", "nonstat_cfg", "weak_cfg", "rms", "test_data", "GEN_K"))
res <- foreach(k = seq_len(nrow(combos)), .combine = "rbind",
  .packages = c("cmdstanr", "posterior")) %dopar% fit_cell(combos$ds[k], combos$arm[k])
stopCluster(cl)

SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
save(res, file = file.path(SP, "ex1_sigma_decomp.RData"))
res <- res[order(match(res$arm, c("nonstat", "weak_rms", "weak_2rms", "weak_est")), res$ds), ]
cat("\nex1 sigma decomposition | stat_weak sigma in {RMS, 2xRMS, est} + nonstat ref | 6 datasets\n")
cat("channel(a) overfit: absz up, lam_norm/sig_x_lam shift | channel(b) prior-widen: eabs up but absz flat | invariants: tau_x_sig, sig_x_lam ~const\n\n")
print(res, row.names = FALSE)
for (a in c("nonstat", "weak_rms", "weak_2rms", "weak_est")) { r <- res[res$arm == a, ]
  cat(sprintf("%-9s: sig=%.1f  tau=%.2f  tau*sig=%.2f  eabs=%.2f  absz=%.2f  lam=%.2f  sig*lam=%.1f  cover=%.3f  div=%d\n",
    a, mean(r$sig,na.rm=T), mean(r$tau_med,na.rm=T), mean(r$tau_x_sig,na.rm=T), mean(r$eabs,na.rm=T),
    mean(r$absz,na.rm=T), mean(r$lam_norm,na.rm=T), mean(r$sig_x_lam,na.rm=T), mean(r$coverage,na.rm=T), sum(r$div,na.rm=T))) }
cat("\nRead: from RMS->2RMS->est, does absz climb (real overfit) or stay flat while eabs climbs (just prior width)? Is tau*sig constant?\n")
