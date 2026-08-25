## Probe (NOT shipped): does the CLEAN unit-norm stat_weak overfit the OBSERVED trajectories more
## than nonstat -- i.e. is its worse delta estimation attributable to trajectory overfit? Uses the
## known truth: prior_pred returns ys (noisy) AND ys_latent (true signal); a posterior fit returns
## y_means (fitted signal = Y_latent). Overfit => fitted signal hugs the noise: small residual-to-obs,
## large residual-to-truth => overfit index rt/ro climbs. Measured on all conditioned points and, more
## pointedly, the treated unit's (col 1) pre-window whose extrapolation sets the counterfactual/delta.
## Same clean model, config, seeds, 6 datasets as ex1_unitnorm_clean.

library(foreach); library(doParallel)
source("../sample_model.r"); source("../pathfinder_init.r")
ife_mod <- cmdstan_model(stan_file = "../ife_named_unitnorm.stan")  # clean unit-norm variant

GEN_K <- 4L; REPS_REF <- 40L; NPT <- 5L
set.seed(40318)
study_units <- sample.int(2 * REPS_REF, size = REPS_REF, replace = FALSE)
pp_seed <- sample.int(.Machine$integer.max, 1)
test_data <- sample_model(overall_scales = rep(1, 8), err_scale = 3, autocor_a = 8, autocor_b = 2,
  nonstationary = TRUE, num_treated = 0, type = "prior_pred", K_latent = GEN_K, iter = 2 * REPS_REF, seed = pp_seed)
DS <- study_units[1:6]

rms <- function(x) sqrt(mean(x^2))
nonstat_cfg <- function(ys) list(
  overall_scales = apply(ys, 2, function(y) sd(diff(y))),
  err_scale = 0, err_scale_mean = 2, err_scale_sd = 2, autocor_a = 8, autocor_b = 2,
  nonstationary = TRUE, m_tau = 2, s_tau = 2, nch = 4L, iter = 500L, warm = 400L)
weak_cfg <- function(ys) list(
  overall_scales = apply(ys, 2, rms),
  err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.1, autocor_a = 97, autocor_b = 3,
  nonstationary = FALSE, m_tau = 0.1, s_tau = 0.1, nch = 6L, iter = 500L, warm = 500L)
arm_cfg <- function(ys, arm) if (arm == "nonstat") nonstat_cfg(ys) else weak_cfg(ys)

fit_cell <- function(ds, arm) {
  tryCatch({
    ys  <- test_data$ys[ds, , ]          # observed (noisy)  [T x N]
    tru <- test_data$ys_latent[ds, , ]   # TRUE latent signal [T x N]
    cfg <- arm_cfg(ys, arm)
    m_tau <- cfg$m_tau; s_tau <- cfg$s_tau; nch <- cfg$nch; it <- cfg$iter; wm <- cfg$warm
    cfg$m_tau <- cfg$s_tau <- cfg$nch <- cfg$iter <- cfg$warm <- NULL
    args <- c(list(N_units = ncol(ys), T_times = nrow(ys), K_latent = GEN_K, data = ys, num_treated = NPT,
      type = "posterior", quiet = TRUE, fit_scales = TRUE, alpha_diag = 0, pathfinder_init = FALSE,
      iter = it, iter_warm = wm, n_chains = nch, seed = 42, parallel_chains = 1, max_treedepth = 12,
      return_draws = c("lp__", "tau", "delta", "sigma")), cfg)
    f <- do.call(sample_model, args); d <- f$draws
    fit <- apply(f$y_means, c(2, 3), mean)   # posterior-mean fitted signal [T x N]

    T_ <- nrow(ys); N_ <- ncol(ys)
    pre <- seq_len(T_ - NPT); post <- (T_ - NPT + 1):T_
    cond <- matrix(TRUE, T_, N_); cond[post, 1] <- FALSE   # everything the fit conditions on
    # trajectory overfit on all conditioned points
    rt_all <- rms((fit - tru)[cond]); ro_all <- rms((fit - ys)[cond])
    # treated unit (col 1) pre-window: the extrapolation basis for the counterfactual
    rt_tr <- rms(fit[pre, 1] - tru[pre, 1]); ro_tr <- rms(fit[pre, 1] - ys[pre, 1])
    # hidden counterfactual (col 1 treated tail): recovery err of what delta is measured against
    rt_cf <- rms(fit[post, 1] - tru[post, 1])
    # roughness of fitted signal vs truth on treated pre-window (2nd diffs); >1 => tracking noise
    rgh <- rms(diff(diff(fit[pre, 1]))) / rms(diff(diff(tru[pre, 1])))

    tau <- as.numeric(posterior::extract_variable_matrix(d, "tau")); tau_med <- median(tau)
    sig <- mean(sapply(seq_len(N_), function(n) mean(as.numeric(posterior::extract_variable_matrix(d, sprintf("sigma[%d]", n))))))
    dm <- sapply(1:NPT, function(j) as.numeric(posterior::extract_variable_matrix(d, sprintf("delta[%d]", j))))
    per_mean <- colMeans(dm); per_sd <- apply(dm, 2, sd)
    yp <- f$y_pred; inc <- 0
    for (n in 1:N_) for (tt in 1:T_) { b <- quantile(yp[, tt, n], c(.025, .975))
      inc <- inc + (ys[tt, n] >= b[1] && ys[tt, n] <= b[2]) }
    sd_ <- f$sampler_diag
    data.frame(ds = ds, arm = arm, sig = round(sig, 2), tau = round(tau_med, 3),
      rt_all = round(rt_all, 3), ro_all = round(ro_all, 3), ovf_all = round(rt_all / ro_all, 3),
      rt_tr = round(rt_tr, 3), ro_tr = round(ro_tr, 3), ovf_tr = round(rt_tr / ro_tr, 3),
      rgh = round(rgh, 2), rt_cf = round(rt_cf, 3),
      eabs = round(mean(abs(per_mean)), 3), absz = round(mean(abs(per_mean) / per_sd), 3),
      cover = round(inc / (T_ * N_), 3), div = sd_$n_div, rhat = round(sd_$rhat_max, 3), row.names = NULL)
  }, error = function(e) data.frame(ds = ds, arm = arm, sig = NA, tau = NA, rt_all = NA, ro_all = NA,
    ovf_all = NA, rt_tr = NA, ro_tr = NA, ovf_tr = NA, rgh = NA, rt_cf = NA, eabs = NA, absz = NA,
    cover = NA, div = NA, rhat = NA, row.names = NULL))
}

combos <- expand.grid(ds = DS, arm = c("nonstat", "weak"), stringsAsFactors = FALSE)
cl <- makeCluster(7, outfile = "")
registerDoParallel(cl); invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, {
  source("../sample_model.r"); source("../pathfinder_init.r")
  ife_mod <- cmdstan_model(stan_file = "../ife_named_unitnorm.stan")
}))
clusterExport(cl, c("fit_cell", "arm_cfg", "nonstat_cfg", "weak_cfg", "rms", "test_data", "GEN_K", "NPT"))
res <- foreach(k = seq_len(nrow(combos)), .combine = "rbind",
  .packages = c("cmdstanr", "posterior")) %dopar% fit_cell(combos$ds[k], combos$arm[k])
stopCluster(cl)

SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
save(res, file = file.path(SP, "ex1_unitnorm_overfit.RData"))
res <- res[order(match(res$arm, c("nonstat", "weak")), res$ds), ]
cat("\nex1 CLEAN unit-norm | TRAJECTORY OVERFIT vs known truth | est sigma, no pathfinder | 6 datasets\n")
cat("rt=resid-to-TRUTH  ro=resid-to-OBS  ovf=rt/ro (higher=hugs noise=overfit)  rgh=fitted/true roughness  rt_cf=counterfactual recovery err\n\n")
print(res, row.names = FALSE)
for (a in c("nonstat", "weak")) { r <- res[res$arm == a, ]
  cat(sprintf("%-8s: sig=%.2f tau=%.2f | ALL rt=%.2f ro=%.2f ovf=%.2f | TREATED-pre rt=%.2f ro=%.2f ovf=%.2f rgh=%.2f | cf_err=%.2f | eabs=%.2f absz=%.2f cover=%.3f div=%d\n",
    a, mean(r$sig,na.rm=T), mean(r$tau,na.rm=T), mean(r$rt_all,na.rm=T), mean(r$ro_all,na.rm=T), mean(r$ovf_all,na.rm=T),
    mean(r$rt_tr,na.rm=T), mean(r$ro_tr,na.rm=T), mean(r$ovf_tr,na.rm=T), mean(r$rgh,na.rm=T), mean(r$rt_cf,na.rm=T),
    mean(r$eabs,na.rm=T), mean(r$absz,na.rm=T), mean(r$cover,na.rm=T), sum(r$div,na.rm=T))) }
cat("\nRead: does weak have higher ovf (hugs noise) + rgh>1 + larger rt_cf than nonstat? That ties overfit -> bad counterfactual -> spurious delta (eabs/absz).\n")
