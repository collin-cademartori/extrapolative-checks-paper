## Probe (NOT shipped): does a fixed sigma = RMS(y) = sqrt(mean(y^2)) reproduce the ~2.6x sd(y)
## scale that estimate-sigma wanted in ex1 (and that ex2's no_ints already uses)? If so, RMS is a
## principled UNIFIED fixed-scale rule. stat_weak at overall_scales = sd(y) vs RMS(y) (both fixed),
## nonstat short reference. Same sampling-fixes config otherwise (alpha=10, K=4, ad=0.98 default).
## Reads: does RMS pull tau back into its N(0.1,0.1) prior (prior_tail up) and lift eabs above nonstat,
## like estimate-sigma did (est picked sig~11)? And is RMS/sd ~ 2.6 in the raw data?

library(foreach); library(doParallel)
source("../sample_model.r"); source("../pathfinder_init.r")

GEN_K <- 4L; REPS_REF <- 40L
set.seed(40318)
study_units <- sample.int(2 * REPS_REF, size = REPS_REF, replace = FALSE)
pp_seed <- sample.int(.Machine$integer.max, 1)
test_data <- sample_model(overall_scales = rep(1, 8), err_scale = 3, autocor_a = 8, autocor_b = 2,
  nonstationary = TRUE, num_treated = 0, type = "prior_pred", K_latent = GEN_K, iter = 2 * REPS_REF, seed = pp_seed)
DS <- study_units[1:6]

nonstat_cfg <- function(ys) list(
  overall_scales = apply(ys, 2, function(y) sd(diff(y))),
  err_scale = 0, err_scale_mean = 2, err_scale_sd = 2, autocor_a = 8, autocor_b = 2,
  nonstationary = TRUE, m_tau = 2, s_tau = 2, nch = 2L, iter = 500L, warm = 300L)
weak_cfg <- function(ys, scale_fn) list(
  overall_scales = apply(ys, 2, scale_fn),
  err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.1, autocor_a = 97, autocor_b = 3,
  nonstationary = FALSE, m_tau = 0.1, s_tau = 0.1, nch = 4L, iter = 1000L, warm = 500L)

rms <- function(y) sqrt(mean(y^2))
arm_cfg <- function(ys, arm) switch(arm,
  nonstat  = nonstat_cfg(ys),
  weak_sd  = weak_cfg(ys, sd),
  weak_rms = weak_cfg(ys, rms))

fit_cell <- function(ds, arm) {
  tryCatch({
    ys <- test_data$ys[ds, , ]
    cfg <- arm_cfg(ys, arm)
    m_tau <- cfg$m_tau; s_tau <- cfg$s_tau; nch <- cfg$nch; it <- cfg$iter; wm <- cfg$warm
    cfg$m_tau <- cfg$s_tau <- cfg$nch <- cfg$iter <- cfg$warm <- NULL
    sd_mean <- mean(apply(ys, 2, sd)); rms_mean <- mean(apply(ys, 2, rms))
    args <- c(list(N_units = ncol(ys), T_times = nrow(ys), K_latent = GEN_K,
      data = ys, num_treated = 5, type = "posterior", quiet = TRUE,
      fit_scales = FALSE, alpha_diag = 10, pathfinder_init = TRUE,
      iter = it, iter_warm = wm, n_chains = nch, seed = 42, parallel_chains = 1,
      return_draws = c("lp__", "tau", "delta")), cfg)
    f <- do.call(sample_model, args)
    d <- f$draws
    tau <- as.numeric(posterior::extract_variable_matrix(d, "tau")); tau_med <- median(tau)
    ptail <- (1 - pnorm(tau_med, m_tau, s_tau)) / (1 - pnorm(0, m_tau, s_tau))
    eff <- rowMeans(sapply(1:5, function(j) as.numeric(posterior::extract_variable_matrix(d, sprintf("delta[%d]", j)))))
    q0 <- mean(eff < 0)
    yp <- f$y_pred; T_ <- nrow(ys); N_ <- ncol(ys); inc <- 0
    for (n in 1:N_) for (tt in 1:T_) { b <- quantile(yp[, tt, n], c(.025, .975))
      inc <- inc + (ys[tt, n] >= b[1] && ys[tt, n] <= b[2]) }
    sd_ <- f$sampler_diag
    data.frame(ds = ds, arm = arm, sd_sc = round(sd_mean, 1), rms_sc = round(rms_mean, 1),
      rms_over_sd = round(rms_mean / sd_mean, 2), sig_used = round(mean(cfg$overall_scales), 1),
      tau_med = round(tau_med, 3), prior_tail = signif(ptail, 2),
      eabs = round(mean(abs(eff)), 3), tail0 = round(min(q0, 1 - q0), 3),
      coverage = round(inc / (T_ * N_), 3), rhat = round(sd_$rhat_max, 3),
      div = sd_$n_div, offmode = sd_$n_offmode, row.names = NULL)
  }, error = function(e) data.frame(ds = ds, arm = arm, sd_sc = NA, rms_sc = NA, rms_over_sd = NA,
    sig_used = NA, tau_med = NA, prior_tail = NA, eabs = NA, tail0 = NA, coverage = NA,
    rhat = NA, div = NA, offmode = NA, row.names = NULL))
}

combos <- expand.grid(ds = DS, arm = c("nonstat", "weak_sd", "weak_rms"), stringsAsFactors = FALSE)
cl <- makeCluster(7, outfile = "")
registerDoParallel(cl); invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, { source("../sample_model.r"); source("../pathfinder_init.r") }))
clusterExport(cl, c("fit_cell", "arm_cfg", "nonstat_cfg", "weak_cfg", "rms", "test_data", "GEN_K"))
res <- foreach(k = seq_len(nrow(combos)), .combine = "rbind",
  .packages = c("cmdstanr", "posterior")) %dopar% fit_cell(combos$ds[k], combos$arm[k])
stopCluster(cl)

SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
save(res, file = file.path(SP, "ex1_rms_scale.RData"))
res <- res[order(res$ds, match(res$arm, c("nonstat", "weak_sd", "weak_rms"))), ]
cat(sprintf("\nex1 RMS-vs-sd fixed-scale probe | stat_weak: sd(y) vs RMS(y), fixed | nonstat short ref | %d datasets\n", length(DS)))
cat("(estimate-sigma earlier picked sig~10-12; target ratio ~2.6x sd)\n\n")
print(res, row.names = FALSE)
cat(sprintf("\nRMS/sd ratio (raw data): mean %.2f\n", mean(res$rms_over_sd, na.rm = TRUE)))
for (a in c("nonstat", "weak_sd", "weak_rms")) { r <- res[res$arm == a, ]
  cat(sprintf("%-9s: sig=%.1f  tau=%.2f  prior_tail=%.2g  eabs=%.2f  coverage=%.3f  div(sum)=%d  rhat=%.3f\n",
    a, mean(r$sig_used, na.rm = TRUE), mean(r$tau_med, na.rm = TRUE), mean(r$prior_tail, na.rm = TRUE),
    mean(r$eabs, na.rm = TRUE), mean(r$coverage, na.rm = TRUE), sum(r$div, na.rm = TRUE), mean(r$rhat, na.rm = TRUE))) }
cat("\nRead: is RMS/sd ~ 2.6? Does weak_rms tau drop into the prior (prior_tail up from ~0) and eabs rise above nonstat, cleanly (div~0)?\n")
