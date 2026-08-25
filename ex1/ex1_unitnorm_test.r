## Probe (NOT shipped): unit-norm-rows Lambda variant (ife_named_unitnorm.stan). sigma now carries
## per-unit magnitude, Lambda is direction-only. ESTIMATE sigma (the point of the variant), NO
## pathfinder, extra chains -> does it sample cleanly intrinsically, does the non-identified row
## magnitude cause trouble, and is alpha_diag still needed? Toggle alpha_diag in {10, 0}. stat_weak
## config + nonstat ref, same 6 datasets as ex1_sigma_decomp. Compare science (absz/eabs/tau/cover)
## to that test's fixed-sigma arms (nonstat absz~0.61, weak_rms 0.85, weak_2rms 0.81, weak_est 0.87).

library(foreach); library(doParallel)
source("../sample_model.r"); source("../pathfinder_init.r")
UNM <- cmdstan_model(stan_file = "../ife_named_unitnorm.stan")  # unit-norm variant (pre-compiled)
ife_mod <- UNM                                                  # sample_model reads global ife_mod

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
  nonstationary = TRUE, alpha = 10, m_tau = 2, s_tau = 2, nch = 4L, iter = 500L, warm = 400L)
weak_cfg <- function(ys, alpha) list(
  overall_scales = apply(ys, 2, rms),          # prior scale for the (now-identified, estimated) sigma
  err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.1, autocor_a = 97, autocor_b = 3,
  nonstationary = FALSE, alpha = alpha, m_tau = 0.1, s_tau = 0.1, nch = 6L, iter = 500L, warm = 500L)

arm_cfg <- function(ys, arm) switch(arm,
  nonstat  = nonstat_cfg(ys),
  weak_a10 = weak_cfg(ys, 10),
  weak_a0  = weak_cfg(ys, 0))

fit_cell <- function(ds, arm) {
  tryCatch({
    ys <- test_data$ys[ds, , ]
    cfg <- arm_cfg(ys, arm)
    m_tau <- cfg$m_tau; s_tau <- cfg$s_tau; nch <- cfg$nch; it <- cfg$iter; wm <- cfg$warm; alpha <- cfg$alpha
    cfg$m_tau <- cfg$s_tau <- cfg$nch <- cfg$iter <- cfg$warm <- cfg$alpha <- NULL
    args <- c(list(N_units = ncol(ys), T_times = nrow(ys), K_latent = GEN_K, data = ys, num_treated = 5,
      type = "posterior", quiet = TRUE, fit_scales = TRUE, alpha_diag = alpha, pathfinder_init = FALSE,
      iter = it, iter_warm = wm, n_chains = nch, seed = 42, parallel_chains = 1, max_treedepth = 12,
      return_draws = c("lp__", "tau", "delta", "Lambda", "sigma")), cfg)
    f <- do.call(sample_model, args); d <- f$draws
    tau <- as.numeric(posterior::extract_variable_matrix(d, "tau")); tau_med <- median(tau)
    ptail <- (1 - pnorm(tau_med, m_tau, s_tau)) / (1 - pnorm(0, m_tau, s_tau))
    sig <- mean(sapply(seq_len(ncol(ys)), function(n) mean(as.numeric(posterior::extract_variable_matrix(d, sprintf("sigma[%d]", n))))))
    dm <- sapply(1:5, function(j) as.numeric(posterior::extract_variable_matrix(d, sprintf("delta[%d]", j))))
    per_mean <- colMeans(dm); per_sd <- apply(dm, 2, sd)
    # min normalized diagonal (collapse indicator; unit-norm => in (0,1))
    min_diag <- min(sapply(1:GEN_K, function(k) mean(as.numeric(posterior::extract_variable_matrix(d, sprintf("Lambda[%d,%d]", k, k))))))
    yp <- f$y_pred; T_ <- nrow(ys); N_ <- ncol(ys); inc <- 0
    for (n in 1:N_) for (tt in 1:T_) { b <- quantile(yp[, tt, n], c(.025, .975))
      inc <- inc + (ys[tt, n] >= b[1] && ys[tt, n] <= b[2]) }
    sd_ <- f$sampler_diag
    data.frame(ds = ds, arm = arm, sig = round(sig, 2), tau_med = round(tau_med, 3),
      tau_x_sig = round(tau_med * sig, 2), eabs = round(mean(abs(per_mean)), 3),
      absz = round(mean(abs(per_mean) / per_sd), 3), min_diag = round(min_diag, 3),
      coverage = round(inc / (T_ * N_), 3), rhat = round(sd_$rhat_max, 3), div = sd_$n_div,
      tree = sd_$n_tree, offmode = sd_$n_offmode, lp_gap = round(sd_$lp_gap_max, 1), row.names = NULL)
  }, error = function(e) data.frame(ds = ds, arm = arm, sig = NA, tau_med = NA, tau_x_sig = NA,
    eabs = NA, absz = NA, min_diag = NA, coverage = NA, rhat = NA, div = NA, tree = NA,
    offmode = NA, lp_gap = NA, row.names = NULL))
}

combos <- expand.grid(ds = DS, arm = c("nonstat", "weak_a10", "weak_a0"), stringsAsFactors = FALSE)
cl <- makeCluster(7, outfile = "")
registerDoParallel(cl); invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, {
  source("../sample_model.r"); source("../pathfinder_init.r")
  ife_mod <- cmdstan_model(stan_file = "../ife_named_unitnorm.stan")
}))
clusterExport(cl, c("fit_cell", "arm_cfg", "nonstat_cfg", "weak_cfg", "rms", "test_data", "GEN_K"))
res <- foreach(k = seq_len(nrow(combos)), .combine = "rbind",
  .packages = c("cmdstanr", "posterior")) %dopar% fit_cell(combos$ds[k], combos$arm[k])
stopCluster(cl)

SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
save(res, file = file.path(SP, "ex1_unitnorm_test.RData"))
res <- res[order(match(res$arm, c("nonstat", "weak_a10", "weak_a0")), res$ds), ]
cat("\nex1 UNIT-NORM variant | estimate sigma, NO pathfinder | alpha_diag {10,0} | 6 datasets\n")
cat("(compare science to fixed-sigma decomp: nonstat absz~0.61, weak_rms 0.85, weak_2rms 0.81)\n\n")
print(res, row.names = FALSE)
for (a in c("nonstat", "weak_a10", "weak_a0")) { r <- res[res$arm == a, ]
  cat(sprintf("%-9s: sig=%.2f tau=%.2f tau*sig=%.2f eabs=%.2f absz=%.2f min_diag=%.3f cover=%.3f | div=%d tree=%d rhat=%.3f offmode=%d\n",
    a, mean(r$sig,na.rm=T), mean(r$tau_med,na.rm=T), mean(r$tau_x_sig,na.rm=T), mean(r$eabs,na.rm=T),
    mean(r$absz,na.rm=T), mean(r$min_diag,na.rm=T), mean(r$coverage,na.rm=T),
    sum(r$div,na.rm=T), sum(r$tree,na.rm=T), max(r$rhat,na.rm=T), sum(r$offmode,na.rm=T))) }
cat("\nRead: (sampling) div/tree/rhat/offmode clean despite non-identified magnitude? (alpha) does a0 collapse min_diag->0 vs a10?\n")
cat("(science) absz vs nonstat comparable to the fixed-sigma story? does estimated sigma behave?\n")
