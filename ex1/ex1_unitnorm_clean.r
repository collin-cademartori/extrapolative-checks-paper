## Probe (NOT shipped): CLEAN unit-norm model -- std_normal on Lambda_raw as the SOLE loading prior
## (normalized-Lambda priors removed; alpha_diag now vestigial). Estimate sigma, NO pathfinder, 6
## chains. Does it (1) sample cleanly now that the magnitude is identified, (2) give an identified
## sigma near the signal scale ~sd(y), (3) avoid the diagonal collapse WITHOUT any zero-avoidance
## (min_diag)? nonstat short ref. Same 6 datasets as the decomposition.

library(foreach); library(doParallel)
source("../sample_model.r"); source("../pathfinder_init.r")
ife_mod <- cmdstan_model(stan_file = "../ife_named_unitnorm.stan")  # clean unit-norm variant

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
  nonstationary = TRUE, m_tau = 2, s_tau = 2, nch = 4L, iter = 500L, warm = 400L)
weak_cfg <- function(ys) list(
  overall_scales = apply(ys, 2, rms),            # prior scale for the estimated (identified) sigma
  err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.1, autocor_a = 97, autocor_b = 3,
  nonstationary = FALSE, m_tau = 0.1, s_tau = 0.1, nch = 6L, iter = 500L, warm = 500L)
arm_cfg <- function(ys, arm) if (arm == "nonstat") nonstat_cfg(ys) else weak_cfg(ys)

sd_ref <- function(ys) mean(apply(ys, 2, sd))

fit_cell <- function(ds, arm) {
  tryCatch({
    ys <- test_data$ys[ds, , ]
    cfg <- arm_cfg(ys, arm)
    m_tau <- cfg$m_tau; s_tau <- cfg$s_tau; nch <- cfg$nch; it <- cfg$iter; wm <- cfg$warm
    cfg$m_tau <- cfg$s_tau <- cfg$nch <- cfg$iter <- cfg$warm <- NULL
    args <- c(list(N_units = ncol(ys), T_times = nrow(ys), K_latent = GEN_K, data = ys, num_treated = 5,
      type = "posterior", quiet = TRUE, fit_scales = TRUE, alpha_diag = 0, pathfinder_init = FALSE,
      iter = it, iter_warm = wm, n_chains = nch, seed = 42, parallel_chains = 1, max_treedepth = 12,
      return_draws = c("lp__", "tau", "delta", "Lambda", "sigma")), cfg)
    f <- do.call(sample_model, args); d <- f$draws
    tau <- as.numeric(posterior::extract_variable_matrix(d, "tau")); tau_med <- median(tau)
    ptail <- (1 - pnorm(tau_med, m_tau, s_tau)) / (1 - pnorm(0, m_tau, s_tau))
    sig <- mean(sapply(seq_len(ncol(ys)), function(n) mean(as.numeric(posterior::extract_variable_matrix(d, sprintf("sigma[%d]", n))))))
    dm <- sapply(1:5, function(j) as.numeric(posterior::extract_variable_matrix(d, sprintf("delta[%d]", j))))
    per_mean <- colMeans(dm); per_sd <- apply(dm, 2, sd)
    min_diag <- min(sapply(1:GEN_K, function(k) mean(as.numeric(posterior::extract_variable_matrix(d, sprintf("Lambda[%d,%d]", k, k))))))
    yp <- f$y_pred; T_ <- nrow(ys); N_ <- ncol(ys); inc <- 0
    for (n in 1:N_) for (tt in 1:T_) { b <- quantile(yp[, tt, n], c(.025, .975))
      inc <- inc + (ys[tt, n] >= b[1] && ys[tt, n] <= b[2]) }
    sd_ <- f$sampler_diag
    data.frame(ds = ds, arm = arm, sd_y = round(sd_ref(ys), 1), sig = round(sig, 2),
      tau_med = round(tau_med, 3), tau_x_sig = round(tau_med * sig, 2), prior_tail = signif(ptail, 2),
      eabs = round(mean(abs(per_mean)), 3), absz = round(mean(abs(per_mean) / per_sd), 3),
      min_diag = round(min_diag, 3), coverage = round(inc / (T_ * N_), 3), rhat = round(sd_$rhat_max, 3),
      div = sd_$n_div, tree = sd_$n_tree, offmode = sd_$n_offmode, row.names = NULL)
  }, error = function(e) data.frame(ds = ds, arm = arm, sd_y = NA, sig = NA, tau_med = NA, tau_x_sig = NA,
    prior_tail = NA, eabs = NA, absz = NA, min_diag = NA, coverage = NA, rhat = NA, div = NA, tree = NA,
    offmode = NA, row.names = NULL))
}

combos <- expand.grid(ds = DS, arm = c("nonstat", "weak"), stringsAsFactors = FALSE)
cl <- makeCluster(7, outfile = "")
registerDoParallel(cl); invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, {
  source("../sample_model.r"); source("../pathfinder_init.r")
  ife_mod <- cmdstan_model(stan_file = "../ife_named_unitnorm.stan")
}))
clusterExport(cl, c("fit_cell", "arm_cfg", "nonstat_cfg", "weak_cfg", "rms", "sd_ref", "test_data", "GEN_K"))
res <- foreach(k = seq_len(nrow(combos)), .combine = "rbind",
  .packages = c("cmdstanr", "posterior")) %dopar% fit_cell(combos$ds[k], combos$arm[k])
stopCluster(cl)

SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
save(res, file = file.path(SP, "ex1_unitnorm_clean.RData"))
res <- res[order(match(res$arm, c("nonstat", "weak")), res$ds), ]
cat("\nex1 CLEAN unit-norm | std_normal(Lambda_raw) sole loading prior, NO zero-avoidance | est sigma, no pathfinder\n")
cat("(fixed-sigma decomp ref: nonstat absz~0.61 sd_y~4.4 | weak_rms absz 0.85 sig 6.9)\n\n")
print(res, row.names = FALSE)
for (a in c("nonstat", "weak")) { r <- res[res$arm == a, ]
  cat(sprintf("%-8s: sd_y=%.1f sig=%.2f tau=%.2f tau*sig=%.2f eabs=%.2f absz=%.2f min_diag=%.3f cover=%.3f | div=%d tree=%d rhat=%.3f offmode=%d\n",
    a, mean(r$sd_y,na.rm=T), mean(r$sig,na.rm=T), mean(r$tau_med,na.rm=T), mean(r$tau_x_sig,na.rm=T),
    mean(r$eabs,na.rm=T), mean(r$absz,na.rm=T), mean(r$min_diag,na.rm=T), mean(r$coverage,na.rm=T),
    sum(r$div,na.rm=T), sum(r$tree,na.rm=T), max(r$rhat,na.rm=T), sum(r$offmode,na.rm=T))) }
cat("\nRead: (sampling) div~0 rhat~1 now magnitude is identified? (collapse) min_diag healthy w/o zero-avoidance?\n")
cat("(sigma) does estimated sig land near sd_y (identified, ~signal scale) vs the confounded ~11? (science) absz vs nonstat?\n")
