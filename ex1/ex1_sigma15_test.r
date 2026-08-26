## Gate test: does sigma_data = 1.5 * RMS(y) pull tau into the plausible range (~0.20-0.25) in ex1
## stat_weak? Matched to the 1.0x RMS run (datasets 53, 27, 1) which gave tau = 0.361 / 0.337 / 0.365.
## If tau scaled as 1/S we would expect ~0.24; ||Lambda|| partially absorbs sigma_data changes (in ex2,
## S rose 1.74x while tau moved only -6%), so the realised drop may be smaller. Reports tau, its
## upper-tail probability under the stat_weak prior N(0.1,0.1)+ (the "is this a contradiction" measure),
## the effective scale, and absz/eabs for context.
library(posterior)
source("../sample_model.r"); source("../pathfinder_init.r")
ife_mod <- cmdstan_model("../ife_named_cartesian.stan")
stopifnot(basename(ife_mod$stan_file()) == "ife_named_cartesian.stan")

GEN_K <- 4L; REPS <- 30L; MULT <- 1.5
set.seed(40318)
study_units <- sample.int(2 * REPS, size = REPS, replace = FALSE)
pp_seed <- sample.int(.Machine$integer.max, 1)
test_data <- sample_model(overall_scales = rep(1, 8), err_scale = 3, alpha_diag = 10,
  autocor_a = 8, autocor_b = 2, nonstationary = TRUE, num_treated = 0, type = "prior_pred",
  K_latent = GEN_K, iter = 2 * REPS, seed = pp_seed)
rms <- function(x) sqrt(mean(x^2))
tail_p <- function(t, m = 0.1, s = 0.1) (1 - pnorm(t, m, s)) / (1 - pnorm(0, m, s))

cat(sprintf("\nex1 stat_weak | Cartesian | sigma_data = %.1f x RMS(y) | matched to the 1.0x run\n", MULT))
cat("baseline 1.0x RMS: tau = 0.361 / 0.337 / 0.365  (tail ~0.6-1%%, a clear contradiction of the prior)\n")
cat("target: tau ~ 0.20-0.25  (tail ~8-19%%, inside the bulk)\n\n")
prev <- c(`53` = 0.361, `27` = 0.337, `1` = 0.365)
for (ds in study_units[1:3]) {
  ys <- test_data$ys[ds, , ]; T_ <- nrow(ys); N_ <- ncol(ys)
  os_base <- apply(ys, 2, rms)
  f <- sample_model(N_units = N_, T_times = T_, K_latent = GEN_K, overall_scales = MULT * os_base,
    err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.1, data = ys,
    autocor_a = 97, autocor_b = 3, nonstationary = FALSE, num_treated = 5,
    fit_scales = FALSE, alpha_diag = 10, pathfinder_init = TRUE, type = "posterior", quiet = TRUE,
    ad = 0.8, iter = 500, iter_warm = 500, n_chains = 4, seed = 42, parallel_chains = 4,
    return_draws = c("tau", "scale_est"))
  d <- f$draws
  tau <- median(as.numeric(extract_variable_matrix(d, "tau")))
  sc <- mean(sapply(1:N_, function(n) mean(as.numeric(extract_variable_matrix(d, sprintf("scale_est[%d]", n))))))
  absz <- mean(abs(f$effect_means / f$effect_sds)); eabs <- mean(abs(f$effect_means))
  p0 <- prev[[as.character(ds)]]
  cat(sprintf("ds %-3d: tau=%.3f (was %.3f, %+.0f%%) | prior tail=%.3f | scale_est=%.2f = %.2f x RMS | absz=%.2f eabs=%.2f | div=%d rhat=%.2f\n",
    ds, tau, p0, 100 * (tau / p0 - 1), tail_p(tau), sc, sc / mean(os_base), absz, eabs,
    f$sampler_diag$n_div, f$sampler_diag$rhat_max))
}
cat("\nRead: tau at/below ~0.25 with tail >~8% => worth pushing on. If tau barely moves, ||Lambda|| is\n")
cat("absorbing the change and scale is not the lever for tau (as in ex2).\n")
