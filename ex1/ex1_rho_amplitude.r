## Quick ex1 probe: does the POSTERIOR rho drift above its 0.97 prior mean, and -- more directly --
## what is the POSTERIOR realized amplitude of Phi? The required sigma multiple is set by the realized
## amplitude, not by rho per se, so measure it rather than infer it:
##   needed multiple ~ 1 / (realized amplitude of the fitted factors)
## Reports, per fit: posterior mean rho by factor, realized sd and RMS of each posterior Phi path,
## tau, and scale_est/RMS(y). Cartesian model with overall_scales = RMS(y) -- the config that was
## found insufficient. Short: 3 datasets, ad=0.8, 4 chains.
library(posterior)
source("../sample_model.r"); source("../pathfinder_init.r")
ife_mod <- cmdstan_model("../ife_named_cartesian.stan")
stopifnot(basename(ife_mod$stan_file()) == "ife_named_cartesian.stan")

GEN_K <- 4L; REPS <- 30L
set.seed(40318)                      # same stream as ex1_sim_study
study_units <- sample.int(2 * REPS, size = REPS, replace = FALSE)
pp_seed <- sample.int(.Machine$integer.max, 1)
test_data <- sample_model(overall_scales = rep(1, 8), err_scale = 3, alpha_diag = 10,
  autocor_a = 8, autocor_b = 2, nonstationary = TRUE, num_treated = 0, type = "prior_pred",
  K_latent = GEN_K, iter = 2 * REPS, seed = pp_seed)

rms <- function(x) sqrt(mean(x^2))
cat("\nex1 stat_weak | Cartesian, overall_scales = RMS(y) | posterior rho + realized Phi amplitude\n")
cat("(prior rho ~ Beta(97,3), mean 0.970; PRIOR realized amplitude at that rho/T: sd 0.381, RMS 0.875)\n\n")
for (ds in study_units[1:3]) {
  ys <- test_data$ys[ds, , ]; T_ <- nrow(ys); N_ <- ncol(ys)
  os <- apply(ys, 2, rms)
  f <- sample_model(N_units = N_, T_times = T_, K_latent = GEN_K, overall_scales = os,
    err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.1, data = ys,
    autocor_a = 97, autocor_b = 3, nonstationary = FALSE, num_treated = 5,
    fit_scales = FALSE, alpha_diag = 10, pathfinder_init = TRUE, type = "posterior", quiet = TRUE,
    ad = 0.8, iter = 500, iter_warm = 500, n_chains = 4, seed = 42, parallel_chains = 4,
    return_draws = c("tau", "rho", "Phi", "scale_est"))
  d <- f$draws
  rho_m <- sapply(1:GEN_K, function(k) mean(as.numeric(extract_variable_matrix(d, sprintf("rho[%d]", k)))))
  # realized amplitude of the POSTERIOR factor paths (per draw, then averaged)
  amp <- t(sapply(1:GEN_K, function(k) {
    P <- sapply(1:T_, function(t) as.numeric(extract_variable_matrix(d, sprintf("Phi[%d,%d]", t, k))))
    c(sd = mean(apply(P, 1, sd)), rms = mean(apply(P, 1, rms)))
  }))
  tau <- median(as.numeric(extract_variable_matrix(d, "tau")))
  sc <- mean(sapply(1:N_, function(n) mean(as.numeric(extract_variable_matrix(d, sprintf("scale_est[%d]", n))))))
  cat(sprintf("dataset %d:  sd(y)=%.2f RMS(y)=%.2f (RMS/sd=%.2f) | tau=%.3f | scale_est/RMS(y)=%.2f | div=%d rhat=%.2f\n",
    ds, mean(apply(ys, 2, sd)), mean(os), mean(os)/mean(apply(ys, 2, sd)), tau, sc/mean(os),
    f$sampler_diag$n_div, f$sampler_diag$rhat_max))
  for (k in 1:GEN_K) cat(sprintf("   factor %d: rho=%.4f | realized sd=%.3f (=> 1/sd = %.2fx) | realized RMS=%.3f (=> 1/RMS = %.2fx)\n",
    k, rho_m[k], amp[k,"sd"], 1/amp[k,"sd"], amp[k,"rms"], 1/amp[k,"rms"]))
  cat("\n")
}
cat("Read: if posterior rho >> 0.97 and realized sd << 0.381, the SD-based multiple needed is much larger\n")
cat("than the prior-predictive calculation suggested -- explaining why RMS alone was insufficient.\n")
