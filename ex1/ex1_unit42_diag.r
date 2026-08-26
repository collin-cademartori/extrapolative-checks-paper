## Localise the rhat 1.497 on unit 42. The study log says n_offmode = 0 and lp_gap = 3.8, i.e. all four
## chains are in the SAME lp basin yet disagree -- so this is not the minor-mode trapping we chased in
## ex2, and adapt_delta / the inverse-gamma would not touch it. Find WHICH parameter fails to mix, and
## whether the chains disagree on its value (separated modes within a basin) or merely wander (slow
## mixing). Reproduces the study's dataset exactly: same seed stream, err_scale = 2, unit 42.
library(posterior)
source("../sample_model.r"); source("../pathfinder_init.r")
ife_mod <- cmdstan_model("../ife_named_cartesian.stan")
rms <- function(x) sqrt(mean(x^2))

GEN_K <- 4L; REPS <- 30L; NPT <- 5L
set.seed(40318)
study_units <- sample.int(2 * REPS, size = REPS, replace = FALSE)
pp_seed <- sample.int(.Machine$integer.max, 1)
td <- sample_model(overall_scales = rep(1, 8), err_scale = 2, alpha_diag = 10, autocor_a = 8,
  autocor_b = 2, nonstationary = TRUE, num_treated = 0, type = "prior_pred", K_latent = GEN_K,
  iter = 2 * REPS, seed = pp_seed)
stopifnot(42 %in% study_units)
ys <- td$ys[42, , ]; N_ <- ncol(ys); T_ <- nrow(ys)
cat(sprintf("\nunit 42: sd(y)=%.2f RMS(y)=%.2f RMS/sd=%.2f\n", mean(apply(ys,2,sd)), mean(apply(ys,2,rms)),
  mean(apply(ys,2,rms))/mean(apply(ys,2,sd))))

for (sd_i in c(42, 7, 2024)) {
  f <- sample_model(N_units = N_, T_times = T_, K_latent = GEN_K,
    overall_scales = 2 * apply(ys, 2, rms),
    err_scale = 0, err_scale_mean = 0.05, err_scale_sd = 0.05, data = ys,
    autocor_a = 97, autocor_b = 3, nonstationary = FALSE, num_treated = NPT,
    fit_scales = FALSE, alpha_diag = 10, pathfinder_init = TRUE, type = "posterior", quiet = TRUE,
    ad = 0.8, iter = 1000, iter_warm = 500, n_chains = 4, seed = sd_i, parallel_chains = 4,
    return_draws = c("lp__", "tau", "rho", "Lambda", "Phi_innovations", "delta", "scale_est", "sigma"))
  d <- f$draws
  vars <- setdiff(dimnames(d)$variable, "lp__")
  ss <- summarise_draws(subset_draws(d, variable = vars), "rhat", "ess_bulk")
  ss <- ss[order(-ss$rhat), ]
  lp <- colMeans(extract_variable_matrix(d, "lp__"))
  cat(sprintf("\n=== seed %d: rhat_max=%.3f  lp per chain=[%s] ===\n", sd_i, max(ss$rhat, na.rm = TRUE),
    paste(sprintf("%.0f", lp), collapse = ",")))
  cat("worst-mixing parameters:\n")
  for (j in 1:6) cat(sprintf("   %-22s rhat=%.3f ess=%.0f\n", ss$variable[j], ss$rhat[j], ss$ess_bulk[j]))
  # for the worst parameter, show the PER-CHAIN means: separated values => distinct modes within the
  # basin; overlapping but drifting => slow mixing
  w <- ss$variable[1]
  m <- extract_variable_matrix(d, w)
  cat(sprintf("   %s per-chain mean: [%s]  (sd within chain: [%s])\n", w,
    paste(sprintf("%+.2f", colMeans(m)), collapse = ","), paste(sprintf("%.2f", apply(m, 2, sd)), collapse = ",")))
  # is the worst parameter a loading? report its row norm and the group rhat by block
  for (b in c("Lambda", "Phi_innovations", "rho", "scale_est", "tau", "delta")) {
    g <- ss[grepl(sprintf("^%s(\\[|$)", b), ss$variable), ]
    if (nrow(g)) cat(sprintf("   block %-16s max rhat=%.3f (min ess=%.0f)\n", b, max(g$rhat), min(g$ess_bulk)))
  }
}
cat("\nRead: if the worst parameter's per-chain means are SEPARATED, chains found distinct configurations\n")
cat("at similar lp -- a label/rotation-type non-identification. If they overlap, it is slow mixing.\n")
