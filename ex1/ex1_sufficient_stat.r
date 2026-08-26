## Test the conditional-independence argument directly.
## In ife_named (error scale = tau*sigma, NO ||Lambda||), the likelihood depends on (Lambda, Phi) only
## through the product M = Lambda_Phi = Phi * Lambda', and delta's prior involves only sigma[1]/omega_sq.
## Hence  delta  _||_  (Lambda, Phi)  |  M, sigma, tau.
## Consequence: movement WITHIN a fiber {(L,P): L P' = M} cannot affect delta, so poor mixing confined
## to fibers is irrelevant -- but a missed mode carrying a different delta would have to carry a
## different M, so R-hat on M is the diagnostic that certifies delta.
## Sharpest test available: the K=3 fits on ds42, where raw R-hat hit 1.26-1.57. If R-hat on M stays
## clean there, the non-mixing is within-fiber. If M is ALSO bad at K=3 but clean at K=4, the diagnostic
## has teeth -- it separates benign rotation (K=4) from genuinely competing solutions (K=3).
library(posterior)
source("../sample_model.r"); source("../pathfinder_init.r")
ife_mod <- cmdstan_model("../ife_named_cartesian.stan")   # (matches the fits that produced the spikes)
rms <- function(x) sqrt(mean(x^2))

GEN_K <- 4L; REPS <- 30L; NPT <- 5L
set.seed(40318)
study_units <- sample.int(2 * REPS, size = REPS, replace = FALSE)
pp_seed <- sample.int(.Machine$integer.max, 1)
td <- sample_model(overall_scales = rep(1, 8), err_scale = 2, alpha_diag = 10, autocor_a = 8,
  autocor_b = 2, nonstationary = TRUE, num_treated = 0, type = "prior_pred", K_latent = GEN_K,
  iter = 2 * REPS, seed = pp_seed)

blk <- function(d, pat, lab) {
  v <- grep(pat, dimnames(d)$variable, value = TRUE)
  if (!length(v)) return(NULL)
  ss <- summarise_draws(subset_draws(d, variable = v), "rhat", "ess_bulk")
  c(rhat = suppressWarnings(max(ss$rhat, na.rm = TRUE)),
    ess  = suppressWarnings(min(ss$ess_bulk, na.rm = TRUE)))
}

cat("\nex1 | R-hat on the SUFFICIENT STATISTIC M = Lambda_Phi vs on Lambda itself\n")
cat("ds42 is the dataset whose K=3 fits spiked to 1.26-1.57; K=4 is the study's setting.\n\n")
cat(sprintf("%-4s %-3s %-7s | %-18s | %-18s | %-16s | %-10s\n",
  "ds", "K", "seed", "Lambda (raw)", "M = Lambda_Phi", "delta", "tau"))
for (ds in c(42, 27)) {
  ys <- td$ys[ds, , ]; N_ <- ncol(ys); T_ <- nrow(ys)
  for (K in c(3L, 4L)) for (sd_i in c(500, 7, 1234)) {
    f <- tryCatch(sample_model(N_units = N_, T_times = T_, K_latent = K,
      overall_scales = 2 * apply(ys, 2, rms),
      err_scale = 0, err_scale_mean = 0.05, err_scale_sd = 0.05, data = ys,
      autocor_a = 97, autocor_b = 3, nonstationary = FALSE, num_treated = NPT,
      fit_scales = FALSE, alpha_diag = 10, pathfinder_init = TRUE, type = "posterior", quiet = TRUE,
      ad = 0.8, iter = 1000, iter_warm = 500, n_chains = 4, seed = sd_i, parallel_chains = 4,
      return_draws = c("Lambda", "Lambda_Phi", "delta", "tau", "cor_sq")), error = function(e) NULL)
    if (is.null(f)) { cat(sprintf("  !! ds%d K%d seed%d FAILED\n", ds, K, sd_i)); next }
    d <- f$draws
    L  <- blk(d, "^Lambda\\[");      M  <- blk(d, "^Lambda_Phi\\[")
    De <- blk(d, "^delta\\[");       Ta <- blk(d, "^tau$");  Cs <- blk(d, "^cor_sq\\[")
    cat(sprintf("%-4d %-3d %-7d | rhat=%.3f ess=%-5.0f | rhat=%.3f ess=%-5.0f | rhat=%.4f ess=%-5.0f | rhat=%.4f  cor_sq rhat=%.3f\n",
      ds, K, sd_i, L["rhat"], L["ess"], M["rhat"], M["ess"], De["rhat"], De["ess"], Ta["rhat"], Cs["rhat"]))
  }
}
cat("\nRead: if M stays clean where Lambda spikes, the non-mixing is WITHIN-FIBER and provably cannot\n")
cat("affect delta. If M degrades too (expected at K=3, where 3 factors cannot represent 4-factor data\n")
cat("and chains settle on genuinely different compromises), the diagnostic correctly flags a case where\n")
cat("delta inference WOULD be suspect -- which is what makes it a real test rather than a rubber stamp.\n")
