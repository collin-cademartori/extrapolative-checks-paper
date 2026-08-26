## Two questions in one run.
## (1) IDENTIFICATION of rho: refit stat_weak with rho ~ Beta(99,1) (mean .99, sd .01; 0.97 is 2sd away
##     and reachable). If the posterior is dragged to ~0.97 the data identify rho; if it stays ~0.99 rho
##     is weakly identified and merely tracks its prior. Control: Beta(97,3) posteriors were .968-.978.
## (2) A BETTER OVERFITTING METRIC. rgh (roughness of 2nd differences, one unit, pre-window) measures
##     wiggliness and saturates. Since we know the truth, use the canonical quantity -- the fraction of
##     the NOISE absorbed into the fitted signal:
##        noise_abs = cov(fit - truth, obs - truth) / var(obs - truth)
##     = slope of the fit's error on the noise; 1 = interpolation, 0 = noise ignored. For a linear
##     smoother E[noise_abs] = tr(H)/n, i.e. EFFECTIVE DF PER OBSERVATION -- the thing that should drive
##     counterfactual error. Also report gen_gap = err-to-truth on the treated POST window relative to
##     the conditioned points (a generalisation gap). In ex1 the DGP has num_treated = 0, so TRUE delta
##     is 0 and absz/eabs ARE the counterfactual error: we can check which metric tracks them.
library(posterior)
source("../sample_model.r"); source("../pathfinder_init.r")
ife_mod <- cmdstan_model("../ife_named_cartesian.stan")
rms <- function(x) sqrt(mean(x^2))

GEN_K <- 4L; REPS <- 30L; NPT <- 5L
set.seed(40318)
study_units <- sample.int(2 * REPS, size = REPS, replace = FALSE)
pp_seed <- sample.int(.Machine$integer.max, 1)
td <- sample_model(overall_scales = rep(1, 8), err_scale = 3, alpha_diag = 10, autocor_a = 8,
  autocor_b = 2, nonstationary = TRUE, num_treated = 0, type = "prior_pred", K_latent = GEN_K,
  iter = 2 * REPS, seed = pp_seed)

fit1 <- function(ys, tru, arm, mult, a_rho, b_rho) {
  T_ <- nrow(ys); N_ <- ncol(ys); ns <- identical(arm, "nonstat")
  os <- if (ns) apply(ys, 2, function(y) sd(diff(y))) else mult * apply(ys, 2, rms)
  f <- sample_model(N_units = N_, T_times = T_, K_latent = GEN_K, overall_scales = os,
    err_scale = 0, err_scale_mean = if (ns) 2 else 0.1, err_scale_sd = if (ns) 2 else 0.1, data = ys,
    autocor_a = a_rho, autocor_b = b_rho, nonstationary = ns, num_treated = NPT,
    fit_scales = FALSE, alpha_diag = 10, pathfinder_init = TRUE, type = "posterior", quiet = TRUE,
    ad = 0.8, iter = 500, iter_warm = 500, n_chains = 4, seed = 42, parallel_chains = 4,
    return_draws = c("tau", "rho"))
  d <- f$draws
  fitm <- apply(f$y_means, c(2, 3), mean)
  post <- (T_ - NPT + 1):T_; cond <- matrix(TRUE, T_, N_); cond[post, 1] <- FALSE
  e_fit <- (fitm - tru)[cond]; noise <- (ys - tru)[cond]
  noise_abs <- as.numeric(cov(e_fit, noise) / var(noise))          # ~ effective df / n
  gen_gap <- rms(fitm[post, 1] - tru[post, 1]) / rms(e_fit)
  ovf <- rms(e_fit) / rms((fitm - ys)[cond])
  rgh <- rms(diff(diff(fitm[seq_len(T_ - NPT), 1]))) / rms(diff(diff(tru[seq_len(T_ - NPT), 1])))
  data.frame(arm = arm, mult = mult, rho_prior = sprintf("B(%g,%g)", a_rho, b_rho),
    rho = round(mean(sapply(1:GEN_K, function(k) mean(as.numeric(extract_variable_matrix(d, sprintf("rho[%d]", k)))))), 4),
    tau = round(median(as.numeric(extract_variable_matrix(d, "tau"))), 3),
    noise_abs = round(noise_abs, 3), gen_gap = round(gen_gap, 2), ovf = round(ovf, 3), rgh = round(rgh, 2),
    absz = round(mean(abs(f$effect_means / f$effect_sds)), 3), eabs = round(mean(abs(f$effect_means)), 3),
    div = f$sampler_diag$n_div, rhat = round(f$sampler_diag$rhat_max, 2))
}

arms <- list(list("nonstat", NA, 8, 2), list("stat_weak", 1.0, 97, 3),
             list("stat_weak", 2.0, 97, 3), list("stat_weak", 1.5, 99, 1))
out <- list()
for (ds in study_units[1:3]) {
  ys <- td$ys[ds, , ]; tru <- td$ys_latent[ds, , ]
  for (a in arms) {
    r <- tryCatch(fit1(ys, tru, a[[1]], as.numeric(a[[2]]), a[[3]], a[[4]]),
      error = function(e) { cat(sprintf("  !! ds%d %s mult=%s B(%g,%g) FAILED\n", ds, a[[1]], a[[2]], a[[3]], a[[4]])); NULL })
    if (is.null(r)) next
    r$ds <- ds; out[[length(out) + 1]] <- r
    cat(sprintf("ds%-3d %-9s mult=%-4s %-8s rho=%.4f tau=%.3f | noise_abs=%.3f gen_gap=%.2f ovf=%.2f rgh=%.2f | absz=%.2f eabs=%.2f\n",
      ds, r$arm, ifelse(is.na(r$mult), "-", format(r$mult)), r$rho_prior, r$rho, r$tau,
      r$noise_abs, r$gen_gap, r$ovf, r$rgh, r$absz, r$eabs))
  }
}
res <- do.call(rbind, out)
saveRDS(res, "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad/ex1_rho_metrics.rds")
cat("\n=== (1) rho identification ===\n")
for (rp in unique(res$rho_prior)) { s <- res[res$rho_prior == rp & res$arm == "stat_weak", ]
  if (nrow(s)) cat(sprintf("  prior %-8s -> posterior rho = %.4f (n=%d)\n", rp, mean(s$rho), nrow(s))) }
cat("\n=== (2) which overfitting metric tracks delta error (absz)? ===\n")
s <- res[is.finite(res$absz), ]
for (m in c("noise_abs", "gen_gap", "ovf", "rgh")) {
  rr <- suppressWarnings(cor(s[[m]], s$absz))
  cat(sprintf("  cor(%-9s, absz) = %+.2f   [range %.2f - %.2f, spread %.2fx]\n",
    m, rr, min(s[[m]]), max(s[[m]]), max(s[[m]]) / max(min(s[[m]]), 1e-6))) }
cat("\nby arm:\n")
for (a in unique(paste(res$arm, res$mult, res$rho_prior))) { s <- res[paste(res$arm, res$mult, res$rho_prior) == a, ]
  cat(sprintf("  %-28s noise_abs=%.3f gen_gap=%.2f ovf=%.2f rgh=%.2f | absz=%.2f\n",
    a, mean(s$noise_abs), mean(s$gen_gap), mean(s$ovf), mean(s$rgh), mean(s$absz))) }
