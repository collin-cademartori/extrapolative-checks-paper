## Bracket the sigma_data multiple for ex1 stat_weak: 2.0x and 2.5x RMS (1.0x gave tau~0.35, 1.5x ~0.27),
## with a nonstat reference so the stat_weak-vs-nonstat GAP is visible -- the thing that actually has to
## work. Uses the known latent truth from the prior-predictive DGP to separate the two claims:
##   ovf = resid-to-TRUTH / resid-to-OBS  (higher = fitting the noise)
##   rgh = roughness of fitted signal / roughness of truth (>1 = tracking noise)
##   absz (standardised) vs eabs (raw): eabs also moves with the delta prior (~sigma[1]), absz does not.
library(posterior)
source("../sample_model.r"); source("../pathfinder_init.r")
ife_mod <- cmdstan_model("../ife_named_cartesian.stan")

GEN_K <- 4L; REPS <- 30L; NPT <- 5L
set.seed(40318)
study_units <- sample.int(2 * REPS, size = REPS, replace = FALSE)
pp_seed <- sample.int(.Machine$integer.max, 1)
test_data <- sample_model(overall_scales = rep(1, 8), err_scale = 3, alpha_diag = 10,
  autocor_a = 8, autocor_b = 2, nonstationary = TRUE, num_treated = 0, type = "prior_pred",
  K_latent = GEN_K, iter = 2 * REPS, seed = pp_seed)
rms <- function(x) sqrt(mean(x^2))
tail_p <- function(t, m, s) (1 - pnorm(t, m, s)) / (1 - pnorm(0, m, s))

fit <- function(ys, tru, arm, mult) {
  T_ <- nrow(ys); N_ <- ncol(ys)
  ns <- identical(arm, "nonstat")
  os <- if (ns) apply(ys, 2, function(y) sd(diff(y))) else mult * apply(ys, 2, rms)
  f <- sample_model(N_units = N_, T_times = T_, K_latent = GEN_K, overall_scales = os,
    err_scale = 0, err_scale_mean = if (ns) 2 else 0.1, err_scale_sd = if (ns) 2 else 0.1, data = ys,
    autocor_a = if (ns) 8 else 97, autocor_b = if (ns) 2 else 3, nonstationary = ns, num_treated = NPT,
    fit_scales = FALSE, alpha_diag = 10, pathfinder_init = TRUE, type = "posterior", quiet = TRUE,
    ad = 0.8, iter = 500, iter_warm = 500, n_chains = 4, seed = 42, parallel_chains = 4,
    return_draws = c("tau", "scale_est"))
  d <- f$draws
  tau <- median(as.numeric(extract_variable_matrix(d, "tau")))
  sc <- mean(sapply(1:N_, function(n) mean(as.numeric(extract_variable_matrix(d, sprintf("scale_est[%d]", n))))))
  fitm <- apply(f$y_means, c(2, 3), mean)
  pre <- seq_len(T_ - NPT); cond <- matrix(TRUE, T_, N_); cond[(T_ - NPT + 1):T_, 1] <- FALSE
  ovf <- rms((fitm - tru)[cond]) / rms((fitm - ys)[cond])
  rgh <- rms(diff(diff(fitm[pre, 1]))) / rms(diff(diff(tru[pre, 1])))
  yp <- f$y_pred; inc <- 0
  for (n in 1:N_) for (tt in 1:T_) { b <- quantile(yp[, tt, n], c(.025, .975)); inc <- inc + (ys[tt, n] >= b[1] && ys[tt, n] <= b[2]) }
  data.frame(arm = arm, mult = mult, tau = round(tau, 3),
    tail = round(tail_p(tau, if (ns) 2 else 0.1, if (ns) 2 else 0.1), 3),
    scale_rms = round(sc / mean(apply(ys, 2, rms)), 2), ovf = round(ovf, 3), rgh = round(rgh, 2),
    absz = round(mean(abs(f$effect_means / f$effect_sds)), 3), eabs = round(mean(abs(f$effect_means)), 3),
    cover = round(inc / (T_ * N_), 3), div = f$sampler_diag$n_div, rhat = round(f$sampler_diag$rhat_max, 2))
}

out <- list()
for (ds in study_units[1:3]) {
  ys <- test_data$ys[ds, , ]; tru <- test_data$ys_latent[ds, , ]
  for (cfg in list(c("nonstat", NA), c("stat_weak", 2.0), c("stat_weak", 2.5))) {
    r <- tryCatch(fit(ys, tru, cfg[1], as.numeric(cfg[2])),
      error = function(e) { cat(sprintf("  !! ds %d %s mult=%s FAILED: %s\n", ds, cfg[1], cfg[2],
        substr(conditionMessage(e), 1, 90))); NULL })
    if (is.null(r)) next
    r$ds <- ds; out[[length(out) + 1]] <- r
    cat(sprintf("ds %-3d %-9s mult=%-4s tau=%.3f tail=%.3f scale=%.2fxRMS ovf=%.2f rgh=%.2f absz=%.2f eabs=%.2f cover=%.3f div=%d\n",
      ds, r$arm, ifelse(is.na(r$mult), "-", format(r$mult)), r$tau, r$tail, r$scale_rms, r$ovf, r$rgh, r$absz, r$eabs, r$cover, r$div))
  }
}
res <- do.call(rbind, out)
saveRDS(res, "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad/ex1_bracket.rds")
cat("\n=== means by arm (1.0x RMS gave tau~0.354; 1.5x gave tau~0.275, scale 1.18xRMS) ===\n")
for (a in unique(paste(res$arm, res$mult))) { s <- res[paste(res$arm, res$mult) == a, ]
  cat(sprintf("%-18s tau=%.3f tail=%.3f scale=%.2fxRMS | ovf=%.2f rgh=%.2f | absz=%.2f eabs=%.2f cover=%.3f\n",
    a, mean(s$tau), mean(s$tail), mean(s$scale_rms), mean(s$ovf), mean(s$rgh), mean(s$absz), mean(s$eabs), mean(s$cover))) }
cat("\nRead: does tau reach ~0.20-0.25 (tail >~8%)? And does stat_weak still MIMIC nonstat (similar cover)\n")
cat("while OVERFITTING more (higher ovf/rgh) and estimating delta worse (higher absz, not just eabs)?\n")
