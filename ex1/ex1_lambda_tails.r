## Are the slow-mixing Lambda marginals HEAVY-TAILED? That is a common cause of exactly this pattern:
## no divergences, all chains in one lp basin (n_offmode = 0), but low ESS and intermittent high rhat,
## because chains explore a long tail at different rates.
## Diagnostics per Lambda entry:
##   ess_tail vs ess_bulk : the standard signature -- tail ESS collapsing relative to bulk ESS means the
##                          extremes are what the sampler struggles with.
##   qratio = (q97.5-q2.5)/(q75-q25) : 2.91 for a NORMAL; larger = heavier tails. Scale-free.
##   excess kurtosis      : 0 for a normal.
## Sweeps seeds on unit 42 (which gave rhat 1.497 in the study) and keeps the WORST fit for dissection,
## so the claim can be checked on a genuinely bad fit rather than a clean one.
library(posterior)
source("../sample_model.r"); source("../pathfinder_init.r")
ife_mod <- cmdstan_model("../ife_named_cartesian.stan")
rms <- function(x) sqrt(mean(x^2))
SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"

GEN_K <- 4L; REPS <- 30L; NPT <- 5L
set.seed(40318)
study_units <- sample.int(2 * REPS, size = REPS, replace = FALSE)
pp_seed <- sample.int(.Machine$integer.max, 1)
td <- sample_model(overall_scales = rep(1, 8), err_scale = 2, alpha_diag = 10, autocor_a = 8,
  autocor_b = 2, nonstationary = TRUE, num_treated = 0, type = "prior_pred", K_latent = GEN_K,
  iter = 2 * REPS, seed = pp_seed)

qratio <- function(x) { q <- quantile(x, c(.025, .25, .75, .975)); (q[4] - q[1]) / (q[3] - q[2]) }
exkurt <- function(x) { z <- (x - mean(x)) / sd(x); mean(z^4) - 3 }

analyse <- function(d, label) {
  lv <- grep("^Lambda\\[", dimnames(d)$variable, value = TRUE)
  ss <- summarise_draws(subset_draws(d, variable = lv), "rhat", "ess_bulk", "ess_tail")
  ss$qratio <- sapply(lv, function(v) as.numeric(qratio(as.numeric(extract_variable_matrix(d, v)))))
  ss$exkurt <- sapply(lv, function(v) exkurt(as.numeric(extract_variable_matrix(d, v))))
  ss$tail_ratio <- ss$ess_tail / ss$ess_bulk
  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("  Lambda entries: rhat max=%.3f | ess_bulk min=%.0f | ess_tail min=%.0f\n",
    max(ss$rhat), min(ss$ess_bulk), min(ss$ess_tail)))
  cat(sprintf("  ess_tail/ess_bulk: median=%.2f min=%.2f   (<<1 => tails are the bottleneck)\n",
    median(ss$tail_ratio), min(ss$tail_ratio)))
  cat(sprintf("  qratio (normal=2.91): median=%.2f max=%.2f | excess kurtosis: median=%+.2f max=%+.2f\n",
    median(ss$qratio), max(ss$qratio), median(ss$exkurt), max(ss$exkurt)))
  o <- ss[order(ss$ess_bulk), ]
  cat("  slowest-mixing entries:\n")
  for (j in 1:5) cat(sprintf("    %-16s rhat=%.3f ess_bulk=%-5.0f ess_tail=%-5.0f qratio=%.2f exkurt=%+.2f\n",
    o$variable[j], o$rhat[j], o$ess_bulk[j], o$ess_tail[j], o$qratio[j], o$exkurt[j]))
  h <- ss[order(-ss$qratio), ]
  cat("  heaviest-tailed entries:\n")
  for (j in 1:3) cat(sprintf("    %-16s rhat=%.3f ess_bulk=%-5.0f qratio=%.2f exkurt=%+.2f\n",
    h$variable[j], h$rhat[j], h$ess_bulk[j], h$qratio[j], h$exkurt[j]))
  cat(sprintf("  cor(qratio, ess_bulk) = %+.2f   cor(exkurt, ess_bulk) = %+.2f  (negative => heavy tails track slow mixing)\n",
    cor(ss$qratio, ss$ess_bulk), cor(ss$exkurt, ss$ess_bulk)))
  invisible(ss)
}

best_rhat <- 0; best <- NULL; best_seed <- NA
for (sd_i in c(1, 42, 7, 101, 500, 1234, 2024, 8675, 31415, 90210, 55555, 777)) {
  ys <- td$ys[42, , ]; N_ <- ncol(ys); T_ <- nrow(ys)
  f <- tryCatch(sample_model(N_units = N_, T_times = T_, K_latent = GEN_K,
    overall_scales = 2 * apply(ys, 2, rms),
    err_scale = 0, err_scale_mean = 0.05, err_scale_sd = 0.05, data = ys,
    autocor_a = 97, autocor_b = 3, nonstationary = FALSE, num_treated = NPT,
    fit_scales = FALSE, alpha_diag = 10, pathfinder_init = TRUE, type = "posterior", quiet = TRUE,
    ad = 0.8, iter = 1000, iter_warm = 500, n_chains = 4, seed = sd_i, parallel_chains = 4,
    return_draws = c("lp__", "tau", "Lambda", "delta")), error = function(e) NULL)
  if (is.null(f)) next
  r <- f$sampler_diag$rhat_max
  cat(sprintf("seed %-6d rhat_max=%.3f offmode=%d\n", sd_i, r, f$sampler_diag$n_offmode))
  if (r > best_rhat) { best_rhat <- r; best <- f$draws; best_seed <- sd_i }
  if (r > 1.15) break                      # a genuinely bad fit -- stop and dissect it
}
cat(sprintf("\nWORST fit: seed %d, rhat_max = %.3f\n", best_seed, best_rhat))
saveRDS(list(draws = best, seed = best_seed, rhat = best_rhat), file.path(SP, "ex1_worst_fit.rds"))
analyse(best, sprintf("unit 42, seed %d (rhat_max = %.3f)", best_seed, best_rhat))
cat("\nRead: ess_tail/ess_bulk << 1, qratio >> 2.91, or a negative cor(qratio, ess_bulk) would show the\n")
cat("tails are what is slow. If qratio ~ 2.9 and ess_tail ~ ess_bulk, the marginals are near-normal and\n")
cat("the slow mixing is a ridge/correlation problem instead, not a tail problem.\n")
