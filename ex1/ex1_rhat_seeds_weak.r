## Reproduce stat_strong's elevated rhat. The earlier A/B used seed = 42 for EVERY fit and saw a max
## rhat of 1.06; ex1_sim_study draws a RANDOM seed per rep (and pathfinder init is seeded from it), so
## if the failures are stochastic minor-mode trapping -- as they were in ex2, where seed 42 trapped one
## dataset while six other seeds ran clean -- a fixed seed would miss them.
## Sweep SEEDS on the same datasets and measure the failure RATE, plus the ex2 mode diagnostics:
##   rhat_all vs rhat_on : rhat over all chains vs only ON-MODE chains (mean lp within 5 of the best).
##                         rhat_on clean while rhat_all is high => minor-mode contamination (a few
##                         chains in a lower basin), not genuine non-convergence.
##   offmode / lp_gap    : how many chains are adrift, and by how much.
## DGP err_scale = 3, matching the overnight run where the failure was observed.
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

SEEDS <- c(42, 7, 101, 500, 1234, 2024)
cat("\nex1 STAT_WEAK | rhat vs SEED | 2x RMS, tau ~ normal(0.1,0.1)+, ad=0.8, 1000/500, 4 chains\n")
cat("DGP err_scale = 3 (the overnight setting). stat_strong (truncated normal) swept clean across 8 seeds; stat_weak is now the worst offender.\n\n")
out <- list()
for (ds in study_units[1:5]) {
  ys <- td$ys[ds, , ]; N_ <- ncol(ys); T_ <- nrow(ys)
  for (sd_i in SEEDS) {
    r <- tryCatch({
      f <- sample_model(N_units = N_, T_times = T_, K_latent = GEN_K,
        overall_scales = 2 * apply(ys, 2, rms),
        err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.1, data = ys,
        autocor_a = 97, autocor_b = 3, nonstationary = FALSE, num_treated = NPT,
        fit_scales = FALSE, alpha_diag = 10, pathfinder_init = TRUE, type = "posterior", quiet = TRUE,
        ad = 0.8, iter = 1000, iter_warm = 500, n_chains = 4, seed = sd_i, parallel_chains = 4,
        return_draws = c("lp__", "tau", "Lambda", "delta", "scale_est"))
      d <- f$draws
      lp <- colMeans(extract_variable_matrix(d, "lp__")); gap <- max(lp) - lp; on <- which(gap <= 5)
      vars <- grep("^(Lambda|delta|scale_est)\\[|^tau$", dimnames(d)$variable, value = TRUE)
      ra <- suppressWarnings(max(summarise_draws(subset_draws(d, variable = vars), "rhat")$rhat, na.rm = TRUE))
      ron <- if (length(on) >= 2) suppressWarnings(max(summarise_draws(
        subset_draws(d, variable = vars, chain = on), "rhat")$rhat, na.rm = TRUE)) else NA
      sd_ <- f$sampler_diag
      data.frame(ds = ds, seed = sd_i, rhat_all = round(ra, 2), rhat_on = round(ron, 2),
        n_on = length(on), offmode = sum(gap > 5), lp_gap = round(max(gap), 1),
        div = sd_$n_div, tau = round(median(as.numeric(extract_variable_matrix(d, "tau"))), 3),
        lps = paste(sprintf("%.0f", sort(lp, decreasing = TRUE)), collapse = ","), row.names = NULL)
    }, error = function(e) { cat(sprintf("  !! ds%d seed%d FAILED\n", ds, sd_i)); NULL })
    if (is.null(r)) next
    out[[length(out) + 1]] <- r
    cat(sprintf("ds%-3d seed%-6d rhat_all=%-5.2f rhat_on=%-5s n_on=%d offmode=%d lp_gap=%-5.1f div=%-4d tau=%.3f  lp=[%s]\n",
      ds, sd_i, r$rhat_all, ifelse(is.na(r$rhat_on), "NA", format(r$rhat_on)), r$n_on, r$offmode,
      r$lp_gap, r$div, r$tau, r$lps))
  }
}
res <- do.call(rbind, out)
saveRDS(res, "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad/ex1_rhat_seeds_weak.rds")
cat(sprintf("\n=== failure rate across %d fits ===\n", nrow(res)))
cat(sprintf("  rhat_all > 1.05 : %d/%d      rhat_all > 1.1 : %d/%d      max = %.2f\n",
  sum(res$rhat_all > 1.05), nrow(res), sum(res$rhat_all > 1.1), nrow(res), max(res$rhat_all)))
cat(sprintf("  of the bad fits, how many have CLEAN on-mode rhat (<1.05)? %d  -> minor-mode contamination\n",
  sum(res$rhat_all > 1.05 & res$rhat_on < 1.05, na.rm = TRUE)))
cat(sprintf("  offmode chains: mean %.2f of 4, max %d | lp_gap max %.1f\n",
  mean(res$offmode), max(res$offmode), max(res$lp_gap)))
for (ds in unique(res$ds)) { s <- res[res$ds == ds, ]
  cat(sprintf("  ds%-3d: %d/%d seeds with rhat>1.05 (max %.2f)\n", ds, sum(s$rhat_all > 1.05), nrow(s), max(s$rhat_all))) }
cat("\nRead: if the failure rate is ~50% across SEEDS on a fixed dataset, it is stochastic trapping, and\n")
cat("rhat_on tells us whether it is minor-mode contamination (fixable by chains/init) or real.\n")
