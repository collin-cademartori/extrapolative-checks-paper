## Does the ||Lambda|| anchor in the error scale cause stat_strong's elevated rhat?
## A/B on the SAME datasets, identical except one line:
##   WITH  : error_precisions = inv(omega_sq * square(tau * sigma * ||Lambda[m,:]||))   [current]
##   WITHOUT: error_precisions = inv(omega_sq * square(tau * sigma))                    [original form]
## The anchor was added so tau is a true noise-to-signal RATIO and stays inside its prior; the
## suspicion is that it also introduces a singularity as a row norm -> 0 (error scale -> 0,
## precision -> inf), which would hit stat_strong hardest because its tau prior is centred at 0.05.
## Reports rhat AND which parameter carries it, divergences, tau (+ prior tail), and the row norms --
## in particular min over rows n > K_latent, which the inverse-gamma zero-avoidance does NOT protect.
library(posterior)
source("../sample_model.r"); source("../pathfinder_init.r")
WITH <- cmdstan_model("../ife_named_cartesian.stan")
WITHOUT <- cmdstan_model("../ife_named_cartesian_nolam.stan")
rms <- function(x) sqrt(mean(x^2))

GEN_K <- 4L; REPS <- 30L; NPT <- 5L
set.seed(40318)                                    # same stream as ex1_sim_study
study_units <- sample.int(2 * REPS, size = REPS, replace = FALSE)
pp_seed <- sample.int(.Machine$integer.max, 1)
td <- sample_model(overall_scales = rep(1, 8), err_scale = 2, alpha_diag = 10, autocor_a = 8,
  autocor_b = 2, nonstationary = TRUE, num_treated = 0, type = "prior_pred", K_latent = GEN_K,
  iter = 2 * REPS, seed = pp_seed)                 # err_scale = 2, matching the current study

tail_p <- function(t, m, s) (1 - pnorm(t, m, s)) / (1 - pnorm(0, m, s))

fit_arm <- function(ys, variant) {
  assign("ife_mod", if (variant == "with") WITH else WITHOUT, envir = globalenv())
  T_ <- nrow(ys); N_ <- ncol(ys)
  f <- sample_model(N_units = N_, T_times = T_, K_latent = GEN_K,
    overall_scales = 2 * apply(ys, 2, rms),        # stat_strong config from ex1_sim_study
    err_scale = 0, err_scale_mean = 0.05, err_scale_sd = 0.05, data = ys,
    autocor_a = 97, autocor_b = 3, nonstationary = FALSE, num_treated = NPT,
    fit_scales = FALSE, alpha_diag = 10, pathfinder_init = TRUE, type = "posterior", quiet = TRUE,
    ad = 0.8, iter = 1000, iter_warm = 500, n_chains = 4, seed = 42, parallel_chains = 4,
    return_draws = c("tau", "Lambda", "scale_est", "sigma"))
  d <- f$draws
  # which parameter carries the worst rhat
  vars <- grep("^(Lambda|tau|scale_est|sigma)\\[|^tau$", dimnames(d)$variable, value = TRUE)
  ss <- summarise_draws(subset_draws(d, variable = vars), "rhat")
  worst <- ss$variable[which.max(ss$rhat)]
  # per-draw row norms; separate rows <= K (inv-gamma protected) from rows > K (NOT protected)
  rn <- sapply(1:N_, function(n) {
    ks <- 1:min(GEN_K, n)
    m <- sapply(ks, function(k) as.numeric(extract_variable_matrix(d, sprintf("Lambda[%d,%d]", n, k))))
    sqrt(rowSums(as.matrix(m)^2))
  })
  tau <- median(as.numeric(extract_variable_matrix(d, "tau")))
  sd_ <- f$sampler_diag
  data.frame(variant = variant, rhat = round(sd_$rhat_max, 2), worst = worst,
    div = sd_$n_div, tree = sd_$n_tree, tau = round(tau, 3), tail = round(tail_p(tau, 0.05, 0.05), 3),
    min_rn_all = round(min(rn), 3),
    min_rn_prot = round(min(rn[, 1:GEN_K]), 3),          # rows 1..K: inv-gamma protected
    min_rn_unprot = round(min(rn[, (GEN_K + 1):N_]), 3), # rows K+1..M: NO zero-avoidance
    q01_rn = round(quantile(rn, 0.01), 3), row.names = NULL)
}

cat("\nex1 stat_strong | error-scale anchor A/B | 2x RMS, normal(0.05,0.05) tau, ad=0.8, 1000/500\n")
cat("WITH = tau*sigma*||Lambda|| (current)   WITHOUT = tau*sigma (original form)\n")
cat(sprintf("K_latent=%d of %d units, so rows %d-8 have NO inv-gamma zero-avoidance.\n\n", GEN_K, 8, GEN_K + 1))
out <- list()
for (ds in study_units[1:5]) {
  ys <- td$ys[ds, , ]
  for (v in c("with", "without")) {
    r <- tryCatch(fit_arm(ys, v), error = function(e) { cat(sprintf("  !! ds%d %s FAILED\n", ds, v)); NULL })
    if (is.null(r)) next
    r$ds <- ds; out[[length(out) + 1]] <- r
    cat(sprintf("ds%-3d %-8s rhat=%-5.2f (%-14s) div=%-4d tree=%-4d tau=%.3f tail=%.3f | min||L|| all=%.3f prot=%.3f UNPROT=%.3f\n",
      ds, v, r$rhat, r$worst, r$div, r$tree, r$tau, r$tail, r$min_rn_all, r$min_rn_prot, r$min_rn_unprot))
  }
}
res <- do.call(rbind, out)
saveRDS(res, "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad/ex1_lamnorm_ab.rds")
cat("\n=== by variant ===\n")
for (v in c("with", "without")) { s <- res[res$variant == v, ]
  cat(sprintf("%-8s: rhat med=%.2f max=%.2f | n(rhat>1.05)=%d/%d | div=%d tree=%d | tau=%.3f tail=%.3f | min||L|| unprot=%.3f\n",
    v, median(s$rhat), max(s$rhat), sum(s$rhat > 1.05), nrow(s), sum(s$div), sum(s$tree),
    mean(s$tau), mean(s$tail), min(s$min_rn_unprot))) }
cat("\nRead: if WITHOUT clears the rhat while WITH does not, the ||Lambda|| anchor is the cause -- and the\n")
cat("cost is visible in tau (WITHOUT should push tau up, out of its prior). If the bad fits also show a\n")
cat("small UNPROT row norm, the singularity is the mechanism and the inv-gamma cannot reach it.\n")
