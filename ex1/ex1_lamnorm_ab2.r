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
mk_td <- function(es) sample_model(overall_scales = rep(1, 8), err_scale = es, alpha_diag = 10,
  autocor_a = 8, autocor_b = 2, nonstationary = TRUE, num_treated = 0, type = "prior_pred",
  K_latent = GEN_K, iter = 2 * REPS, seed = pp_seed)
TD <- list(`2` = mk_td(2), `3` = mk_td(3))   # 3 = the overnight run's DGP, where the failure appeared

tail_p <- function(t, m, s) (1 - pnorm(t, m, s)) / (1 - pnorm(0, m, s))

fit_arm <- function(ys, variant, arm) {
  assign("ife_mod", if (variant == "with") WITH else WITHOUT, envir = globalenv())
  T_ <- nrow(ys); N_ <- ncol(ys)
  f <- sample_model(N_units = N_, T_times = T_, K_latent = GEN_K,
    overall_scales = 2 * apply(ys, 2, rms),        # stat_strong config from ex1_sim_study
    err_scale = 0, err_scale_mean = if (arm == "stat_strong") 0.05 else 0.1,
    err_scale_sd = if (arm == "stat_strong") 0.05 else 0.1, data = ys,
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
  data.frame(arm = arm, variant = variant, rhat = round(sd_$rhat_max, 2), worst = worst,
    div = sd_$n_div, tree = sd_$n_tree, tau = round(tau, 3), tail = round(tail_p(tau, if (arm == "stat_strong") 0.05 else 0.1, if (arm == "stat_strong") 0.05 else 0.1), 3),
    min_rn_all = round(min(rn), 3),
    min_rn_prot = round(min(rn[, 1:GEN_K]), 3),          # rows 1..K: inv-gamma protected
    min_rn_unprot = round(min(rn[, (GEN_K + 1):N_]), 3), # rows K+1..M: NO zero-avoidance
    q01_rn = round(quantile(rn, 0.01), 3), row.names = NULL)
}

cat("\nex1 stat_strong | error-scale anchor A/B | 2x RMS, normal(0.05,0.05) tau, ad=0.8, 1000/500\n")
cat("WITH = tau*sigma*||Lambda|| (current)   WITHOUT = tau*sigma (original form)\n")
cat(sprintf("K_latent=%d of %d units, so rows %d-8 have NO inv-gamma zero-avoidance.\n\n", GEN_K, 8, GEN_K + 1))
out <- list()
for (es in c("2", "3")) for (ds in study_units[1:4]) {
  ys <- TD[[es]]$ys[ds, , ]
  for (arm in c("stat_strong", "stat_weak")) for (v in c("with", "without")) {
    r <- tryCatch(fit_arm(ys, v, arm), error = function(e) { cat(sprintf("  !! es%s ds%d %s %s FAILED\n", es, ds, arm, v)); NULL })
    if (is.null(r)) next
    r$ds <- ds; r$es <- es; out[[length(out) + 1]] <- r
    cat(sprintf("es%s ds%-3d %-11s %-8s rhat=%-5.2f (%-13s) div=%-4d tau=%.3f tail=%.3f min||L||=%.3f\n",
      es, ds, arm, v, r$rhat, r$worst, r$div, r$tau, r$tail, r$min_rn_all))
  }
}
res <- do.call(rbind, out)
saveRDS(res, "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad/ex1_lamnorm_ab2.rds")
cat("\n=== by variant ===\n")
for (es in c("2","3")) for (arm in c("stat_strong","stat_weak")) for (v in c("with", "without")) {
  s <- res[res$variant == v & res$arm == arm & res$es == es, ]; if (!nrow(s)) next
  cat(sprintf("DGP err=%s %-11s %-8s: rhat med=%.2f max=%.2f n(>1.05)=%d/%d | div=%d | tau=%.3f tail=%.3f\n",
    es, arm, v, median(s$rhat), max(s$rhat), sum(s$rhat > 1.05), nrow(s), sum(s$div), mean(s$tau), mean(s$tail))) }
cat("\nRead: if WITHOUT clears the rhat while WITH does not, the ||Lambda|| anchor is the cause -- and the\n")
cat("cost is visible in tau (WITHOUT should push tau up, out of its prior). If the bad fits also show a\n")
cat("small UNPROT row norm, the singularity is the mechanism and the inv-gamma cannot reach it.\n")
