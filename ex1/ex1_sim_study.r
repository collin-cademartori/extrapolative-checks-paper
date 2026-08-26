## This file runs the simulation study for the nonstationary example, fitting the
## nonstationary model and two stationary models (with weak and strong priors on the
## iid error scale) to samples from the prior predictive distribution of the
## nonstationary model.

library(foreach)
library(doParallel)
library(doRNG)
library(purrr)

# Worker count: pass as the first CLI arg (e.g. `Rscript ex1_sim_study.r 4`),
# otherwise a conservative default. (Apple Silicon has no hyperthreading, so base
# this on physical cores and do not oversubscribe.)
requested_cores <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)[1]))
n_cores <- if (!is.na(requested_cores) && requested_cores >= 1) {
  requested_cores
} else {
  max(1, round(detectCores() / 2) - 1)
}
cl <- makeCluster(n_cores, outfile = "")
registerDoParallel(cl)

# PSOCK workers may start in a different working directory than the master; sync
# them so worker-side relative paths (progress.log, ggsave to ../figs) resolve the
# same as here.
invisible(clusterCall(cl, setwd, getwd()))

# Pre-attach the workers' packages quietly, so their startup banners don't clutter
# the console (outfile = "" surfaces all worker output).
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(forcats)
  library(dplyr)
  library(ggplot2)
  library(doRNG)
  library(purrr)
})))

source("../sample_model.r")
source("../pathfinder_init.r")
source("../plotting.r")

# Uses sample_model.r's default model, ../ife_named.stan, whose error scale is tau*sigma and does NOT
# involve ||Lambda||. That matters beyond convention: it is what makes the likelihood depend on
# (Lambda, Phi) only through the product M = Lambda_Phi, giving
#     delta  _||_  (Lambda, Phi)  |  M, sigma, tau,
# so the rhat_M diagnostic reported by sample_model certifies delta. Anchoring the error scale to
# sigma*||Lambda|| (the Cartesian variant) breaks that conditional independence, and also measured a
# WORSE tau: at 2x RMS it pushed tau up (stat_weak prior tail 0.22 vs 0.42, stat_strong 0.017 vs 0.093).
stopifnot("Lambda" %in% PF_PARAM_BASES)

# Per-worker progress: each worker appends a running count of completed tasks to a
# shared log. cat(append = TRUE) flushes every write, so it shows up live via
# `tail -f progress.log` -- unlike message()/outfile, which block-buffers a
# redirected stream (nothing appears until the worker exits).
worker_progress <- function(label, logfile = "progress.log") {
  n <- get0(".worker_done", envir = globalenv(), ifnotfound = 0L) + 1L
  assign(".worker_done", n, envir = globalenv())
  cat(sprintf(
    "[%s] [worker %d] %3d done | %s\n",
    format(Sys.time(), "%H:%M"), Sys.getpid(), n, label
  ), file = logfile, append = TRUE)
}

run_sim_stat <- function(test_data, i, K_latent, post_check = FALSE, progress_log = NULL) {
  test_ys <- test_data$ys[i, , ]
  N_units <- ncol(test_ys)
  T_times <- nrow(test_ys)
  # Length of the treatment window; the three fits below all pass num_treated = this.
  num_treated_ex1 <- 5
  # Prior scale for the stationary fits: a multiple of RMS(y), not sd(y).
  #
  # RMS rather than sd because these models have no intercept and no factor means, so the level of
  # each series must be produced by the factors themselves; the second moment about zero is what the
  # factors have to reproduce.
  #
  # The multiple corrects for the realised amplitude of a near-unit-root factor over a short window.
  # The factors have unit LONG-RUN variance, but at rho ~ 0.97 over T = 20 a realised path's sample sd
  # averages only ~0.38 (closed form: E[s^2] = [(T-1) - (2/T)*sum_k (T-k) rho^k] / (T-1)), because the
  # sample mean absorbs the low-frequency wandering. Without the correction the factors cannot generate
  # the observed excursion, the shortfall is booked as error, and tau is dragged far above its
  # N(0.1, 0.1) prior -- at which point stat_weak mimics nonstat by inflating noise rather than by
  # fitting structure, which is not the phenomenon the example is meant to show.
  #
  # Measured on these datasets: 1.0x -> tau 0.35 (0.7% prior tail, a clear contradiction);
  # 1.5x -> 0.275 (3.7-7.1%); 2.0x -> 0.225 (12.8%, inside the bulk); 2.5x -> 0.186 with no further
  # change in the overfitting or delta metrics. The prior predictive is consistent with this range
  # (it centres on the observed sd at ~1.7x) but is too diffuse to discriminate on its own -- a
  # near-unit-root process has ~3.3x spread in realised sd even at FIXED rho, so the prior-posterior
  # consistency of tau is what selects the multiple.
  stat_scale_multiple <- 2
  overall_scales_stat <- stat_scale_multiple * apply(test_ys, 2, function(y) sqrt(mean(y^2)))
  # For the nonstationary fit, sigma scales the *differenced* series (the model fits
  # on first-differences), so estimate its scale from sd(diff(y)) -- using sd(y) would
  # be the wrong, inflating scale for integrated data.
  overall_scales_nonstat <- apply(test_ys, 2, function(y) sd(diff(y)))

  # Draw all three Stan seeds up front, before any sample_model() call: cmdstanr's
  # $sample() advances R's global RNG, so a seed drawn after a fit would not be
  # reproducible. Invariant: never derive a seed after a fit.
  fit_seeds <- sample.int(.Machine$integer.max, 3)

  fits <- list()

  fits$nonstat <- sample_model(
    alpha_diag = 10,
    overall_scales = overall_scales_nonstat, err_scale = 0,
    err_scale_mean = 2,
    err_scale_sd = 2,
    data = test_ys,
    autocor_a = 8, autocor_b = 2,
    nonstationary = TRUE, num_treated = 5,
    fit_scales = FALSE,
    type = "posterior", K_latent = K_latent, ad = 0.8,
    iter = 1000, iter_warm = 500,
    n_chains = 4, seed = fit_seeds[1], pathfinder_init = TRUE,
    log_file = progress_log, log_label = sprintf("unit %d nonstat", i)
  )
  fits$nonstat$name <- "nonstat"

  fits$stat_weak <- sample_model(
    alpha_diag = 10,
    overall_scales = overall_scales_stat, err_scale = 0,
    err_scale_mean = 0.1,
    err_scale_sd = 0.1,
    data = test_ys,
    autocor_a = 97, autocor_b = 3,
    nonstationary = FALSE, num_treated = 5,
    fit_scales = FALSE,
    type = "posterior", K_latent = K_latent, ad = 0.8,
    iter = 1000, iter_warm = 500,
    n_chains = 4, seed = fit_seeds[2], pathfinder_init = TRUE,
    log_file = progress_log, log_label = sprintf("unit %d stat_weak", i)
  )
  fits$stat_weak$name <- "stat_weak"

  fits$stat_strong <- sample_model(
    alpha_diag = 10,
    # Stronger error prior: the same truncated-normal form as stat_weak, with the location and
    # scale halved (0.1 -> 0.05). Previously tau was FIXED at 0.1, which made stat_strong differ
    # from stat_weak in kind (no error-scale uncertainty at all) rather than in degree; estimating
    # tau under a tighter prior isolates the strength of the prior as the only difference.
    overall_scales = overall_scales_stat, err_scale = 0,
    err_scale_mean = 0.05,
    err_scale_sd = 0.05,
    data = test_ys,
    autocor_a = 97, autocor_b = 3,
    nonstationary = FALSE, num_treated = 5,
    fit_scales = FALSE,
    type = "posterior", K_latent = K_latent, ad = 0.8,
    iter = 1000, iter_warm = 500,
    n_chains = 4, seed = fit_seeds[3], pathfinder_init = TRUE,
    log_file = progress_log, log_label = sprintf("unit %d stat_strong", i)
  )
  fits$stat_strong$name <- "stat_strong"

  res <- fits |>
    map(function(pfit) {
      res <- list()

      stat_y_pred <- pfit$y_pred
      pred_inc <- matrix(NA, nrow = T_times, ncol = N_units)
      pred_width <- matrix(NA, nrow = T_times, ncol = N_units)
      for (n in 1:N_units) {
        for (t in 1:T_times) {
          y_bounds <- quantile(stat_y_pred[, t, n], c(0.025, 0.975))
          pred_inc[t, n] <-
            (test_ys[t, n] >= y_bounds[1]) &&
              (test_ys[t, n] <= y_bounds[2])
          pred_width[t, n] <- (y_bounds[2] - y_bounds[1]) / (max(test_ys[, n]) - min(test_ys[, n]))
        }
      }
      res$pred_perc <- mean(pred_inc)
      res$pred_width <- mean(pred_width)

      # tau, the noise-to-signal ratio: the quantity that decides whether stat_weak collapses onto
      # nonstat. Under the Cartesian parametrization the error sd is tau * sigma_data[n]*||Lambda[n,:]||,
      # so tau ~ N / S with S the effective scale -- record it (and how far into its own prior's upper
      # tail it sits) so scale changes can be judged directly instead of inferred from the plots.
      tau_draws <- pfit$err_scale
      res$tau_med <- median(tau_draws)
      res$tau_q95 <- unname(quantile(tau_draws, 0.95))
      # prior upper-tail probability of the posterior median, under this fit's own truncated tau
      # prior. All three now estimate tau, so all three get a tail (stat_strong's used to be NA
      # because its tau was fixed at 0.1).
      tau_prior <- switch(pfit$name,
        nonstat     = c(2,    2),
        stat_weak   = c(0.1,  0.1),
        stat_strong = c(0.05, 0.05)
      )
      res$tau_prior_tail <- (1 - pnorm(res$tau_med, tau_prior[1], tau_prior[2])) /
        (1 - pnorm(0, tau_prior[1], tau_prior[2]))

      res$time_cor_pval <- pfit$time_cor_pval

      absz <- abs(pfit$effect_means / pfit$effect_sds)
      res[paste0("absz_", seq_along(absz))] <- absz

      pmean <- pfit$effect_means
      res[paste0("mean_", seq_along(pmean))] <- pmean

      psds <- pfit$effect_sds
      res[paste0("sd_", seq_along(psds))] <- psds

      pred_mad <- pfit$mean_abs_diffs
      res$pred_mad <- pred_mad

      # Overfitting of the treated unit's pre-treatment window -- the basis the counterfactual is
      # extrapolated from, and the only fit that feeds the delta estimate.
      #   noise_abs_tr = cov(fitted - truth, observed - truth) / var(observed - truth)
      # i.e. the fraction of THIS unit's noise absorbed into its fitted signal: 1 = interpolation,
      # 0 = noise ignored. For a linear smoother this estimates tr(H)/n, the effective degrees of
      # freedom per observation. It cannot be computed in generated quantities because Stan never sees
      # the truth; here the DGP's latent signal is available as test_data$ys_latent.
      #
      # Chosen over mean|fitted - observed| (pred_mad) and over mean|fitted - truth|: across the
      # sigma/rho configurations tried, noise_abs_tr tracked the delta error (absz) at r ~ +0.77
      # (within-arm +0.94), versus +0.35 for scale-normalised mean|fitted - truth| and ~0 for the raw
      # magnitudes. The magnitude says how far the basis is from the truth; only the covariance says
      # whether that error is ALIGNED WITH THE NOISE, and only noise-aligned error corrupts the
      # extrapolation.
      true_ys <- test_data$ys_latent[i, , ]
      fit_means <- apply(pfit$y_means, c(2, 3), mean)
      pre_times <- seq_len(T_times - num_treated_ex1)
      fit_err <- fit_means[pre_times, 1] - true_ys[pre_times, 1]
      unit_noise <- test_ys[pre_times, 1] - true_ys[pre_times, 1]
      res$noise_abs_tr <- as.numeric(cov(fit_err, unit_noise) / var(unit_noise))

      return(res)
    }) |>
    list_flatten()

  pns_means <- apply(fits$nonstat$y_means, c(2, 3), mean)
  p2_means <- apply(fits$stat_weak$y_means, c(2, 3), mean)
  p1_means <- apply(fits$stat_strong$y_means, c(2, 3), mean)

  # Show a single unit's fits per rep
  plot_unit <- 2
  fit_plot <- plot_post_fits_stat(test_ys, pns_means, p2_means, p1_means, unit = plot_unit)
  ggsave(
    fit_plot,
    file = paste0("../figs/sim_stat_figs/post_fit_plot_u", plot_unit, "_", i, ".png"),
    create.dir = TRUE
  )

  if (post_check) {
    check_plot <- plot_data_matrix_post(test_ys, fits$stat_weak$y_pred)
    ggsave(
      check_plot,
      file = paste0("../figs/sim_stat_figs/check_plot_", i, ".pdf"),
      device = "pdf", height = 4, width = 8, create.dir = TRUE
    )
  }

  return(res)
}

run_sim_study_stat <- function(K_latent = 3, reps, seed, post_check = FALSE) {
  # Seed the master-process RNG, and draw the dataset selection and the prior-
  # predictive Stan seed *before* generating test_data, so they are not perturbed
  # by cmdstanr's $sample() (which advances R's RNG).
  set.seed(seed)
  study_units <- sample.int(2 * reps, size = reps, replace = FALSE)
  pp_seed <- sample.int(.Machine$integer.max, 1)

  test_data <- sample_model(
    # err_scale = 2 matches the nonstationary model actually being fit (its tau prior is
    # centred at err_scale_mean = 2), so the DGP is not generating data noisier than any of the
    # fitted models expects.
    overall_scales = rep(1, 8), err_scale = 2,
    autocor_a = 8, autocor_b = 2,
    nonstationary = TRUE, num_treated = 0,
    type = "prior_pred", K_latent = K_latent,
    iter = 2 * reps, seed = pp_seed
  )

  exp_vars <- c("run_sim_stat", "worker_progress", "sample_model", "ife_mod", "plot_post_fits_stat", "plot_data_matrix_post",
    "pathfinder_inits", "draw_to_init", "PF_PARAM_BASES")
  exp_packages <- c("cmdstanr", "posterior", "forcats", "dplyr", "ggplot2")
  cat(sprintf(
    paste0(
      "\n=== Example 1 simulation study (nonstationary) ===\n",
      "  reps    : %d\n",
      "  tasks   : %d  (3 model fits each)\n",
      "  workers : %d   seed: %d\n\n"
    ),
    reps, reps, getDoParWorkers(), seed
  ))
  # Absolute log path, so workers write it where the master expects regardless of
  # their working directory.
  progress_log <- file.path(getwd(), "progress.log")
  cat("", file = progress_log) # truncate: fresh per-worker progress log per run
  t0 <- Sys.time()

  # Single (non-nested) foreach, so %dorng% gives each task a reproducible RNG
  # stream invariant to worker count. Invariant: keep this a single, non-nested loop.
  study_res <-
    foreach(
      s = study_units, iter = seq(reps),
      .combine = "rbind", .export = exp_vars, .packages = exp_packages,
      .options.RNG = seed
    ) %dorng% {
      unit_res <- as.data.frame(run_sim_stat(test_data, s, K_latent, post_check, progress_log = progress_log))
      worker_progress(sprintf("iteration %d (unit %d)", iter, s), logfile = progress_log)
      unit_res
    }

  cat(sprintf(
    "--- study complete: %d tasks in %.1f min ---\n",
    reps, as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ))

  return(study_res)
}

study_reps <- 30
sim_study_stat <- run_sim_study_stat(K_latent = 4, study_reps, seed = 40318)

stopCluster(cl)

# Save the raw study results; numeric summaries and plots are produced by
# ex1_sim_study_summary.r (run it to view the results).
save(sim_study_stat, file = "sim_study_ns.RData")
cat("Results saved to sim_study_ns.RData -- run ex1_sim_study_summary.r to summarize.\n")
