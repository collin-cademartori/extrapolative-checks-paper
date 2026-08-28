## This file runs the simulation study for the nonstationary example, fitting the
## nonstationary model and two stationary models (with weak and strong priors on the
## iid error scale) to samples from the prior predictive distribution of the
## nonstationary model.

library(foreach)
library(doParallel)
library(doRNG)
library(purrr)

# (Apple Silicon has no hyperthreading, so base the worker count on physical cores and do not
# oversubscribe.)
# Command line: Rscript ex1_sim_study.r [n_cores] [mode] [reps]
#   n_cores  worker count (default: half the physical cores, less one)
#   mode     "full" (default) or "fast"
#   reps     overrides the mode's default rep count
#
# FAST MODE is for iterating on the model specification -- e.g. searching for a workable trio of tau
# priors -- where the question is whether the arms separate, not whether any single fit has
# converged. It shortens the chains, disables escalation entirely (a zero-length ladder, so
# fit_with_escalation breaks after the first fit), and cuts the rep count. Diagnostics are still
# recorded per fit, so a fast run will show elevated rhat_M; that is expected and is NOT evidence
# about the specification. Never report fast-mode numbers.
.args <- commandArgs(trailingOnly = TRUE)
requested_cores <- suppressWarnings(as.integer(.args[1]))
STUDY_MODE <- if (length(.args) >= 2 && tolower(.args[2]) %in% c("fast", "f")) "fast" else "full"
.reps_arg <- suppressWarnings(as.integer(.args[3]))
n_cores <- if (!is.na(requested_cores) && requested_cores >= 1) {
  requested_cores
} else {
  max(1, round(detectCores() / 2) - 1)
}
# Keep CmdStan's CSVs and R's temp staging OFF /tmp. On the study machine /tmp is a TMPFS, so
# anything written there is held in RAM against a cap of roughly half of physical memory. Two long
# runs died on this: one as an OOM kill, and one -- after the $draws() fix cut R's own heap use --
# as tmpfs hitting its size cap, which surfaced as data.table's "disk is full in the temporary
# directory" and then "No chains finished successfully" once CmdStan could not write at all.
#
# Two separate consumers, and output_dir only covers the first:
#   CMDSTAN_OUTPUT_DIR  the per-chain sampling CSVs, live for the whole duration of a fit
#                       (~264 MB at the round-2 rung, ~396 MB at round 3, times the worker count).
#   TMPDIR              cmdstanr reads draws with data.table::fread(cmd = "grep -v '^#' ..."), which
#                       stages the command's output through a file in tempdir() before parsing it.
#                       That is set from TMPDIR at R STARTUP, so it must be exported BEFORE
#                       makeCluster() -- the PSOCK workers inherit the environment, and setting it
#                       afterwards would not move their tempdir().
# Override either by exporting it before launching the script.
scratch_root <- Sys.getenv("STUDY_SCRATCH", file.path(getwd(), ".scratch"))
dir.create(scratch_root, showWarnings = FALSE, recursive = TRUE)
if (!nzchar(Sys.getenv("CMDSTAN_OUTPUT_DIR")))
  Sys.setenv(CMDSTAN_OUTPUT_DIR = file.path(scratch_root, "cmdstan"))
if (!nzchar(Sys.getenv("STUDY_KEEP_TMPDIR"))) {
  tmp_dir <- file.path(scratch_root, "rtmp")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir)
}
cat(sprintf("  scratch: CMDSTAN_OUTPUT_DIR=%s  TMPDIR=%s\n",
  Sys.getenv("CMDSTAN_OUTPUT_DIR"), Sys.getenv("TMPDIR")))

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

# Convergence escalation ladder for ex1. The mechanism lives in sample_model.r
# (escalation_ladder / fit_with_escalation, with the argument that it is adaptive computation rather
# than selection, and why rhat_M is the criterion); the rungs below are ex1-specific and each is set
# by measurement from the first overnight 200-rep run.
#
# ITERATIONS -- round 2 jumps straight to 8000 rather than doubling. That run's ladder was
# 2000 -> 4000 -> 8000, and the intermediate 4000-iteration rung improved rhat_M in only 12 of 21
# trajectories while making it WORSE in 6, twice badly (57 stat_weak 1.011 -> 1.065, 314 stat_weak
# 1.013 -> 1.052). That is not sampler noise: each round draws a FRESH seed, so a marginal case
# sitting near the 1.01 threshold gets re-rolled rather than continued. Dropping the marginal rung
# removes the re-roll and spends the compute where it demonstrably paid -- every R2 -> R3 step in
# that run improved rhat_M (8 of 8, median 1.017 -> 1.005).
#
# WARMUP -- escalates in step. The base configuration warms up for 500 and samples for 2000; leaving
# warmup at 500 while sampling ran to 20000 would draw tens of thousands of iterations from a metric
# adapted on a quarter of a short run. Warming up longer at FIXED sampling length was measured on six
# of the problem datasets, and in all three stat_strong cases -- the worst arm -- rhat_M improved
# monotonically from warmup alone:
#   265 stat_strong  1.092 -> 1.038 -> 1.038   (warm 500 -> 2000 -> 4000)
#   278 stat_strong  1.017 -> 1.010 -> 1.004
#   120 stat_strong  1.014 ->   .   -> 1.009
# with ess_delta1 moving in step (120: 454 -> 734, 278: 1741 -> 1986, 187 stat_weak: 2152 -> 2678).
# The stat_weak cases were flat, but were already clean. Do NOT expect it to move E-BFMI: the same
# test refuted that (265 stat_strong 0.27 -> 0.30 -> 0.27, 278 stat_strong 0.47 -> 0.51 -> 0.50,
# 120 stat_strong 0.33 -> 0.33, 187 stat_weak 0.49 -> 0.44 -> 0.54, 203 stat_weak 0.40 -> 0.43 ->
# 0.39). See the note on ebfmi_min below.
#
# ADAPT_DELTA FLOOR -- every escalated round runs at >= 0.95 whatever triggered it. Raising
# adapt_delta to 0.95 drove divergences to exactly ZERO in 8 of the 8 fits where it was tried
# (initial rates 0.12%-0.63%), whereas fits that bought iterations at adapt_delta = 0.8 saw
# divergences GROW with chain length (278 stat_weak 2 -> 8 -> 166, i.e. 0.03% -> 0.69%; 183 stat_weak
# 0 -> 0 -> 34). One caveat it also absorbs: 265 stat_strong picked up 8 divergences at warm = 4000
# having had none at 500 or 2000, so longer warmup is not unconditionally benign.
#
# ESS -- retained as a cheap tripwire, but it never bound in that run: the minimum ess_delta1 over
# all 127 logged fits was 729, against the threshold of 400. Every escalation was triggered by
# rhat_M (13), the divergence rate (4), or both (3).
#
# Round 3 is 12000, not the 20000 first tried. Two reasons, both from the run that was OOM-killed:
#   * Memory. Even after the $draws() fix in sample_model.r, Y_latent and Y_pred are legitimately
#     materialized per fit, so cost still scales with iterations: roughly 200 MB per fit at 20000
#     against 120 MB at 12000, in every worker that reaches round 3, concurrently.
#   * Diminishing returns. Round 2 at 8000 resolved 14 of the 17 escalated fits, and the three still
#     due a round 3 sat at rhat_M 1.011, 1.016 and 1.058 -- the first two marginal. Nothing in that
#     run suggests 20000 rather than 12000 is what separates them; 137 stat_weak (1.058) mixed
#     slowly with treedepth saturating 27 times, which is a geometry problem that more iterations
#     do not fix.
EX1_LADDER <- if (STUDY_MODE == "fast") escalation_ladder(integer(0), integer(0)) else
  escalation_ladder(iter = c(8000L, 12000L), warm = c(2000L, 3000L))
ESCALATE_MAX <- EX1_LADDER$max_rounds   # seeds are drawn one per fit per round

# Base sampling length, per mode. All three arms share these.
EX1_ITER <- if (STUDY_MODE == "fast") 500L else 2000L
EX1_WARM <- if (STUDY_MODE == "fast") 500L else 500L

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
  stat_scale_multiple <- 1.5
  overall_scales_stat <- stat_scale_multiple * apply(test_ys, 2, function(y) sqrt(mean(y^2)))
  # For the nonstationary fit, sigma scales the *differenced* series (the model fits
  # on first-differences), so estimate its scale from sd(diff(y)) -- using sd(y) would
  # be the wrong, inflating scale for integrated data.
  overall_scales_nonstat <- apply(test_ys, 2, function(y) sd(diff(y)))

  # Draw every Stan seed up front, before any sample_model() call: cmdstanr's $sample() advances R's
  # global RNG, so a seed drawn after a fit would not be reproducible. Invariant: never derive a seed
  # after a fit. One seed per model PER ESCALATION ROUND, so a refit is reproducible too.
  fit_seeds <- matrix(sample.int(.Machine$integer.max, 3L * ESCALATE_MAX), nrow = 3L)

  fits <- list()

  fits$nonstat <- fit_with_escalation(
    list(
      alpha_diag = 20,
      N_units = N_units, T_times = T_times,
      overall_scales = overall_scales_nonstat, err_scale = 0,
      err_scale_mean = 2,
      err_scale_sd = 2,
      data = test_ys,
      autocor_a = 8, autocor_b = 2,
      nonstationary = TRUE, num_treated = 5,
      fit_scales = FALSE,
      type = "posterior", K_latent = K_latent, ad = 0.8,
      # iter = 2000, not 500. At 500 the nonstationary fit degrades sharply and starts escalating:
      # across the two long runs, rhat_M median went 1.002 -> 1.008 and its maximum 1.003 -> 1.027,
      # ess_delta1 median fell 4919 -> 1285, and the share of fits over rhat_M = 1.01 went from
      # 0 of 21 to 13 of 71 (18%). Those escalations are largely INVISIBLE in progress.log -- a
      # round-2 fit is only logged if it trips a threshold, and nonstat's rhat_loadings is small
      # enough that a clean refit leaves no trace -- so the cost showed up as memory pressure rather
      # than as log lines. Paying 4x on the base fit is cheaper than refitting a fifth of them at
      # 8000, and it restores the arm to the near-pristine behaviour it had at 2000.
      iter = EX1_ITER, iter_warm = EX1_WARM,
      n_chains = 3, pathfinder_init = TRUE
    ),
    seeds = fit_seeds[1, ], label = sprintf("unit %d nonstat", i), progress_log = progress_log, ladder = EX1_LADDER
  )
  fits$nonstat$name <- "nonstat"

  fits$stat_weak <- fit_with_escalation(
    list(
      alpha_diag = 20,
      N_units = N_units, T_times = T_times,
      overall_scales = overall_scales_stat, err_scale = 0,
      err_scale_mean = 0.1,
      err_scale_sd = 0.1,
      data = test_ys,
      autocor_a = 97, autocor_b = 3,
      nonstationary = FALSE, num_treated = 5,
      fit_scales = FALSE,
      type = "posterior", K_latent = K_latent + 1, ad = 0.8,
      iter = EX1_ITER, iter_warm = EX1_WARM,
      n_chains = 3, pathfinder_init = TRUE
    ),
    seeds = fit_seeds[2, ], label = sprintf("unit %d stat_weak", i), progress_log = progress_log, ladder = EX1_LADDER
  )
  fits$stat_weak$name <- "stat_weak"

  fits$stat_strong <- fit_with_escalation(
    list(
      alpha_diag = 20,
      N_units = N_units, T_times = T_times,
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
      type = "posterior", K_latent = K_latent + 1, ad = 0.8,
      iter = EX1_ITER, iter_warm = EX1_WARM,
      n_chains = 3, pathfinder_init = TRUE
    ),
    seeds = fit_seeds[3, ], label = sprintf("unit %d stat_strong", i), progress_log = progress_log, ladder = EX1_LADDER
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

      # Convergence diagnostics recorded for EVERY fit. progress.log only keeps fits that trip a
      # threshold, so on its own it cannot show that rhat_M stays reasonable across the run -- the
      # clean majority would be invisible. rhat_M is the one that certifies delta (see sample_model).
      sdg <- pfit$sampler_diag
      res$rhat_max <- sdg$rhat_max
      res$rhat_M <- sdg$rhat_M
      res$rhat_estimands <- sdg$rhat_estimands
      res$rhat_loadings <- sdg$rhat_loadings
      res$rhat_cor_sq <- sdg$rhat_cor_sq
      res$n_div <- sdg$n_div
      res$lp_gap_max <- sdg$lp_gap_max
      res$ess_delta1 <- sdg$ess_delta1
      # How much computation this fit needed: 1 = met the criterion on the first attempt.
      res$n_rounds <- pfit$n_rounds
      res$final_iter <- pfit$final_iter
      res$final_warm <- pfit$final_warm
      res$final_ad <- pfit$final_ad
      # E-BFMI is recorded but is NOT part of the escalation criterion, because no amount of
      # computation moves it: it is a function of the adapted metric and step size, so it is inert to
      # sampling length by construction, and was measured to be nearly inert to warmup as well. It is
      # also concentrated in the stationary arms (median 0.41, 10 of 76 fits below 0.3, versus median
      # 0.78 and none below 0.3 for nonstat) yet essentially unrelated to the mixing of the estimands:
      # across the stationary round-1 fits, cor(ebfmi, rhat_M) = -0.12 and cor(ebfmi, ess_delta1) =
      # +0.07, against -0.34 with lp_gap_max and -0.31 with the divergence rate. So it tracks the
      # energy geometry, not delta. Recording it per fit is what makes that checkable rather than
      # asserted: the summary can compare delta error and noise_abs_tr between the low-E-BFMI fits and
      # the rest, exactly as for rhat_M.
      res$ebfmi_min <- sdg$ebfmi_min
      res$n_tree <- sdg$n_tree

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
    "pathfinder_inits", "draw_to_init", "PF_PARAM_BASES",
    "fit_with_escalation", "escalation_ladder", "ESCALATE_MAX", "EX1_LADDER",
    "EX1_ITER", "EX1_WARM")
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
  t0 <- Sys.time()

  # Per-rep checkpoints. A single failing task used to abort the whole study via %dorng% and take
  # every completed rep with it -- which is how ~1 hour of a 200-rep run was lost to a disk-quota
  # error on the temp directory. Each rep now writes its own row as soon as it finishes, and a rerun
  # picks up the rows that already exist instead of recomputing them. Two consequences worth knowing:
  #   * A rep that errors records a `failed = TRUE` row and the study continues; the error text is
  #     kept in the row and echoed to progress.log.
  #   * The checkpoint directory is keyed to the study seed and rep count, so changing either starts
  #     a fresh set rather than silently reusing rows from a different configuration. Delete the
  #     directory by hand to force a full recompute under the same configuration.
  # STUDY_MODE is in the directory name so a fast run's checkpoints can never be resumed into a
  # full run -- the two hold results from different sampler configurations at the same rep index.
  ckpt_dir <- file.path(getwd(),
    sprintf("ckpt_ns_%s_seed%d_reps%d", STUDY_MODE, seed, reps))
  dir.create(ckpt_dir, showWarnings = FALSE)
  n_resume <- length(list.files(ckpt_dir, pattern = "^rep_.*\\.rds$"))
  if (n_resume > 0) {
    cat(sprintf("  resuming: %d of %d reps already checkpointed in %s\n\n",
      n_resume, reps, basename(ckpt_dir)))
  } else {
    cat("", file = progress_log) # truncate: fresh per-worker progress log per run
  }

  # Single (non-nested) foreach, so %dorng% gives each task a reproducible RNG
  # stream invariant to worker count. Invariant: keep this a single, non-nested loop.
  study_res <-
    foreach(
      s = study_units, iter = seq(reps),
      # bind_rows rather than rbind: a failed rep contributes a short row, and rbind would error on
      # the column mismatch instead of recording the failure.
      .combine = function(...) dplyr::bind_rows(...),
      .export = exp_vars, .packages = exp_packages,
      .options.RNG = seed
    ) %dorng% {
      # Rep index AND unit in the filename, so a checkpoint can never be reused for a different
      # dataset even if study_units were to change under a reused directory.
      ckpt_file <- file.path(ckpt_dir, sprintf("rep_%04d_unit_%04d.rds", iter, s))
      if (file.exists(ckpt_file)) {
        readRDS(ckpt_file)
      } else {
        unit_res <- tryCatch(
          as.data.frame(run_sim_stat(test_data, s, K_latent, post_check, progress_log = progress_log)),
          error = function(e) {
            cat(sprintf("[%s] iteration %d (unit %d) FAILED: %s\n",
              format(Sys.time(), "%H:%M"), iter, s, conditionMessage(e)),
              file = progress_log, append = TRUE)
            data.frame(error = conditionMessage(e))
          }
        )
        unit_res$rep <- iter
        unit_res$unit <- s
        if (is.null(unit_res$error)) unit_res$error <- NA_character_
        unit_res$failed <- !is.na(unit_res$error)
        saveRDS(unit_res, ckpt_file)
        worker_progress(sprintf("iteration %d (unit %d)", iter, s), logfile = progress_log)
        unit_res
      }
    }

  n_failed <- sum(study_res$failed, na.rm = TRUE)
  cat(sprintf(
    "--- study complete: %d tasks in %.1f min%s ---\n",
    reps, as.numeric(difftime(Sys.time(), t0, units = "mins")),
    if (n_failed > 0) sprintf("; %d FAILED (see the `error` column and progress.log)", n_failed) else ""
  ))

  return(study_res)
}

study_reps <- if (!is.na(.reps_arg)) .reps_arg else if (STUDY_MODE == "fast") 25L else 200L
cat(sprintf("\n=== mode: %s | reps: %d | iter/warm: %d/%d | escalation: %s ===\n",
  STUDY_MODE, study_reps, EX1_ITER, EX1_WARM,
  if (EX1_LADDER$max_rounds > 1) sprintf("up to %d rounds", EX1_LADDER$max_rounds) else "OFF"))
if (STUDY_MODE == "fast")
  cat("    FAST MODE -- for specification search only. Elevated rhat_M is expected here and says\n",
      "   nothing about the specification. Do not report these numbers.\n", sep = "")
sim_study_stat <- run_sim_study_stat(K_latent = 4, study_reps, seed = 40318)

stopCluster(cl)

# Save the raw study results; numeric summaries and plots are produced by
# ex1_sim_study_summary.r (run it to view the results).
# Mode-keyed output, so a quick fast-mode check cannot silently destroy the full study's results.
out_file <- if (STUDY_MODE == "fast") "sim_study_ns_fast.RData" else "sim_study_ns.RData"
save(sim_study_stat, file = out_file)
cat(sprintf("Results saved to %s -- run `Rscript ex1_sim_study_summary.r %s` to summarize.\n",
  out_file, out_file))
