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
# Workers source the shared config themselves rather than receiving each constant through
# .export: adding a constant later would otherwise need remembering to export it, and a worker
# running a different value than the master is exactly the drift ex1_config.r exists to prevent.
invisible(clusterCall(cl, source, "ex1_config.r"))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(forcats)
  library(dplyr)
  library(ggplot2)
  library(doRNG)
  library(purrr)
})))

source("ex1_config.r")
source("../sample_model.r")
source("../pathfinder_init.r")
source("../plotting.r")

# NOTE ON THE MEASUREMENTS BELOW: the escalation notes in this header were recorded when ex1 fitted
# THREE arms -- nonstat plus two stationary arms, stat_weak and stat_strong, differing only in their
# error-scale prior. Only one stationary arm remains, now called `stat`, carrying what stat_weak had.
# The old names are left in those records rather than rewritten, because rewriting them would
# misreport which configuration was actually measured.
#
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
# How many reps get a per-rep fit figure (see the gate in run_sim_stat).
PLOT_REPS <- 25L

EX1_ITER <- if (STUDY_MODE == "fast") 500L else 2000L
EX1_WARM <- if (STUDY_MODE == "fast") 500L else 500L

# Data-driven anchor ordering for the triangular (Cholesky) loadings. Copied verbatim from
# ex2_sim_study.r. Keeps the treated unit (column 1) first, then greedily picks the untreated
# columns most orthogonal to those already chosen so the leading K x K loading block is full rank --
# giving every factor a distinct anchor and avoiding near-zero Cholesky diagonals. This is a
# reparameterization only: the fitted means, treatment effect, and (permutation-invariant) check
# statistics do not depend on the choice.
#
# ex2's companion unpermute_untreated() is NOT ported: it exists to map per-untreated-unit outputs
# (cor_sq, abs_cors_err) back to the original order, and ex1 stores no per-unit fit output.
# Checkpoint compatibility fingerprint. The directory key (mode, seed, rep count) does NOT capture
# the model CONFIGURATION, so a checkpoint set written under a different one is silently resumed and
# mixed in. That has already cost a run: after ex1 dropped from three arms to two, an old set was
# resumed, bind_rows filled NA for the columns the old rows lacked, those rows carried
# failed = FALSE so they survived the filter, and the summary died on
# `quantile(stat_pred_perc, 0.05)` with "missing values and NaN's not allowed". The crash was the
# lucky outcome -- had the column names lined up, two configurations would have been averaged
# together silently.
#
# The fingerprint covers what actually determines the column set and the arms' meaning: the arm
# names and the contents of the shared config. It deliberately does NOT hash the whole study script,
# so editing a comment mid-run does not throw away hours of completed tasks.
ckpt_fingerprint <- function(arms, config_file) {
  paste(c(paste(sort(arms), collapse = ","),
          unname(tools::md5sum(config_file))), collapse = " ")
}

# Refuse to resume a set written under a different configuration, rather than silently mixing it in.
ckpt_check_fingerprint <- function(ckpt_dir, arms, config_file) {
  fp_file <- file.path(ckpt_dir, "FINGERPRINT")
  fp <- ckpt_fingerprint(arms, config_file)
  if (file.exists(fp_file)) {
    old <- readLines(fp_file, warn = FALSE)[1]
    if (!identical(old, fp)) {
      stop("checkpoint directory ", basename(ckpt_dir), " was written under a DIFFERENT model ",
           "configuration (arms or ", basename(config_file), " have changed).\n",
           "  stored:  ", old, "\n  current: ", fp, "\n",
           "Resuming it would mix configurations. Delete the directory to start a fresh set.")
    }
  } else if (length(list.files(ckpt_dir, pattern = "\\.rds$"))) {
    stop("checkpoint directory ", basename(ckpt_dir), " holds results but no FINGERPRINT, so it ",
         "predates this check and its configuration cannot be verified. Delete it to start fresh.")
  } else {
    writeLines(fp, fp_file)
  }
  invisible(fp)
}

anchor_order <- function(y, K) {
  N <- ncol(y)
  yc <- scale(y, center = TRUE, scale = FALSE)
  sel <- 1L
  remaining <- setdiff(seq_len(N), sel)
  while (length(sel) < K && length(remaining) > 0) {
    Q <- qr.Q(qr(yc[, sel, drop = FALSE]))
    resid <- yc[, remaining, drop = FALSE] - Q %*% crossprod(Q, yc[, remaining, drop = FALSE])
    pick <- remaining[which.max(colSums(resid^2))]
    sel <- c(sel, pick)
    remaining <- setdiff(remaining, pick)
  }
  c(sel, remaining)
}

run_sim_stat <- function(test_data, i, K_latent, post_check = FALSE, progress_log = NULL) {
  test_ys <- test_data$ys[i, , ]
  N_units <- ncol(test_ys)
  T_times <- nrow(test_ys)

  # Anchor-order the columns before fitting, exactly as ex2 does. Two ex1-specific points:
  #   * The permutation is computed at K_latent + 1, one more than any arm now fits -- all three
  #     use K_latent. That is deliberate slack, not a leftover: anchor_order's selection is greedy
  #     and therefore prefix-consistent (anchor_order(y, 5)[1:4] equals anchor_order(y, 4)[1:4]),
  #     so asking for one extra leaves the leading K_latent block exactly as it would have been
  #     while keeping the ordering stable if an arm is ever moved back to K_latent + 1. Changing
  #     the argument would change the permutation and hence every fitted dataset, so it is left
  #     as is.
  #   * perm[1] == 1 by construction, so the treated unit stays in column 1. Everything downstream
  #     that indexes the treated unit by position (noise_abs_tr, delta) is therefore unaffected;
  #     the stopifnot makes that dependency explicit rather than implicit.
  perm <- anchor_order(test_ys, K_latent + 1)
  stopifnot(perm[1] == 1)
  fit_ys <- test_ys[, perm]
  # The DGP's latent signal is permuted alongside, so it stays column-aligned with the fit output
  # that noise_abs_tr compares it against.
  true_ys_perm <- test_data$ys_latent[i, , perm]
  # Length of the treatment window. The three fits below pass this, and noise_abs_tr's pre_times
  # is T_times minus this, so the overfitting measurement window and the fitted treatment window
  # cannot drift apart. They used to: the fits hard-coded 5 while pre_times read this variable.
  num_treated_ex1 <- NUM_TREATED
  # ---- scales -------------------------------------------------------------------------------
  # Every scale constant below is a fraction of ONE anchor, RMS(y_n). RMS rather than sd(y) because
  # these models have no intercept and no factor means, so the level of each series must be produced
  # by the factors themselves; the second moment about zero is what they have to reproduce. RMS
  # rather than sd(diff(y)) because sd(diff(y)) is 93% observation error by variance in this DGP
  # (mean sd(diff SIGNAL) = 0.67 against mean sd(diff OBSERVED) = 2.94, the iid level error alone
  # contributing sqrt(2)*eta = 2.83) -- it measures eta, not sigma. RMS is the other way round, 86%
  # signal by variance, which is why it can anchor both arms.
  #
  # One anchor means the constants are directly comparable, to each other and to the truth
  # (values from ex1_config.r; ETA_FRAC_STAT is still provisional):
  #
  #     arm            sigma        eta / RMS(y)
  #     truth (DGP)    RMS/6.8         0.295
  #     nonstat        RMS/7           0.286    estimated, prior centred on the truth
  #     stat           2*RMS           0.200    estimated
  #
  # That table is the mechanism this example demonstrates. The stationary arm is handed an error
  # smaller than the truth, and that single number does both jobs: it is what produces
  # the overfitting, and it is what the prior predictive check on |cor(y, t)| requires, since a large
  # error attenuates a linear trend and only a small one lets a stationary model look nonstationary
  # enough to be plausible. The example is exactly that trade-off -- under a wrong model, a prior
  # that survives its own predictive check is a prior that overfits.
  #
  # The two SIGMA multiples are different KINDS of constant and will break differently:
  #   * SIGMA_MULT_NONSTAT is a units conversion. In the nonstationary branch sigma denotes the
  #     scale of the DIFFERENCED signal (Y_means_0 = sigma * Lambda_Phi, then cumulative_sum), so
  #     an anchor measured on the level has to be converted, and integrating a T = 20 window costs
  #     roughly an order of magnitude. It is the self-consistency fixed point: generate at
  #     sigma = RMS/7 and the prior predictive gives that RMS back (derived to -0.4% by
  #     ex1_derive_scales.r -- run it to see where the number comes from). It moves with T and with
  #     K_latent, and is nearly inert to the error size: the required multiple runs 1/6.5, 1/6.6,
  #     1/7.0, 1/8.0 as err_sd goes 0.5, 1, 2, 4.
  #
  #     This arm shares the DGP's rho prior, Beta(8, 2), so the fixed point does double duty: the
  #     multiple that reproduces the observed RMS also hands the arm the DGP's own sigma and eta,
  #     1.00 and 2.00. That is the reference arm being handed the truth -- deliberate, and it
  #     belongs in the write-up. Were the two priors allowed to diverge the criteria would separate
  #     (at Beta(7, 3) self-consistency gives 1/6.2, about 10% above the truth) and the derivation
  #     would have to say which one it solves.
  #   * SIGMA_MULT_STAT is NOT a conversion either, but for a different reason: under this
  #     misspecification there is no multiple that matches the data's RMS and its SD at once, so it
  #     chooses which to match. The two fixed points are far apart -- sigma = 1.28 * RMS reproduces
  #     the observed RMS, sigma = 2.02 * RMS reproduces the observed SD -- because a near-unit-root
  #     AR(1) realises far less dispersion over a short window than its long-run marginal sd: at
  #     T = 20 the closed form E[s^2] = [(T-1) - (2/T) sum_k (T-k) rho^k] / (T-1) gives E[s]/sigma_lr
  #     = 0.429 at rho = 0.97, and 0.333 averaged over rho ~ Beta(98, 2).
  #
  #     Matching the SD is the conservative choice, and that is why the code uses 2. The SD is what
  #     controls how far the fitted factors can wiggle with the error term, so it is the quantity
  #     that governs the overfitting this arm exists to exhibit. At 2.0x the prior predictive
  #     reproduces the data's SD to 99% (4.13 against 4.18) while exceeding its RMS by 1.57x --
  #     a prior that can reach the truth and then some. At the RMS fixed point it would instead
  #     reach only 60% of the data's SD (2.52 against 4.18), PROHIBITING the realised dispersion the
  #     data actually show. Overshooting a moment is a weak assumption; excluding one is a strong
  #     assumption made silently.
  #
  #     This is a prior-only argument and needs no posterior evidence, but the posterior agrees:
  #     at 1.0x the error-scale posterior landed in the 0.7% tail of its own prior (a clear
  #     prior-data conflict), at 1.5x the 3.7-7.1% tail, at 2.0x near 13% -- inside the bulk, with
  #     no further movement in the overfitting or delta metrics by 2.5x. Without the correction the
  #     factors cannot generate the observed excursion, the shortfall is booked as error, and
  #     stat_weak mimics nonstat by inflating noise rather than by fitting structure -- which is not
  #     the phenomenon this example is meant to show.
  #
  # ETA_FRAC_NONSTAT is the DGP's own level-error sd (2.0 against an innovation scale of 1.0), so
  # the correctly specified arm is handed the truth while the stationary arms are not. That is a
  # deliberate choice for the reference arm and belongs in the write-up.
  #
  # eta is on the LEVEL scale in BOTH branches: ife_named.stan applies the sqrt(2) differencing
  # inflation itself (errors_cov is normalized by 2v and error_precisions is halved), so these
  # fractions take no differencing correction and a caller-side one would double-count. An earlier
  # revision anchored the nonstationary arm on sd(diff(y)) instead, which silently multiplied both
  # its sigma and its eta by ~3; the reference arm then failed its own S1 check (p = 0.018 / 0.000 /
  # 0.001 on three datasets, against 0.944 / 0.232 / 0.807 here) with 95% predictive intervals
  # 200-260% of the observed data range. Hence the single anchor, and hence this comment.
  rms_y <- apply(fit_ys, 2, function(y) sqrt(mean(y^2)))
  sd_y <- apply(fit_ys, 2, sd)

  # SIGMA_MULT_*, ETA_FRAC_*, DELTA_FRAC, RHO_*, ALPHA_DIAG and the DGP constants all come from
  # ex1_config.r, which the figures and ex1_derive_scales.r read too. Nothing scale-related is
  # defined here any more; only the per-dataset quantities below, which carry this dataset's units.
  #
  # Default FALSE: the nonstationary arm ESTIMATES eta like the stationary arm. See the note at its
  # fit below for why, and what setting this changes.
  NONSTAT_FIX_ETA <- nzchar(Sys.getenv("NONSTAT_FIX_ETA"))


  # Both variables ARE sigma -- the vector the model receives -- so no call site has to remember
  # which multiple was applied where. The error scales read from rms_y, never from these; keeping
  # the two apart is what makes "eta is a fraction of the data scale" true as written rather than
  # true of one arm and not the other.
  overall_scales_stat    <- SIGMA_MULT_STAT * rms_y
  overall_scales_nonstat <- SIGMA_MULT_NONSTAT * rms_y
  eta_anchor <- mean(rms_y)
  delta_scale_ex1 <- DELTA_FRAC * mean(sd_y)


  # Draw every Stan seed up front, before any sample_model() call: cmdstanr's $sample() advances R's
  # global RNG, so a seed drawn after a fit would not be reproducible. Invariant: never derive a seed
  # after a fit. One seed per model PER ESCALATION ROUND, so a refit is reproducible too.
  fit_seeds <- matrix(sample.int(.Machine$integer.max, 2L * ESCALATE_MAX), nrow = 2L)

  fits <- list()

  fits$nonstat <- fit_with_escalation(
    list(
      alpha_diag = ALPHA_DIAG,
      N_units = N_units, T_times = T_times,
      # eta is ESTIMATED here, as it is for the stationary arms, under the same truncated-normal
      # form (location = scale) located at ETA_FRAC_NONSTAT * eta_anchor. Two reasons.
      #
      # First, fairness: this arm used to be handed a fixed error scale while the stationary arms
      # had to estimate theirs, so part of the contrast between them was fixed-vs-estimated rather
      # than stationary-vs-not. Estimating everywhere makes the comparison about stationarity
      # alone, which is what the example claims to be about.
      #
      # Second, and less obviously, the fixed version was NOT fixed at the truth. It was fixed at
      # ETA_FRAC_NONSTAT * mean(RMS(y_n)), which equals the DGP's err_sd = 2.0 only ON AVERAGE;
      # per dataset it followed the realised RMS. Measured over six datasets it ran 3.81 / 2.92 /
      # 1.17 / 1.71 / 2.85 / 1.65 against a truth that is exactly 2.0 every time -- a 3.3x spread
      # injected into the reference arm purely by the anchor. Estimating it recovers the truth
      # almost exactly instead: 1.99 / 2.04 / 2.00 / 2.03 / 1.82 / 2.00, with prior tails of
      # 0.28-0.81, i.e. sitting mid-prior rather than in a tail. Divergences fell from 1.17 to 0.17
      # per rep, and the four properties this study depends on all held: delta error 0.120 against
      # the stationary arms' 2.19, noise_abs 0.045 against their 0.229, coverage 0.975, S1 0.349.
      # (Six reps in fast mode -- indicative, not the full run.)
      #
      # NONSTAT_FIX_ETA=1 restores the fixed version, kept as an escape hatch for checking that a
      # result does not depend on this choice.
      overall_scales = overall_scales_nonstat,
      err_scale = if (NONSTAT_FIX_ETA) ETA_FRAC_NONSTAT * eta_anchor else 0,
      absolute_error = TRUE,
      err_scale_mean = if (NONSTAT_FIX_ETA) 0 else ETA_FRAC_NONSTAT * eta_anchor,
      err_scale_sd   = if (NONSTAT_FIX_ETA) 0 else ETA_FRAC_NONSTAT * eta_anchor,
      data = fit_ys,
      autocor_a = RHO_NONSTAT[1], autocor_b = RHO_NONSTAT[2],
      nonstationary = TRUE, num_treated = num_treated_ex1, delta_scale = delta_scale_ex1,
      fit_scales = FALSE,
      type = "posterior", K_latent = K_latent, ad = 0.8,
      iter = EX1_ITER, iter_warm = EX1_WARM,
      n_chains = 3, pathfinder_init = TRUE
    ),
    seeds = fit_seeds[1, ], label = sprintf("unit %d nonstat", i), progress_log = progress_log, ladder = EX1_LADDER
  )
  fits$nonstat$name <- "nonstat"

  fits$stat <- fit_with_escalation(
    list(
      alpha_diag = ALPHA_DIAG,
      N_units = N_units, T_times = T_times,
      overall_scales = overall_scales_stat, err_scale = 0, absolute_error = TRUE,
      err_scale_mean = ETA_FRAC_STAT * eta_anchor,
      err_scale_sd = ETA_FRAC_STAT * eta_anchor,
      data = fit_ys,
      autocor_a = RHO_STAT[1], autocor_b = RHO_STAT[2],
      nonstationary = FALSE, num_treated = num_treated_ex1, delta_scale = delta_scale_ex1,
      fit_scales = FALSE,
      type = "posterior", K_latent = K_latent, ad = 0.8,
      iter = EX1_ITER, iter_warm = EX1_WARM,
      n_chains = 3, pathfinder_init = TRUE
    ),
    seeds = fit_seeds[2, ], label = sprintf("unit %d stat", i), progress_log = progress_log, ladder = EX1_LADDER
  )
  fits$stat$name <- "stat"


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
            (fit_ys[t, n] >= y_bounds[1]) &&
              (fit_ys[t, n] <= y_bounds[2])
          pred_width[t, n] <- (y_bounds[2] - y_bounds[1]) / (max(fit_ys[, n]) - min(fit_ys[, n]))
        }
      }
      res$pred_perc <- mean(pred_inc)
      res$pred_width <- mean(pred_width)

      # The estimated observation-error scale, recorded so scale changes can be judged directly
      # instead of inferred from the plots, together with how far into its own prior's upper tail
      # the posterior median sits -- a prior-data conflict here is what says a sigma multiple is
      # too small (see SIGMA_MULT_STAT above).
      #
      # Reported as eta, on the DATA's scale, because that is the scale the priors are written in.
      # pfit$err_scale holds draws of tau[1] = eta / sigma[1], so multiply back by this arm's own
      # sigma[1]. The previous version compared tau directly against the eta priors, which is off
      # by a factor of sigma[1] (~13) and so reported a near-constant number carrying no
      # information about the fit.
      # Each arm's OWN sigma[1]. stat_weak is on sd_y, not overall_scales_stat, so a two-way
      # nonstat/other split would misreport its eta by a factor of 2 * RMS / sd (~3.4x here).
      sigma_1 <- if (pfit$name == "nonstat") overall_scales_nonstat[1] else overall_scales_stat[1]
      eta_draws <- pfit$err_scale * sigma_1
      res$eta_med <- median(eta_draws)
      res$eta_q95 <- unname(quantile(eta_draws, 0.95))
      # NA only when the arm's eta is FIXED, and so has no prior to sit in a tail of.
      eta_prior_loc <- if (pfit$name == "nonstat") {
        if (NONSTAT_FIX_ETA) NA_real_ else ETA_FRAC_NONSTAT * eta_anchor
      } else ETA_FRAC_STAT * eta_anchor
      # location and scale are equal for every arm, by construction above
      res$eta_prior_tail <- if (is.na(eta_prior_loc)) NA_real_ else
        (1 - pnorm(res$eta_med, eta_prior_loc, eta_prior_loc)) /
          (1 - pnorm(0, eta_prior_loc, eta_prior_loc))

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
      true_ys <- true_ys_perm
      fit_means <- apply(pfit$y_means, c(2, 3), mean)
      pre_times <- seq_len(T_times - num_treated_ex1)
      fit_err <- fit_means[pre_times, 1] - true_ys[pre_times, 1]
      unit_noise <- fit_ys[pre_times, 1] - true_ys[pre_times, 1]
      res$noise_abs_tr <- as.numeric(cov(fit_err, unit_noise) / var(unit_noise))

      return(res)
    }) |>
    list_flatten()

  pns_means <- apply(fits$nonstat$y_means, c(2, 3), mean)
  pstat_means <- apply(fits$stat$y_means, c(2, 3), mean)

  # Show a single unit's fits per rep
  plot_unit <- 2
  # One posterior predictive replicate from the stationary arm, overlaid in red. Everything here is
  # in ANCHOR-PERMUTED column order -- fit_ys, both posterior means, and y_pred alike. y_pred is
  # Y_latent plus observation noise, which is what S1 sees and what the mean lines do not show.
  #
  # plot_unit = 2 is AN untreated unit, not a distinguished one: S1 now averages over every
  # untreated unit, so no single column is the one the statistic is "about". Note this figure
  # therefore says nothing about the treated unit's counterfactual, which is where the delta error
  # lives -- set plot_unit = 1 to look at that instead.
  # Set pred_rep = NULL to drop the overlay.
  # Gate the per-rep figure. Unconditional, a 200-rep run wrote 200 PNGs -- I/O and clutter per rep
  # for a diagnostic nobody inspects in bulk; ex2 gates this the same way. Raise PLOT_REPS for more.
  #
  # Note `i` is the DATASET index, not the rep counter, so this keeps the low-numbered datasets
  # rather than the first PLOT_REPS reps. That is deliberate: the same datasets are then plotted
  # across runs with different rep counts, which makes two runs comparable figure-by-figure.
  if (i <= PLOT_REPS) {
    fit_plot <- plot_post_fits_stat(fit_ys, pns_means, pstat_means, unit = plot_unit,
      pred_rep = fits$stat$y_pred[1, , plot_unit])
    ggsave(
      fit_plot,
      file = paste0("../figs/sim_stat_figs/post_fit_plot_u", plot_unit, "_", i, ".png"),
      create.dir = TRUE
    )
  }

  if (post_check) {
    check_plot <- plot_data_matrix_post(fit_ys, fits$stat$y_pred)
    ggsave(
      check_plot,
      file = paste0("../figs/sim_stat_figs/check_plot_", i, ".pdf"),
      device = "pdf", height = 4, width = 8, create.dir = TRUE
    )
  }

  return(res)
}

run_sim_study_stat <- function(K_latent = K_LATENT, reps, seed, post_check = FALSE) {
  # Seed the master-process RNG, and draw the dataset selection and the prior-
  # predictive Stan seed *before* generating test_data, so they are not perturbed
  # by cmdstanr's $sample() (which advances R's RNG).
  set.seed(seed)
  study_units <- sample.int(2 * reps, size = reps, replace = FALSE)
  pp_seed <- sample.int(.Machine$integer.max, 1)

  test_data <- sample_model(
    # sigma = 1 and err_sd = 2 are the truth this study measures its arms against. The nonstat arm
    # is handed both (see SIGMA_MULT_NONSTAT above), so "correctly specified" here means the true
    # scales, the true rho prior Beta(8, 2) and the same alpha_diag -- not merely the true
    # functional form.
    #
    # absolute_error = TRUE though sigma is constant, which makes it numerically identical to the
    # ratio mode this used to be in (verified bit-identical). It is declared anyway so that the
    # DGP states the same parametrization as every fitted arm: in ratio mode err_scale would be
    # read as a MULTIPLE of overall_scales, so making overall_scales data-dependent here would
    # silently rescale the error. That is exactly how the nonstationary arm's scales went ~3x
    # wrong once already.
    overall_scales = rep(DGP_SIGMA, N_UNITS), err_scale = DGP_ETA, absolute_error = TRUE,
    alpha_diag = ALPHA_DIAG,
    autocor_a = DGP_RHO[1], autocor_b = DGP_RHO[2],
    nonstationary = TRUE, num_treated = 0,
    type = "prior_pred", K_latent = K_latent,
    iter = 2 * reps, seed = pp_seed
  )

  exp_vars <- c("run_sim_stat", "worker_progress", "sample_model", "ife_mod", "plot_post_fits_stat", "plot_data_matrix_post",
    "pathfinder_inits", "draw_to_init", "PF_PARAM_BASES", "anchor_order",
    "fit_with_escalation", "escalation_ladder", "ESCALATE_MAX", "EX1_LADDER",
    "EX1_ITER", "EX1_WARM", "PLOT_REPS")
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
  ckpt_check_fingerprint(ckpt_dir, c("nonstat", "stat"), "ex1_config.r")
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
sim_study_stat <- run_sim_study_stat(K_latent = K_LATENT, study_reps, seed = 40318)

stopCluster(cl)

# Save the raw study results; numeric summaries and plots are produced by
# ex1_sim_study_summary.r (run it to view the results).
# Mode-keyed output, so a quick fast-mode check cannot silently destroy the full study's results.
out_file <- if (STUDY_MODE == "fast") "sim_study_ns_fast.RData" else "sim_study_ns.RData"
save(sim_study_stat, file = out_file)
cat(sprintf("Results saved to %s -- run `Rscript ex1_sim_study_summary.r %s` to summarize.\n",
  out_file, out_file))
