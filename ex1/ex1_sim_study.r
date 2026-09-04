## Simulation study for the nonstationary example. Fits a correctly specified nonstationary model
## and a misspecified stationary one to datasets drawn from the nonstationary model's prior
## predictive distribution.

library(foreach)
library(doParallel)
library(doRNG)
library(purrr)

# Command line: Rscript ex1_sim_study.r [n_cores] [mode] [reps]
#   n_cores  worker count (default: half the physical cores, less one)
#   mode     "full" (default) or "fast"
#   reps     overrides the mode's default rep count
#
# In fast mode, chains are shorter and never escalated (re-run when sampler diagnostics unhappy). 
.args <- commandArgs(trailingOnly = TRUE)
requested_cores <- suppressWarnings(as.integer(.args[1]))
STUDY_MODE <- if (length(.args) >= 2 && tolower(.args[2]) %in% c("fast", "f")) "fast" else "full"
.reps_arg <- suppressWarnings(as.integer(.args[3]))
n_cores <- if (!is.na(requested_cores) && requested_cores >= 1) {
  requested_cores
} else {
  max(1, round(detectCores() / 2) - 1)
}
# Redefine temporary directories to ensure temporary files are written to file system rather
# than RAM, decreasing memory usage in runs with many parallel samplers.
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

invisible(clusterCall(cl, setwd, getwd()))
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

# Per-worker progress: each worker appends a running count of completed tasks to a
# shared log.
worker_progress <- function(label, logfile = "progress.log") {
  n <- get0(".worker_done", envir = globalenv(), ifnotfound = 0L) + 1L
  assign(".worker_done", n, envir = globalenv())
  cat(sprintf(
    "[%s] [worker %d] %3d done | %s\n",
    format(Sys.time(), "%H:%M"), Sys.getpid(), n, label
  ), file = logfile, append = TRUE)
}

# Fits that miss the convergence criterion are refit with more iterations and a higher
# adapt_delta, up to two extra rounds. The mechanism is in sample_model.r; these are ex1's rungs.
EX1_LADDER <- if (STUDY_MODE == "fast") escalation_ladder(integer(0), integer(0)) else
  escalation_ladder(iter = c(8000L, 12000L), warm = c(2000L, 3000L))
ESCALATE_MAX <- EX1_LADDER$max_rounds   # seeds are drawn one per fit per round

# How many reps get a per-rep fit figure.
PLOT_REPS <- 25L

EX1_ITER <- if (STUDY_MODE == "fast") 500L else 2000L
EX1_WARM <- if (STUDY_MODE == "fast") 500L else 500L

# Column ordering for the triangular (Cholesky) loadings: the treated unit stays first, then the
# untreated columns most orthogonal to those already chosen. A reparameterization only -- the fitted
# means, treatment effect and check statistics do not depend on it.

# Identifies the model configuration a checkpoint set was written under, so a set from a different
# one is not resumed into this run.
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

  # Order the columns for the triangular loadings. The treated unit stays in column 1, which
  # everything downstream indexing it by position relies on.
  perm <- anchor_order(test_ys, K_latent + 1)
  stopifnot(perm[1] == 1)
  fit_ys <- test_ys[, perm]
  true_ys_perm <- test_data$ys_latent[i, , perm]
  num_treated_ex1 <- NUM_TREATED

  # ---- scales -------------------------------------------------------------------------------
  # Every scale constant is a fraction of one anchor, RMS(y_n); the values are in ex1_config.r and
  # their derivation in ex1_derive_scales.r.
  #
  # What the example demonstrates. The nonstationary arm shares the DGP's rho prior and so receives
  # the DGP's own sigma and eta: the reference arm is handed the truth. The stationary arm is handed
  # an error smaller than the truth, and that one number does both jobs -- it produces the
  # overfitting, and it is what the prior predictive check on |cor(y, t)| requires, since a large
  # error attenuates a linear trend and only a small one lets a stationary model look nonstationary
  # enough to be plausible. Under a wrong model, a prior that survives its own predictive check is a
  # prior that overfits.
  rms_y <- apply(fit_ys, 2, function(y) sqrt(mean(y^2)))
  sd_y <- apply(fit_ys, 2, sd)

  # Set NONSTAT_FIX_ETA=1 to fix the nonstationary arm's error scale instead of estimating it.
  NONSTAT_FIX_ETA <- nzchar(Sys.getenv("NONSTAT_FIX_ETA"))

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
      # eta is estimated here, as it is for the stationary arm, so the contrast between them is
      # stationarity alone rather than fixed-versus-estimated.
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

      # The estimated observation-error scale, reported as eta on the data's scale, with how far
      # into its own prior's upper tail the posterior median sits.
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
      res$ebfmi_min <- sdg$ebfmi_min
      res$n_tree <- sdg$n_tree

      # Overfitting of the treated unit's pre-treatment window, the basis the counterfactual is
      # extrapolated from:
      #   noise_abs_tr = cov(fitted - truth, observed - truth) / var(observed - truth)
      # the fraction of that unit's noise absorbed into its fitted signal; 1 = interpolation,
      # 0 = noise ignored.
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
  # One posterior predictive replicate from the stationary arm, overlaid in red; y_pred is the
  # latent signal plus observation noise, which is what S1 sees and the mean lines do not show.
  # plot_unit selects an untreated unit; set it to 1 for the treated unit's counterfactual, and
  # pred_rep = NULL to drop the overlay.
  # Only the first PLOT_REPS datasets get a figure. `i` is the dataset index, not the rep counter,
  # so runs with different rep counts plot the same datasets.
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
    # sigma = 1 and err_sd = 2 are the truth the arms are measured against. The nonstationary arm
    # is handed both, along with the true rho prior and the same alpha_diag, so "correctly
    # specified" here means more than the true functional form.
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
      "  tasks   : %d  (2 model fits each)\n",
      "  workers : %d   seed: %d\n\n"
    ),
    reps, reps, getDoParWorkers(), seed
  ))
  # Absolute log path, so workers write it where the master expects regardless of
  # their working directory.
  progress_log <- file.path(getwd(), "progress.log")
  t0 <- Sys.time()

  # Each rep writes its own row as it finishes, and a rerun picks up the rows that already exist.
  # A rep that errors records failed = TRUE and the study continues. The directory is keyed to mode,
  # seed and rep count; delete it to force a full recompute.
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
