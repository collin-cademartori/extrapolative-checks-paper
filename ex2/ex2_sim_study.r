## This file runs the simulation study for the intercepts example,
## fitting the models with and without intercepts to data generated from an
## adversarial process whereby some units are simulated to have different
## long run means than the treated but spuriously correlated in the pre-
## treatment period only.

library(foreach)
library(doParallel)
library(doRNG)
library(purrr)

# Command line: Rscript ex2_sim_study.r [n_cores] [mode] [reps]
#   n_cores  worker count (default: half the physical cores, less one)
#   mode     "full" (default) or "fast"
#   reps     number of simulation replications to perform
#
# In fast mode, chains are shorter and never escalated.
.args <- commandArgs(trailingOnly = TRUE)
requested_cores <- suppressWarnings(as.integer(.args[1]))
STUDY_MODE <- if (length(.args) >= 2 && tolower(.args[2]) %in% c("fast", "f")) "fast" else "full"
.reps_arg <- suppressWarnings(as.integer(.args[3]))
n_cores <- if (!is.na(requested_cores) && requested_cores >= 1) {
  requested_cores
} else {
  max(1, round(detectCores() / 2) - 1)
}
# Redirect temporary files to filesystem to reduce pressure on RAM from
# maintaining many independent R instances when running many parallel workers.
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

# Set PSOCK workers' working directories to current working directory so
# relative paths resolve correctly.
invisible(clusterCall(cl, setwd, getwd()))

# Attach the workers' packages, suppressing load messages.
invisible(clusterCall(cl, source, "ex2_config.r"))
invisible(clusterCall(cl, source, "ex2_dgp.r"))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(ggplot2)
  library(doRNG)
  library(purrr)
})))

source("ex2_config.r")
source("ex2_dgp.r")
source("../sample_model.r")
source("../pathfinder_init.r")
source("../plotting.r")

# Workers call this function to write progress to a common log file
worker_progress <- function(label, logfile = "progress.log") {
  n <- get0(".worker_done", envir = globalenv(), ifnotfound = 0L) + 1L
  assign(".worker_done", n, envir = globalenv())
  cat(sprintf(
    "[%s] [worker %d] %3d done | %s\n",
    format(Sys.time(), "%H:%M"), Sys.getpid(), n, label
  ), file = logfile, append = TRUE)
}

# Runs with bad diagnostics (high R-hat or divergences) are refit with more iterations and 
# a higher adapt_delta, up to two extra rounds. 
# These constant defines the number of increased iterations.
EX2_LADDER <- if (STUDY_MODE == "fast") escalation_ladder(integer(0), integer(0)) else
  escalation_ladder(iter = c(3000L, 6000L), warm = c(1500L, 2000L))
ESCALATE_MAX <- EX2_LADDER$max_rounds   # seeds are drawn one per fit per round

# Base sampling length.
EX2_ITER <- if (STUDY_MODE == "fast") 500L else 1500L
EX2_WARM <- if (STUDY_MODE == "fast") 500L else 1000L

# Identifies the model configuration a checkpoint set was written under, so a set from a different
# one is not resumed into this run.
ckpt_fingerprint <- function(arms, config_file) {
  paste(c(paste(sort(arms), collapse = ","),
          unname(tools::md5sum(config_file))), collapse = " ")
}

# Refuse to resume a set written under a different configuration.
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

# Column ordering for the triangular (Cholesky) loadings: the treated unit stays first, then the
# untreated columns most orthogonal to those already chosen.
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

# Map vector of per-unit summaries output from Stan model back to original
# order (before reordering applied by [anchor_order] above).
unpermute_untreated <- function(v, perm) {
  out <- numeric(length(v))
  out[perm[-1] - 1] <- v
  out
}


run_sim_intercepts <- function(N_comp, level, K_latent = K_LATENT, rep_i = NA, plot_iters = 0,
                               progress_log = NULL) {

  N_unc <- DGP_N_UNITS - 1 - DGP_N_COMP_TRUE - N_comp
  gen <- sim_model_intercepts(
    N_unc = N_unc, N_comp_true = DGP_N_COMP_TRUE, N_comp_spur = N_comp, K_unc = DGP_K_UNC,
    sim = DGP_SIM, level_offset = level, T_times = DGP_T_TIMES
  )
  test_ys <- gen$Y

  N_units <- ncol(test_ys)
  T_times <- nrow(test_ys)

  # Re-order columns (treated stays first) so the leading K x K loading
  # block is full rank.
  perm <- anchor_order(test_ys, K_latent)
  stopifnot(perm[1] == 1)
  fit_ys <- test_ys[, perm]

  fit_seeds <- matrix(sample.int(.Machine$integer.max, 2L * ESCALATE_MAX), nrow = 2L)

  fits <- list()

  # ---- scales ---------------------------------------------------------------------------------
  pre_y <- fit_ys[seq_len(nrow(fit_ys) - NUM_TREATED), , drop = FALSE]
  sd_y <- apply(pre_y, 2, sd)
  eta_anchor <- mean(sd_y)
  eta_loc <- ETA_FRAC_EX2 * eta_anchor
  eta_scale <- ETA_CV_EX2 * eta_loc
  delta_scale_ex2 <- DELTA_FRAC_EX2 * sd(pre_y[, 1])
  int_loc_ex2 <- mean(pre_y)
  int_scale_ex2 <- INT_FRAC * sd(colMeans(pre_y))

  overall_scales <- apply(pre_y, 2, \(x) sqrt(mean(x ^ 2)))
  fits$no_ints <- fit_with_escalation(
    list(
      N_units = ncol(fit_ys), T_times = nrow(fit_ys), K_latent = K_latent,
      overall_scales = overall_scales,
      err_scale = 0,
      err_scale_mean = eta_loc, err_scale_sd = eta_scale,
      data = fit_ys,
      autocor_a = RHO_EX2[1], autocor_b = RHO_EX2[2],
      nonstationary = FALSE, num_treated = NUM_TREATED, delta_scale = delta_scale_ex2,
      include_factor_means = TRUE,
      alpha_diag = ALPHA_DIAG, pathfinder_init = TRUE,
      type = "posterior", quiet = TRUE, ad = 0.8,
      iter = EX2_ITER, iter_warm = EX2_WARM,
      n_chains = 3
    ),
    seeds = fit_seeds[1, ],
    label = sprintf("level %g num_comp %d rep %d no_ints", level, N_comp, rep_i),
    progress_log = progress_log, ladder = EX2_LADDER
  )
  fits$no_ints$name <- "no_ints"

  overall_sds <- sd_y
  fits$ints <- fit_with_escalation(
    list(
      N_units = ncol(fit_ys), T_times = nrow(fit_ys), K_latent = K_latent,
      overall_scales = overall_sds,
      err_scale = 0,
      err_scale_mean = eta_loc, err_scale_sd = eta_scale,
      data = fit_ys,
      autocor_a = RHO_EX2[1], autocor_b = RHO_EX2[2],
      nonstationary = FALSE, num_treated = NUM_TREATED, delta_scale = delta_scale_ex2,
      include_ints = TRUE, int_scale = int_scale_ex2, int_loc = int_loc_ex2,
      alpha_diag = ALPHA_DIAG, pathfinder_init = TRUE,
      type = "posterior", quiet = TRUE, ad = 0.8,
      iter = EX2_ITER, iter_warm = EX2_WARM,
      n_chains = 3
    ),
    seeds = fit_seeds[2, ],
    label = sprintf("level %g num_comp %d rep %d ints", level, N_comp, rep_i),
    progress_log = progress_log, ladder = EX2_LADDER
  )
  fits$ints$name <- "ints"

  for (m in c("no_ints", "ints")) {
    fits[[m]]$cor_sq <- unpermute_untreated(fits[[m]]$cor_sq, perm)
    fits[[m]]$abs_cors_err <- unpermute_untreated(fits[[m]]$abs_cors_err, perm)
  }

  res <- fits |>
    map(function(pfit) {
      res <- list()

      # Posterior-predictive interval coverage, per (time, unit). Compared against fit_ys.
      stat_y_pred <- pfit$y_pred
      pred_inc <- matrix(NA, nrow = T_times, ncol = N_units)
      pred_width <- matrix(NA, nrow = T_times, ncol = N_units)
      for (n in 1:N_units) {
        for (t in 1:T_times) {
          y_bounds <- quantile(stat_y_pred[, t, n], c(0.005, 0.995))
          pred_inc[t, n] <- (fit_ys[t, n] >= y_bounds[1]) && (fit_ys[t, n] <= y_bounds[2])
          pred_width[t, n] <- (y_bounds[2] - y_bounds[1]) / (max(fit_ys[, n]) - min(fit_ys[, n]))
        }
      }
      res$pred_perc <- mean(pred_inc)
      res$pred_width <- mean(pred_width)

      sigma_1 <- if (pfit$name == "no_ints") overall_scales[1] else overall_sds[1]
      eta_draws <- pfit$err_scale * sigma_1
      res$eta_med <- median(eta_draws)
      res$eta_prior_tail <- (1 - pnorm(res$eta_med, eta_loc, eta_scale)) /
        (1 - pnorm(0, eta_loc, eta_scale))

      res$loc_cor_pval <- pfit$loc_cor_pval

      absz <- abs(pfit$effect_means / pfit$effect_sds)
      res[paste0("absz_", seq_along(absz))] <- absz

      pmean <- pfit$effect_means
      res[paste0("mean_", seq_along(pmean))] <- pmean

      psds <- pfit$effect_sds
      res[paste0("sd_", seq_along(psds))] <- psds

      # Exact posterior P(delta_t > 0) and two-sided tail area of the true (zero) effect.
      res[paste0("ppos_", seq_along(psds))] <- pfit$effect_p_pos
      res[paste0("tail_", seq_along(psds))] <- pfit$effect_tail

      cor_sq <- pfit$cor_sq
      res[paste0("cor_sq_", seq_along(cor_sq))] <- cor_sq

      acors_err <- pfit$abs_cors_err
      res[paste0("acor_err_", seq_along(acors_err))] <- acors_err

      pred_mad <- pfit$mean_abs_diffs
      res$pred_mad <- pred_mad

      # Convergence diagnostics recorded for every fit, so the run can be audited.
      # rhat_M checks convergence of M = Lambda %*% Phi, the product of loadings and
      # factors, through which the likelihood depends on the loadings and factors.
      sdg <- pfit$sampler_diag
      res$rhat_max <- sdg$rhat_max
      res$rhat_M <- sdg$rhat_M
      res$rhat_estimands <- sdg$rhat_estimands
      res$rhat_loadings <- sdg$rhat_loadings
      res$rhat_cor_sq <- sdg$rhat_cor_sq
      res$n_div <- sdg$n_div
      res$n_tree <- sdg$n_tree
      res$ebfmi_min <- sdg$ebfmi_min
      res$lp_gap_max <- sdg$lp_gap_max
      res$ess_delta1 <- sdg$ess_delta1
      # How much computation this fit needed: 1 = met the criterion on the first attempt.
      res$n_rounds <- pfit$n_rounds
      res$final_iter <- pfit$final_iter
      res$final_warm <- pfit$final_warm
      res$final_ad <- pfit$final_ad

      return(res)
    }) |>
    list_flatten()

  if (!is.na(rep_i) && rep_i <= plot_iters) {
    for (m in c("no_ints", "ints")) {
      int_plot <- plot_intercepts_fits(
        test_ys,
        cor_sq = fits[[m]]$cor_sq,
        groups = gen$groups, num_treated = NUM_TREATED
      )
      ggsave(
        int_plot,
        file = sprintf("../figs/sim_int_figs/int_fit_%s_lv%g_nc%d_rep%d.png", m, level, N_comp, rep_i),
        width = 5, height = 4, create.dir = TRUE
      )
    }
  }

  return(res)
}

run_sim_study_intercepts <- function(K_latent = K_LATENT, reps, N_comps, levels, seed, plot_iters = 3) {

  exp_vars <- c(
    "run_sim_intercepts",
    "worker_progress", "sample_model", "ife_mod", "plot_intercepts_fits",
    "anchor_order", "unpermute_untreated",
    "pathfinder_inits", "draw_to_init", "PF_PARAM_BASES",
    "fit_with_escalation", "escalation_ladder", "ESCALATE_MAX", "EX2_LADDER",
    "EX2_ITER", "EX2_WARM"
  )
  exp_packages <- c("cmdstanr", "posterior", "ggplot2", "dplyr")

  grid <- expand.grid(rep = seq_len(reps), N_comp = N_comps, level = levels)

  cat(sprintf(
    paste0(
      "\n=== Example 2 simulation study ===\n",
      "  conditions : level {%s} x num_comp {%s}  (%d)\n",
      "  reps/cond  : %d\n",
      "  tasks      : %d  (2 model fits each)\n",
      "  workers    : %d   seed: %d\n\n"
    ),
    paste(levels, collapse = ", "), paste(N_comps, collapse = ", "),
    length(levels) * length(N_comps), reps, nrow(grid),
    getDoParWorkers(), seed
  ))

  progress_log <- file.path(getwd(), "progress.log")
  t0 <- Sys.time()

  # Each task writes its own row as it finishes.
  ckpt_dir <- file.path(getwd(),
    sprintf("ckpt_ints_%s_seed%d_n%d", STUDY_MODE, seed, nrow(grid)))
  dir.create(ckpt_dir, showWarnings = FALSE)
  ckpt_check_fingerprint(ckpt_dir, c("no_ints", "ints"), "ex2_config.r")
  n_resume <- length(list.files(ckpt_dir, pattern = "^task_.*\\.rds$"))
  if (n_resume > 0) {
    cat(sprintf("  resuming: %d of %d tasks already checkpointed in %s\n\n",
      n_resume, nrow(grid), basename(ckpt_dir)))
  } else {
    cat("", file = progress_log)
  }

  study_res <-
    foreach(
      rep_i = grid$rep, N_comp = grid$N_comp, level = grid$level, task_i = seq_len(nrow(grid)),
      .combine = function(...) dplyr::bind_rows(...),
      .export = exp_vars, .packages = exp_packages,
      .options.RNG = seed
    ) %dorng% {
      ckpt_file <- file.path(ckpt_dir,
        sprintf("task_%05d_lv%g_nc%d_rep%04d.rds", task_i, level, N_comp, rep_i))
      if (file.exists(ckpt_file)) {
        readRDS(ckpt_file)
      } else {
        unit_res <- tryCatch(
          as.data.frame(c(
            run_sim_intercepts(
              N_comp = N_comp, level = level, K_latent = K_latent,
              rep_i = rep_i, plot_iters = plot_iters, progress_log = progress_log
            ),
            list(level = level, num_comp = N_comp)
          )),
          error = function(e) {
            cat(sprintf("[%s] level %g num_comp %d rep %d FAILED: %s\n",
              format(Sys.time(), "%H:%M"), level, N_comp, rep_i, conditionMessage(e)),
              file = progress_log, append = TRUE)
            data.frame(level = level, num_comp = N_comp, error = conditionMessage(e))
          }
        )
        unit_res$rep <- rep_i
        if (is.null(unit_res$error)) unit_res$error <- NA_character_
        unit_res$failed <- !is.na(unit_res$error)
        saveRDS(unit_res, ckpt_file)
        worker_progress(sprintf("level %g  num_comp %d  rep %d", level, N_comp, rep_i), logfile = progress_log)
        unit_res
      }
    }

  n_failed <- sum(study_res$failed, na.rm = TRUE)
  cat(sprintf(
    "--- study complete: %d tasks in %.1f min%s ---\n",
    nrow(grid), as.numeric(difftime(Sys.time(), t0, units = "mins")),
    if (n_failed > 0) sprintf("; %d FAILED (see the `error` column and progress.log)", n_failed) else ""
  ))

  return(study_res)
}

study_reps <- if (!is.na(.reps_arg)) .reps_arg else if (STUDY_MODE == "fast") 6L else 200L
cat(sprintf("\n=== mode: %s | reps/cond: %d | iter/warm: %d/%d | escalation: %s ===\n",
  STUDY_MODE, study_reps, EX2_ITER, EX2_WARM,
  if (EX2_LADDER$max_rounds > 1) sprintf("up to %d rounds", EX2_LADDER$max_rounds) else "OFF"))
if (STUDY_MODE == "fast")
  cat("    FAST MODE -- for specification search only. Elevated rhat_M is expected here and says\n",
      "   nothing about the specification. Do not report these numbers.\n", sep = "")
sim_study_ints <- run_sim_study_intercepts(
  reps = study_reps,
  N_comps = DGP_N_COMP_SPUR,
  levels = DGP_LEVEL,
  K_latent = K_LATENT,
  seed = 52918,
  plot_iters = 50
)

stopCluster(cl)

# Save the raw study results. Numeric summaries and plots are produced by
# ex2_sim_study_summary.r
out_file <- if (STUDY_MODE == "fast") "sim_study_ints_fast.RData" else "sim_study_ints.RData"
save(sim_study_ints, file = out_file)
cat(sprintf("Results saved to %s -- run `Rscript ex2_sim_study_summary.r %s` to summarize.\n",
  out_file, out_file))
