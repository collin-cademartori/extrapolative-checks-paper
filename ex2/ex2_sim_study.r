## This file runs the simulation study for the intercepts example, fitting the
## fitting the models with and without intercepts to data generated from an
## adversarial process whereby some units are simulated to have different
## long run means than the treated but spuriously correlation in the pre-
## treatment period only.

library(foreach)
library(doParallel)
library(doRNG)
library(purrr)

# Command line: Rscript ex2_sim_study.r [n_cores] [mode] [reps]
#   n_cores  worker count (default: half the physical cores, less one)
#   mode     "full" (default) or "fast"
#   reps     overrides the mode's default rep count
#
# In fast mode, chains are shorter and never escalated. Elevated rhat_M is expected there and says
# nothing about the specification; never report fast-mode numbers.
.args <- commandArgs(trailingOnly = TRUE)
requested_cores <- suppressWarnings(as.integer(.args[1]))
STUDY_MODE <- if (length(.args) >= 2 && tolower(.args[2]) %in% c("fast", "f")) "fast" else "full"
.reps_arg <- suppressWarnings(as.integer(.args[3]))
n_cores <- if (!is.na(requested_cores) && requested_cores >= 1) {
  requested_cores
} else {
  max(1, round(detectCores() / 2) - 1)
}
# Redirect temporary files to the file system rather than RAM, so that a long run with many
# parallel samplers cannot exhaust a tmpfs /tmp. Override by exporting either variable before
# launching.
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
# running a different value than the master is exactly the drift ex2_config.r exists to prevent.
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

# Fits that miss the convergence criterion are refit with more iterations and a higher
# adapt_delta, up to two extra rounds. The mechanism is shared with ex1 and lives in
# sample_model.r; only the rungs below are per-example. rhat_M is the criterion that binds in
# practice here -- the ess and divergence thresholds have never triggered.
EX2_LADDER <- if (STUDY_MODE == "fast") escalation_ladder(integer(0), integer(0)) else
  escalation_ladder(iter = c(6000L, 9000L), warm = c(2500L, 4000L))
ESCALATE_MAX <- EX2_LADDER$max_rounds   # seeds are drawn one per fit per round

# Base sampling length, per mode. Both arms share these.
EX2_ITER <- if (STUDY_MODE == "fast") 500L else 1500L
EX2_WARM <- if (STUDY_MODE == "fast") 500L else 1000L

# Column ordering for the triangular (Cholesky) loadings: the treated unit stays first, then the
# untreated columns most orthogonal to those already chosen. A reparameterization only -- the
# fitted means, treatment effect and check statistics do not depend on it, though the per-unit
# outputs need mapping back (see unpermute_untreated).
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

# Map a per-untreated-unit vector returned from a fit on anchor-permuted data back to
# the original unit order. `perm` is the column permutation (perm[1] == 1, the treated
# unit); element j of `v` belongs to permuted column j + 1 = original column perm[j+1].
unpermute_untreated <- function(v, perm) {
  out <- numeric(length(v))
  out[perm[-1] - 1] <- v
  out
}


run_sim_intercepts <- function(N_comp, level, K_latent = K_LATENT, rep_i = NA, plot_iters = 0,
                               progress_log = NULL) {
  # Fixed total of 8 units; N_comp spurious comparators trade off against the
  # uncorrelated fillers (1 treated + 2 true + N_comp spurious + N_unc = 8), so the
  # swept quantity is the *share* of units spuriously correlated with the treated.
  N_unc <- DGP_N_UNITS - 1 - DGP_N_COMP_TRUE - N_comp
  gen <- sim_model_intercepts(
    N_unc = N_unc, N_comp_true = DGP_N_COMP_TRUE, N_comp_spur = N_comp, K_unc = DGP_K_UNC,
    sim = DGP_SIM, level_offset = level, T_times = DGP_T_TIMES
  )
  test_ys <- gen$Y

  N_units <- ncol(test_ys)
  T_times <- nrow(test_ys)

  # Fit on anchor-ordered columns (treated stays first) so the leading K x K loading
  # block is full rank; per-unit outputs are mapped back to the original order after the
  # fits. The plots below use the original test_ys.
  perm <- anchor_order(test_ys, K_latent)
  # perm[1] == 1 by construction. Everything downstream that indexes the treated unit by position
  # (delta, and unpermute_untreated's assumption that element j of a per-untreated vector belongs to
  # permuted column j + 1) depends on it, so make the dependency explicit rather than implicit.
  stopifnot(perm[1] == 1)
  fit_ys <- test_ys[, perm]

  # Draw every Stan seed up front, before any sample_model() call: cmdstanr's
  # $sample() advances R's global RNG, so a seed drawn after a fit would not be
  # reproducible. Invariant: never derive a seed after a fit. One seed per model PER ESCALATION
  # ROUND, so a refit is reproducible too.
  fit_seeds <- matrix(sample.int(.Machine$integer.max, 2L * ESCALATE_MAX), nrow = 2L)

  fits <- list()

  # ---- scales ---------------------------------------------------------------------------------
  # sigma differs between the arms by design: no_ints must produce each unit's level from its factor
  # means, so it is anchored on RMS(y_n), while ints has gamma for the level and is anchored on
  # sd(y_n). ex1's 2 x RMS correction does not apply here, since both fits have a level mechanism.
  #
  # eta does NOT differ between the arms, which is the point of this study: it is one shared
  # observation-error sd on the data's own scale, so the justified sigma difference cannot leak into
  # the error scale, where nothing would justify it. Values come from ex2_config.r and their
  # derivation from ex2_derive_scales.r.
  sd_y <- apply(fit_ys, 2, sd)
  eta_anchor <- mean(sd_y)
  eta_loc <- ETA_FRAC_EX2 * eta_anchor
  eta_scale <- ETA_CV_EX2 * eta_loc
  # Shared effect-prior scale, for the same reason eta is shared: see ex2_config.r.
  delta_scale_ex2 <- DELTA_FRAC_EX2 * eta_anchor
  # Unit-intercept prior, anchored on the data rather than fixed: see ex2_config.r.
  int_loc_ex2 <- mean(fit_ys)
  int_scale_ex2 <- INT_FRAC * sd(colMeans(fit_ys))

  overall_scales <- apply(fit_ys, 2, \(x) sqrt(mean(x ^ 2)))
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
      n_chains = 4
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
      n_chains = 4
    ),
    seeds = fit_seeds[2, ],
    label = sprintf("level %g num_comp %d rep %d ints", level, N_comp, rep_i),
    progress_log = progress_log, ladder = EX2_LADDER
  )
  fits$ints$name <- "ints"

  # Map per-unit model outputs (indexed by permuted untreated columns) back to the
  # original unit order, so downstream storage and the plots align with test_ys /
  # gen$groups. Scalar outputs (delta, loc_cor_pval, pred_perc) are permutation-invariant.
  for (m in c("no_ints", "ints")) {
    fits[[m]]$cor_sq <- unpermute_untreated(fits[[m]]$cor_sq, perm)
    fits[[m]]$abs_cors_err <- unpermute_untreated(fits[[m]]$abs_cors_err, perm)
  }

  res <- fits |>
    map(function(pfit) {
      res <- list()

      # Posterior-predictive interval coverage, per (time, unit) as ex1 does. Compared against
      # fit_ys, which is in the same anchor-permuted column order as y_pred.
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

      # The estimated error scale, on the data's own scale. In absolute mode err_sd is a single
      # shared eta and draws("tau") returns eta / sigma[n], so multiply back by this arm's sigma[1].
      # Recorded with its prior tail so a prior-data conflict on the error scale is visible.
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

      cor_sq <- pfit$cor_sq
      res[paste0("cor_sq_", seq_along(cor_sq))] <- cor_sq

      acors_err <- pfit$abs_cors_err
      res[paste0("acor_err_", seq_along(acors_err))] <- acors_err

      pred_mad <- pfit$mean_abs_diffs
      res$pred_mad <- pred_mad

      # Convergence diagnostics recorded for every fit, so the run can be audited rather than only
      # its failures. rhat_M is the one that certifies delta: the likelihood sees (Lambda, Phi) only
      # through M = Lambda_Phi, with the intercepts entering additively and the factor means folding
      # into Phi.
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
  # Worker-called functions must be exported explicitly (foreach only auto-exports
  # locals like K_latent); posterior is attached for sample_model()'s unqualified
  # extract_variable_array() call, ggplot2 for the per-condition figures.
  exp_vars <- c(
    "run_sim_intercepts",
    "worker_progress", "sample_model", "ife_mod", "plot_intercepts_fits",
    "anchor_order", "unpermute_untreated",
    "pathfinder_inits", "draw_to_init", "PF_PARAM_BASES",
    "fit_with_escalation", "escalation_ladder", "ESCALATE_MAX", "EX2_LADDER",
    "EX2_ITER", "EX2_WARM"
  )
  exp_packages <- c("cmdstanr", "posterior", "ggplot2", "dplyr")

  # Flatten the sim x N_comp x rep design into a single (non-nested) foreach, so
  # %dorng% gives each task a reproducible RNG stream invariant to worker count.
  # Invariant: keep this a single, non-nested loop.
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
  # Absolute log path, so workers write it where the master expects regardless of
  # their working directory.
  progress_log <- file.path(getwd(), "progress.log")
  t0 <- Sys.time()

  # Each task writes its own row as it finishes, and a rerun picks up the rows that already exist.
  # A task that errors records failed = TRUE and the study continues. The directory is keyed to
  # mode, seed and grid size; delete it to force a full recompute.
  ckpt_dir <- file.path(getwd(),
    sprintf("ckpt_ints_%s_seed%d_n%d", STUDY_MODE, seed, nrow(grid)))
  dir.create(ckpt_dir, showWarnings = FALSE)
  ckpt_check_fingerprint(ckpt_dir, c("no_ints", "ints"), "ex2_config.r")
  n_resume <- length(list.files(ckpt_dir, pattern = "^task_.*\\.rds$"))
  if (n_resume > 0) {
    cat(sprintf("  resuming: %d of %d tasks already checkpointed in %s\n\n",
      n_resume, nrow(grid), basename(ckpt_dir)))
  } else {
    cat("", file = progress_log) # truncate: fresh per-worker progress log per run
  }

  study_res <-
    foreach(
      rep_i = grid$rep, N_comp = grid$N_comp, level = grid$level, task_i = seq_len(nrow(grid)),
      # bind_rows rather than rbind: a failed task contributes a short row, and rbind would error on
      # the column mismatch instead of recording the failure.
      .combine = function(...) dplyr::bind_rows(...),
      .export = exp_vars, .packages = exp_packages,
      .options.RNG = seed
    ) %dorng% {
      # The full condition triple goes in the filename, so a checkpoint can never be reused for a
      # different cell of the grid even if the grid ordering were to change.
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
  N_comps = c(2, 3),
  levels = DGP_LEVELS,
  K_latent = K_LATENT,
  seed = 52918,
  plot_iters = 50
)

stopCluster(cl)

# Save the raw study results; numeric summaries and plots are produced by
# ex2_sim_study_summary.r (run it to view the results).
# Mode-keyed output, so a quick fast-mode check cannot silently destroy the full study's results.
out_file <- if (STUDY_MODE == "fast") "sim_study_ints_fast.RData" else "sim_study_ints.RData"
save(sim_study_ints, file = out_file)
cat(sprintf("Results saved to %s -- run `Rscript ex2_sim_study_summary.r %s` to summarize.\n",
  out_file, out_file))
