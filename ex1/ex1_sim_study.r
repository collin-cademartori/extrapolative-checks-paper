## This file runs the simulation study for the nonstationary example, fitting the
## nonstationary model and two stationary models (with weak and strong priors on the 
## iid error scale) to samples from the prior predictive distribution of the
## nonstationary model.

library(foreach)
library(doParallel)
library(doRNG)
library(purrr)

cl <- makeCluster(round(detectCores()/2) - 1, outfile = "")
registerDoParallel(cl)

# Quietly pre-attach the packages the workers need. `outfile = ""` surfaces all
# worker output, so otherwise each worker's package startup banners would clutter
# the console on first task dispatch. Pre-attaching them here silences that -- an
# already-attached package is re-attached with no banner.
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(forcats)
  library(dplyr)
  library(ggplot2)
  library(doRNG)
})))

source("../sample_model.r")
source("../plotting.r")

# Per-worker progress: each worker keeps a running count of the tasks it has
# completed (persisted in its own session via the global `.worker_done`) and
# prints a tagged line, so the workers' pace can be compared at a glance. No
# inter-process communication -- each line is self-contained.
worker_progress <- function(label) {
  n <- get0(".worker_done", envir = globalenv(), ifnotfound = 0L) + 1L
  assign(".worker_done", n, envir = globalenv())
  message(sprintf("[worker %d] %3d done | %s", Sys.getpid(), n, label))
}

run_sim_stat <- function(test_data, i, K_latent, post_check = FALSE) {
  test_ys <- test_data$ys[i, , ]
  N_units <- ncol(test_ys)
  T_times <- nrow(test_ys)
  overall_scales_stat <- apply(test_ys, 2, sd)
  overall_scales_nonstat <- apply(test_ys, 2, function(y) sd(diff(y)))

  # Draw all three Stan seeds up front, before any sample_model() call. cmdstanr's
  # $sample() advances R's RNG by a worker-count-dependent amount, so a seed drawn
  # after a fit would be misaligned across worker counts, breaking reproducibility.
  fit_seeds <- sample.int(.Machine$integer.max, 3)

  fits <- list()

  fits$nonstat <- sample_model(
    overall_scales = overall_scales_nonstat, err_scale = 0,
    err_scale_mean = 2,
    err_scale_sd = 2,
    data = test_ys,
    autocor_a = 8, autocor_b = 2,
    nonstationary = TRUE, num_treated = 5,
    type = "posterior", K_latent = K_latent,
    iter = 400,
    n_chains = 1, seed = fit_seeds[1]
  )
  fits$nonstat$name <- "nonstat"

  fits$stat2 <- sample_model(
    overall_scales = overall_scales_stat, err_scale = 0,
    err_scale_mean = 0.1,
    err_scale_sd = 0.1,
    data = test_ys,
    autocor_a = 97, autocor_b = 3,
    nonstationary = FALSE, num_treated = 5,
    type = "posterior", K_latent = K_latent,
    iter = 400,
    n_chains = 1, seed = fit_seeds[2]
  )
  fits$stat2$name <- "stat2"

  fits$stat1 <- sample_model(
    overall_scales = overall_scales_stat, err_scale = 0.05,
    data = test_ys,
    autocor_a = 97, autocor_b = 3,
    nonstationary = FALSE, num_treated = 5,
    type = "posterior", K_latent = K_latent,
    iter = 400,
    n_chains = 1, seed = fit_seeds[3]
  )
  fits$stat1$name <- "stat1"

  res <- fits |> map(function(pfit) {

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
    
    res$time_cor_pval <- pfit$time_cor_pval

    absz <- abs(pfit$effect_means / pfit$effect_sds)
    res[paste0("absz_", seq_along(absz))] <- absz

    pmean <- pfit$effect_means
    res[paste0("mean_", seq_along(pmean))] <- pmean

    psds <- pfit$effect_sds
    res[paste0("sd_", seq_along(psds))] <- psds

    pred_mad <- pfit$mean_abs_diffs
    res$pred_mad <- pred_mad

    return(res)

  }) |> list_flatten()

  pns_means <- apply(fits$nonstat$y_means, c(2,3), mean)
  p2_means <- apply(fits$stat2$y_means, c(2,3), mean)
  p1_means <- apply(fits$stat1$y_means, c(2,3), mean)

  fit_plot <- plot_post_fits_all(test_ys, pns_means, p2_means, p1_means)
  ggsave(fit_plot, file=paste0("../figs/sim_stat_figs/post_fit_plot_", i, ".png"), create.dir = TRUE)

  if (post_check) {
    check_plot <- plot_data_matrix_post(test_ys, fits$stat2$y_pred)
    ggsave(
        check_plot, file=paste0("../figs/sim_stat_figs/check_plot_", i, ".pdf"),
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

  test_data <- sample_model(overall_scales = rep(1, 8), err_scale = 3,
                            autocor_a = 8, autocor_b = 2,
                            nonstationary = TRUE, num_treated = 0,
                            type = "prior_pred", K_latent = K_latent,
                            iter = 2 * reps, seed = pp_seed)

  exp_vars <- c('run_sim_stat', 'worker_progress', 'sample_model', 'ife_mod', 'plot_post_fits_all', 'plot_data_matrix_post')
  exp_packages <- c('cmdstanr', 'posterior', 'forcats', 'dplyr', 'ggplot2')
  cat(sprintf(
    paste0("\n=== Example 1 simulation study (nonstationary) ===\n",
           "  reps    : %d\n",
           "  tasks   : %d  (3 model fits each)\n",
           "  workers : %d   seed: %d\n\n"),
    reps, reps, getDoParWorkers(), seed))
  t0 <- Sys.time()

  # Single (non-nested) foreach, so %dorng% gives each task a reproducible RNG
  # substream from `seed`, invariant to the number of workers.
  study_res <-
    foreach(
      s = study_units, iter = seq(reps),
      .combine = 'rbind', .export = exp_vars, .packages = exp_packages,
      .options.RNG = seed
    ) %dorng% {
      unit_res <- as.data.frame(run_sim_stat(test_data, s, K_latent, post_check))
      worker_progress(sprintf("iteration %d (unit %d)", iter, s))
      unit_res
    }

  cat(sprintf("--- study complete: %d tasks in %.1f min ---\n",
              reps, as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  return(study_res)
}

study_reps <- 1000
sim_study_stat <- run_sim_study_stat(K_latent = 4, study_reps, seed = 40318)

stopCluster(cl)

# Save the raw study results; numeric summaries and plots are produced by
# ex1_sim_study_summary.r (run it to view the results).
save(sim_study_stat, file="sim_study_ns.RData")
cat("Results saved to sim_study_ns.RData -- run ex1_sim_study_summary.r to summarize.\n")
