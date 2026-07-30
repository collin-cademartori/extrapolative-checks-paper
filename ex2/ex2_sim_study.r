## This file runs the simulation study for the intercepts example, fitting the
## fitting the models with and without intercepts to data generated from an
## adversarial process whereby some units are simulated to have different
## long run means than the treated but spuriously correlation in the pre-
## treatment period only.

library(foreach)
library(doParallel)
library(doRNG)
library(purrr)

cl <- makeCluster(round(detectCores()/2) - 1, outfile = "")
registerDoParallel(cl)

# Pre-attach the workers' packages quietly, so their startup banners don't clutter
# the console (outfile = "" surfaces all worker output).
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(doRNG)
})))

source("../sample_model.r")
source("../plotting.r")

ruv <- function(d) {
  v <- rnorm(d)
  uv <- v / sqrt(sum(v * v))
  return(uv)
}

# Per-worker progress: each worker prints a running count of the tasks it has
# completed, so the workers' pace can be compared at a glance.
worker_progress <- function(label) {
  n <- get0(".worker_done", envir = globalenv(), ifnotfound = 0L) + 1L
  assign(".worker_done", n, envir = globalenv())
  message(sprintf("[worker %d] %3d done | %s", Sys.getpid(), n, label))
}

# Adversarial DGP for the intercepts example (paper Section 5): untreated units
# split into three groups differing in their pre-treatment correlation with the
# treated unit and their long-run average, so location and correlation are
# entangled in a way the unit-intercepts model wrongly assumes independent.
sim_model_intercepts <- function(
  N_unc = 5, N_comp_true = 2, N_comp_spur = 2, T_times = 20, T_treated = 5, K_unc = 3, sim = 0.9) {

    f_treat <- arima.sim(model = list(ar = 0.96), n = T_times)

    # f_alt matches the treated factor pre-treatment, then diverges downward over
    # the treatment window -- the driver of the "spurious" comparators.
    f_alt <- f_treat +
      c(rep(0, T_times - T_treated), (-1 / 2) * seq(T_treated))

    # Reject until the "uncorrelated" factors are genuinely uncorrelated (< 0.4)
    # with the treated factor.
    f_unc <- matrix(nrow = K_unc, ncol = T_times)
    cor_unc <- Inf
    while(cor_unc > 0.4) {
      for(k in 1:K_unc) {
        f_unc[k, ] <- arima.sim(model = list(ar = 0.96), n = T_times)
      }
      cor_unc <- max(abs(cor(t(f_unc), t(t(f_treat)))))
    }

    facs <- rbind(f_treat, f_alt, f_unc)

    N_units <- 1 + N_comp_true + N_comp_spur + N_unc
    K_gen <- 2 + K_unc
    loads <- matrix(nrow = N_units, ncol = K_gen)

    loads[1, ] <- c(1, rep(0, K_gen - 1))

    # True comparators load on the treated factor (genuine correlation throughout).
    for(n in 1:N_comp_true) {
      loads[1 + n, ] <- c(sqrt(0.9), 0, sqrt(1 - 0.9) * ruv(K_gen - 2))
    }

    # Spurious comparators load on f_alt (pre-treatment correlation only).
    for(n in 1:N_comp_spur) {
      loads[1 + N_comp_true + n, ] <- c(0, sqrt(sim), sqrt(1 - sim) * ruv(K_gen - 2))
    }

    # Uncorrelated units.
    for(n in 1:N_unc) {
      loads[1 + N_comp_true + N_comp_spur + n, ] <- c(0, 0 , ruv(K_gen - 2))
    }

    lat <- loads %*% facs
    # Intercepts are the units' long-run averages: the treated and true comparators
    # sit high, the spurious comparators low -- so correlation with the treated does
    # not determine location, contrary to the unit-intercepts model's assumption.
    intercepts <- c(
      5, rep(5, N_comp_true),
      rep(1, N_comp_spur),
      rnorm(N_unc, mean = 2.5, sd = 1)
    )
    
    Y <- lat + intercepts + rnorm(N_units * T_times, sd = 0.02)

    return(t(Y))

}

run_sim_intercepts <- function(N_comp, sim, K_latent = 5) {
  test_ys <- sim_model_intercepts(N_unc = 7 - N_comp, N_comp_spur = N_comp, K_unc = 3, sim = sim)

  N_units <- ncol(test_ys)
  T_times <- nrow(test_ys)
  overall_scales <- apply(test_ys, 2, sd)

  # Draw both Stan seeds up front, before any sample_model() call: cmdstanr's
  # $sample() advances R's global RNG, so a seed drawn after a fit would not be
  # reproducible. Invariant: never derive a seed after a fit.
  fit_seeds <- sample.int(.Machine$integer.max, 2)

  fits <- list()

  fits$no_ints <- sample_model(N_units = 10, T_times = 20, K_latent = K_latent,
                            overall_scales = overall_scales, err_scale = 0.2,
                            data = test_ys,
                            autocor_a = 90, autocor_b = 10,
                            nonstationary = FALSE, num_treated = 5,
                            type = "posterior", quiet = TRUE, ad = 0.8, iter = 500,
                            n_chains = 1, seed = fit_seeds[1])
  fits$no_ints$name <- "no_ints"

  fits$ints <- sample_model(N_units = 10, T_times = 20, K_latent = K_latent,
                            overall_scales = overall_scales, err_scale = 0.2,
                            data = test_ys,
                            autocor_a = 90, autocor_b = 10,
                            nonstationary = FALSE, num_treated = 5,
                            include_ints = TRUE, int_scale = 10,
                            type = "posterior", quiet = TRUE, iter = 500,
                            n_chains = 1, seed = fit_seeds[2])
  fits$ints$name <- "ints"

  res <- fits |> map(function(pfit) {

    res <- list()

    stat_y_pred <- pfit$y_pred
    pred_inc <- matrix(NA, nrow = T_times, ncol = N_units)
    for (n in 1:N_units) {
      for (t in 1:T_times) {
        y_bounds <- quantile(stat_y_pred, c(0.005, 0.995))
        pred_inc[t, n] <-
          (test_ys[t, n] >= y_bounds[1]) &&
          (test_ys[t, n] <= y_bounds[2])
      }
    }
    res$pred_perc <- mean(pred_inc)

    res$loc_cor_pval <- pfit$loc_cor_pval

    absz <- abs(pfit$effect_means / pfit$effect_sds)
    res[paste0("absz_", seq_along(absz))] <- absz

    pmean <- pfit$effect_means
    res[paste0("mean_", seq_along(pmean))] <- pmean

    psds <- pfit$effect_sds
    res[paste0("sd_", seq_along(psds))] <- psds

    acors <- pfit$abs_cors
    res[paste0("acor_", seq_along(acors))] <- acors

    acors_err <- pfit$abs_cors_err
    res[paste0("acor_err_", seq_along(acors_err))] <- acors_err

    pred_mad <- pfit$mean_abs_diffs
    res$pred_mad <- pred_mad

    return(res)

  }) |> list_flatten()

  return(res)
}

run_sim_study_intercepts <- function(K_latent = 6, reps, N_comps, sims, seed) {
  # Worker-called functions must be exported explicitly (foreach only auto-exports
  # locals like K_latent); posterior is attached for sample_model()'s unqualified
  # extract_variable_array() call.
  exp_vars <- c('run_sim_intercepts', 'sim_model_intercepts', 'ruv',
                'worker_progress', 'sample_model', 'ife_mod')
  exp_packages <- c('cmdstanr', 'posterior')

  # Flatten the sim x N_comp x rep design into a single (non-nested) foreach, so
  # %dorng% gives each task a reproducible RNG stream invariant to worker count.
  # Invariant: keep this a single, non-nested loop.
  grid <- expand.grid(rep = seq_len(reps), N_comp = N_comps, sim = sims)

  cat(sprintf(
    paste0("\n=== Example 2 simulation study ===\n",
           "  conditions : sim {%s} x num_comp {%s}  (%d)\n",
           "  reps/cond  : %d\n",
           "  tasks      : %d  (2 model fits each)\n",
           "  workers    : %d   seed: %d\n\n"),
    paste(sims, collapse = ", "), paste(N_comps, collapse = ", "),
    length(sims) * length(N_comps), reps, nrow(grid),
    getDoParWorkers(), seed))
  t0 <- Sys.time()

  study_res <-
    foreach(
      rep_i = grid$rep, N_comp = grid$N_comp, sim = grid$sim,
      .combine = 'rbind', .export = exp_vars, .packages = exp_packages,
      .options.RNG = seed
    ) %dorng% {
      unit_res <- run_sim_intercepts(N_comp = N_comp, sim = sim, K_latent = K_latent)
      worker_progress(sprintf("sim %.2g  num_comp %d  rep %d", sim, N_comp, rep_i))
      as.data.frame(c(unit_res, list(sim = sim, num_comp = N_comp)))
    }

  cat(sprintf("--- study complete: %d tasks in %.1f min ---\n",
              nrow(grid), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  return(study_res)
}

sim_study_ints <- run_sim_study_intercepts(
  reps = 300,
  N_comps = c(2, 3),
  sims = c(0.7, 0.9),
  seed = 52918
)

stopCluster(cl)

# Save the raw study results; numeric summaries and plots are produced by
# ex2_sim_study_summary.r (run it to view the results).
save(sim_study_ints, file="sim_study_ints.RData")
cat("Results saved to sim_study_ints.RData -- run ex2_sim_study_summary.r to summarize.\n")
