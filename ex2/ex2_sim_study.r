## This file runs the simulation study for the intercepts example, fitting the
## fitting the models with and without intercepts to data generated from an
## adversarial process whereby some units are simulated to have different
## long run means than the treated but spuriously correlation in the pre-
## treatment period only.

library(foreach)
library(doParallel)
library(doRNG)
library(purrr)

# Worker count: pass as the first CLI arg (e.g. `Rscript ex2_sim_study.r 4`),
# otherwise a conservative default.
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
  library(ggplot2)
  library(doRNG)
  library(purrr)
})))

source("../sample_model.r")
source("../pathfinder_init.r")
source("../plotting.r")

ruv <- function(d) {
  v <- rnorm(d)
  uv <- v / sqrt(sum(v * v))
  return(uv)
}

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

# Data-driven anchor ordering for the triangular (Cholesky) loadings. Keeps the treated
# unit (column 1) first, then greedily picks the untreated columns most orthogonal to
# those already chosen so the leading K x K loading block is
# full rank -- giving every factor a distinct anchor and avoiding near-zero
# Cholesky diagonals. This is a reparameterization only: the
# fitted means, treatment effect, and (permutation-invariant) check statistics do not
# depend on the choice; only the per-unit outputs need mapping back (see below).
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

# Adversarial DGP for the intercepts example (paper Section 5): untreated units
# split into three groups differing in their pre-treatment correlation with the
# treated unit and their long-run average, so location and correlation are
# entangled in a way the unit-intercepts model wrongly assumes independent.
sim_model_intercepts <- function(
    N_unc = 2, N_comp_true = 2, N_comp_spur = 2, T_times = 30, T_treated = 5,
    K_unc = 1, sim = 0.9) {
  N_units <- 1 + N_comp_true + N_comp_spur + N_unc
  K_gen <- 2 + K_unc

  f_treat <- 6 + arima.sim(model = list(ar = 0.9), n = T_times)
  f_treat_sd <- 1.9

  # f_alt matches the treated factor pre-treatment, then diverges downward over the
  # treatment window -- the driver of the "spurious" comparators.
  f_alt <- (f_treat - 6) +
    c(rep(0, T_times - T_treated), rep(-f_treat_sd, T_treated))

  # Reject until the "uncorrelated" factors are genuinely uncorrelated with the
  # treated factor.
  f_unc <- matrix(nrow = K_unc, ncol = T_times)
  cor_unc <- Inf
  while (cor_unc > 0.01) {
    for (k in 1:K_unc) {
      f_unc[k, ] <- rnorm(1, 0, 2) + arima.sim(model = list(ar = 0.9), n = T_times)
    }
    cor_unc <- max(abs(cor(t(f_unc), t(t(f_treat)))))
  }

  facs <- rbind(f_treat, f_alt, f_unc)
  loads <- matrix(nrow = N_units, ncol = K_gen)
  loads[1, ] <- c(1, rep(0, K_gen - 1))
  # True comparators load on the treated factor (genuine correlation throughout).
  # seq_len (not 1:N) so a group size of 0 yields an empty loop, not 1:0 = c(1, 0).
  for (n in seq_len(N_comp_true)) {
    loads[1 + n, ] <- c(sqrt(sim), 0, sqrt(1 - sim) * ruv(K_gen - 2))
  }
  # Spurious comparators load on f_alt (pre-treatment correlation only).
  for (n in seq_len(N_comp_spur)) {
    loads[1 + N_comp_true + n, ] <- c(0, sqrt(sim), sqrt(1 - sim) * ruv(K_gen - 2))
  }
  # Uncorrelated units (fill to a fixed total; carry no shared-factor signal).
  for (n in seq_len(N_unc)) {
    loads[1 + N_comp_true + N_comp_spur + n, ] <- c(0, 0, ruv(K_gen - 2))
  }

  # Treated and true comparators sit high (via f_treat), spurious low (via f_alt), so
  # correlation with the treated does not determine location -- contrary to the
  # unit-intercepts model's assumption. Y is returned T x N (the fit orientation).
  lat <- loads %*% facs
  Y <- t(lat + rnorm(nrow(lat) * ncol(lat), sd = 0.1 * max(apply(lat, 1, sd))))

  # Ground-truth group of each unit (column), in generating order: the treated unit,
  # the true comparators, the spurious comparators, then the uncorrelated units. The
  # plotting code uses this to color the post-treatment segments by known truth.
  groups <- c(
    "treated",
    rep("true", N_comp_true),
    rep("spurious", N_comp_spur),
    rep("uncorrelated", N_unc)
  )

  return(list(Y = Y, groups = groups))
}

run_sim_intercepts <- function(N_comp, sim, K_latent = 3, rep_i = NA, plot_iters = 0,
                               progress_log = NULL) {
  # Fixed total of 8 units; N_comp spurious comparators trade off against the
  # uncorrelated fillers (1 treated + 2 true + N_comp spurious + N_unc = 8), so the
  # swept quantity is the *share* of units spuriously correlated with the treated.
  N_unc <- 5 - N_comp
  gen <- sim_model_intercepts(
    N_unc = N_unc, N_comp_true = 2, N_comp_spur = N_comp, K_unc = 1,
    sim = sim, T_times = 30
  )
  test_ys <- gen$Y

  N_units <- ncol(test_ys)
  T_times <- nrow(test_ys)

  # Fit on anchor-ordered columns (treated stays first) so the leading K x K loading
  # block is full rank; per-unit outputs are mapped back to the original order after the
  # fits. The plots below use the original test_ys.
  perm <- anchor_order(test_ys, K_latent)
  fit_ys <- test_ys[, perm]

  # Draw both Stan seeds up front, before any sample_model() call: cmdstanr's
  # $sample() advances R's global RNG, so a seed drawn after a fit would not be
  # reproducible. Invariant: never derive a seed after a fit.
  fit_seeds <- sample.int(.Machine$integer.max, 2)

  fits <- list()

  overall_scales <- apply(fit_ys, 2, \(x) sqrt(mean(x ^ 2)))
  fits$no_ints <- sample_model(
    N_units = ncol(fit_ys), T_times = nrow(fit_ys), K_latent = K_latent,
    overall_scales = overall_scales,
    err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.05,
    data = fit_ys,
    autocor_a = 90, autocor_b = 10,
    nonstationary = FALSE, num_treated = 5,
    include_factor_means = TRUE,
    fit_scales = 0, alpha_diag = 0, pathfinder_init = TRUE,
    type = "posterior", quiet = TRUE, ad = 0.8,
    iter = 500, iter_warm = 500, max_treedepth = 12,
    n_chains = 4, seed = fit_seeds[1],
    log_file = progress_log,
    log_label = sprintf("sim %.2g num_comp %d rep %d no_ints", sim, N_comp, rep_i)
  )
  fits$no_ints$name <- "no_ints"

  overall_sds <- apply(fit_ys, 2, sd)
  fits$ints <- sample_model(
    N_units = ncol(fit_ys), T_times = nrow(fit_ys), K_latent = K_latent,
    overall_scales = overall_sds,
    err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.05,
    data = fit_ys,
    autocor_a = 90, autocor_b = 10,
    nonstationary = FALSE, num_treated = 5,
    include_ints = TRUE, int_scale = 3, int_loc = 4,
    fit_scales = 0, alpha_diag = 0, pathfinder_init = TRUE,
    type = "posterior", quiet = TRUE, ad = 0.8,
    iter = 500, iter_warm = 500, max_treedepth = 12,
    n_chains = 4, seed = fit_seeds[2],
    log_file = progress_log,
    log_label = sprintf("sim %.2g num_comp %d rep %d ints", sim, N_comp, rep_i)
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

      cor_sq <- pfit$cor_sq
      res[paste0("cor_sq_", seq_along(cor_sq))] <- cor_sq

      acors_err <- pfit$abs_cors_err
      res[paste0("acor_err_", seq_along(acors_err))] <- acors_err

      pred_mad <- pfit$mean_abs_diffs
      res$pred_mad <- pred_mad

      return(res)
    }) |>
    list_flatten()

  if (!is.na(rep_i) && rep_i <= plot_iters) {
    for (m in c("no_ints", "ints")) {
      int_plot <- plot_intercepts_fits(
        test_ys,
        cor_sq = fits[[m]]$cor_sq,
        groups = gen$groups, num_treated = 5
      )
      ggsave(
        int_plot,
        file = sprintf("../figs/sim_int_figs/int_fit_%s_sim%.2g_nc%d_rep%d.png", m, sim, N_comp, rep_i),
        width = 5, height = 4, create.dir = TRUE
      )
    }
  }

  return(res)
}

run_sim_study_intercepts <- function(K_latent = 3, reps, N_comps, sims, seed, plot_iters = 3) {
  # Worker-called functions must be exported explicitly (foreach only auto-exports
  # locals like K_latent); posterior is attached for sample_model()'s unqualified
  # extract_variable_array() call, ggplot2 for the per-condition figures.
  exp_vars <- c(
    "run_sim_intercepts", "sim_model_intercepts", "ruv",
    "worker_progress", "sample_model", "ife_mod", "plot_intercepts_fits",
    "anchor_order", "unpermute_untreated",
    "pathfinder_inits", "draw_to_init", "PF_PARAM_BASES"
  )
  exp_packages <- c("cmdstanr", "posterior", "ggplot2")

  # Flatten the sim x N_comp x rep design into a single (non-nested) foreach, so
  # %dorng% gives each task a reproducible RNG stream invariant to worker count.
  # Invariant: keep this a single, non-nested loop.
  grid <- expand.grid(rep = seq_len(reps), N_comp = N_comps, sim = sims)

  cat(sprintf(
    paste0(
      "\n=== Example 2 simulation study ===\n",
      "  conditions : sim {%s} x num_comp {%s}  (%d)\n",
      "  reps/cond  : %d\n",
      "  tasks      : %d  (2 model fits each)\n",
      "  workers    : %d   seed: %d\n\n"
    ),
    paste(sims, collapse = ", "), paste(N_comps, collapse = ", "),
    length(sims) * length(N_comps), reps, nrow(grid),
    getDoParWorkers(), seed
  ))
  # Absolute log path, so workers write it where the master expects regardless of
  # their working directory.
  progress_log <- file.path(getwd(), "progress.log")
  cat("", file = progress_log) # truncate: fresh per-worker progress log per run
  t0 <- Sys.time()

  study_res <-
    foreach(
      rep_i = grid$rep, N_comp = grid$N_comp, sim = grid$sim,
      .combine = "rbind", .export = exp_vars, .packages = exp_packages,
      .options.RNG = seed
    ) %dorng% {
      unit_res <- run_sim_intercepts(
        N_comp = N_comp, sim = sim, K_latent = K_latent,
        rep_i = rep_i, plot_iters = plot_iters, progress_log = progress_log
      )
      worker_progress(sprintf("sim %.2g  num_comp %d  rep %d", sim, N_comp, rep_i), logfile = progress_log)
      as.data.frame(c(unit_res, list(sim = sim, num_comp = N_comp)))
    }

  cat(sprintf(
    "--- study complete: %d tasks in %.1f min ---\n",
    nrow(grid), as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ))

  return(study_res)
}

sim_study_ints <- run_sim_study_intercepts(
  reps = 200,
  N_comps = c(2, 3),
  sims = c(0.7, 0.9),
  K_latent = 3,
  seed = 52918,
  plot_iters = 50
)

stopCluster(cl)

# Save the raw study results; numeric summaries and plots are produced by
# ex2_sim_study_summary.r (run it to view the results).
save(sim_study_ints, file = "sim_study_ints.RData")
cat("Results saved to sim_study_ints.RData -- run ex2_sim_study_summary.r to summarize.\n")
