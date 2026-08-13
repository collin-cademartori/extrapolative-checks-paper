library(cmdstanr)
library(posterior)
library(ggplot2)
library(forcats)
library(dplyr)

ife_mod <- cmdstan_model(stan_file = "../ife_named.stan")

sample_model <- function(
    N_units = 8, T_times = 20, K_latent = 4,
    data = NULL, overall_scales = NULL, err_scale = 0.05,
    err_scale_mean = 0, err_scale_sd = 0,
    autocor_a, autocor_b, nonstationary, int_scale = 1, int_loc = 0, include_ints = FALSE, include_factor_means = FALSE,
    num_treated, type = "prior_pred", iter = 1000, quiet = TRUE, ad = 0.98,
    seed = NULL, n_chains = 4, log_file = NULL, log_label = NULL) {
  stopifnot(type %in% c("prior_pred", "posterior"))
  stopifnot(0 < autocor_a)
  stopifnot(0 < autocor_b)
  stopifnot(K_latent < N_units)
  stopifnot(num_treated >= 0 && num_treated < T_times)
  stopifnot(err_scale > 0 || (err_scale_mean > 0 && err_scale_sd > 0))

  # Draw the shuffle index *before* $sample(): cmdstanr's $sample() advances R's
  # global RNG, so drawing it afterward would make the draw ordering
  # irreproducible. Invariant: this must precede the fit.
  sample_index <- sample.int(iter, size = iter)

  stat_data <- list(
    M_units = N_units,
    T_times = T_times,
    K_latent = K_latent,
    Y = if (is.null(data)) matrix(0, nrow = T_times, ncol = N_units) else data,
    a_rho = autocor_a,
    b_rho = autocor_b,
    tau_val = err_scale,
    m_tau = err_scale_mean,
    s_tau = err_scale_sd,
    sigma_data =
      if (is.null(overall_scales)) rep(0, N_units) else overall_scales,
    fit_overall_scales = if (type == "prior_pred") 0 else 1,
    nonstationary = nonstationary,
    unit_intercepts = include_ints,
    factor_means = include_factor_means,
    sample_posterior = (type == "posterior"),
    num_treated = num_treated,
    gamma_scale = int_scale,
    gamma_loc = int_loc
  )

  model_sample <- ife_mod$sample(
    data = stat_data,
    parallel_chains = n_chains,
    chains = n_chains,
    iter_warmup = iter,
    iter_sampling = iter,
    adapt_delta = ad,
    refresh = if (quiet) 0 else 100,
    show_exceptions = !quiet,
    show_messages = !quiet,
    seed = seed
  )

  # Optionally note this fit's sampler diagnostics (divergences, max-treedepth hits,
  # min E-BFMI, and the bulk ESS of the first treatment-effect element as a cheap
  # mixing canary) to a shared log, tagged with the caller's label, so Stan warnings
  # can be cross-correlated with the simulation rep. Only problematic fits are logged,
  # so a clean run leaves the log to the per-rep progress lines.
  if (!is.null(log_file)) {
    ds <- model_sample$diagnostic_summary(quiet = TRUE)
    n_div <- sum(ds$num_divergent)
    n_tree <- sum(ds$num_max_treedepth)
    ebfmi_min <- suppressWarnings(min(ds$ebfmi))
    ess_delta1 <- if (num_treated > 0) {
      model_sample$summary("delta", "ess_bulk")$ess_bulk[1]
    } else {
      NA_real_
    }
    if (n_div > 0 || n_tree > 0 || (is.finite(ebfmi_min) && ebfmi_min < 0.3) ||
      (is.finite(ess_delta1) && ess_delta1 < 100)) {
      cat(sprintf(
        "[%s] %s  STAN div=%d treedepth=%d ebfmi_min=%.2f ess_delta1=%.0f\n",
        format(Sys.time(), "%H:%M"), log_label, n_div, n_tree, ebfmi_min, ess_delta1
      ), file = log_file, append = TRUE)
    }
  }

  if (type == "prior_pred") {
    ys_prior_all <-
      extract_variable_array(model_sample$draws("Y_prior"), "Y_prior")
    ys_prior <- ys_prior_all[sample_index, 1, , ]
    ys_latent_all <-
      extract_variable_array(model_sample$draws("Y_latent"), "Y_latent")
    ys_latent <- ys_latent_all[sample_index, 1, , ]

    return(list(
      ys = ys_prior,
      ys_latent = ys_latent
    ))
  } else if (type == "posterior") {
    y_means_all <-
      extract_variable_array(model_sample$draws("Y_latent"), "Y_latent")
    y_means_post <- y_means_all[sample_index, 1, , ]

    y_pred_all <-
      extract_variable_array(model_sample$draws("Y_pred"), "Y_pred")
    y_pred_post <- y_pred_all[sample_index, 1, , ]

    effects <-
      extract_variable_array(model_sample$draws(), "delta")[, 1, ]
    effect_means <- colMeans(effects)
    effect_sds <- apply(effects, 2, sd)

    err_scale <- as.numeric(model_sample$draws("tau"))
    mad <- mean(as.numeric(model_sample$draws("mean_abs_diffs")))

    cor_sq <- extract_variable_array(model_sample$draws(), "cor_sq")[, 1, ]
    cor_sq_mean <- colMeans(cor_sq)

    abs_cor_pred <- as.numeric(model_sample$draws("time_cor_pred"))
    abs_cor_data <- abs(cor(data[, 2], seq(nrow(data))))
    time_cor_pval <- mean(abs_cor_pred > abs_cor_data)

    # Statistic S2 (paper Section 5) on the pre-treatment window: the correlation
    # across untreated units between their correlation with the treated unit and
    # their mean level. The p-value compares S2 on the predictive replicates to
    # S2 on the observed data.
    loc_cor_pred <- as.numeric(model_sample$draws("loc_cor_pred"))
    pre_times <- seq_len(nrow(data) - num_treated)
    # Untreated units over the pre-treatment window; drop = FALSE keeps this a
    # matrix so colMeans() works (relevant only in the degenerate single-untreated
    # case, N_units == 2, which does not arise in the examples).
    untreated_pre <- data[pre_times, -1, drop = FALSE]
    cor_with_treated <- as.numeric(cor(untreated_pre, data[pre_times, 1]))
    unit_location <- colMeans(untreated_pre)
    loc_cor_data <- cor(cor_with_treated, unit_location)
    loc_cor_pval <- mean(loc_cor_pred > loc_cor_data)

    y_cor <- abs(cor(data)[1, 2:ncol(data)])
    cor_err_mean <- rowMeans(abs(y_cor - t(cor_sq)))

    return(list(
      y_means = y_means_post,
      y_pred = y_pred_post,
      effect_means = effect_means,
      effect_sds = effect_sds,
      mean_abs_diffs = mad,
      cor_sq = cor_sq_mean,
      abs_cors_err = cor_err_mean,
      err_scale = err_scale,
      time_cor_pval = time_cor_pval,
      loc_cor_pval = loc_cor_pval
    ))
  }
}
