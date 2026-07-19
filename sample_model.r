library(cmdstanr)
library(posterior)
library(ggplot2)
library(forcats)
library(dplyr)

ife_mod <- cmdstan_model(stan_file = "../ife.stan")

sample_model <- function(
  N_units = 8, T_times = 20, K_latent = 4,
  data = NULL, overall_scales = NULL, err_scale = 0.05,
  err_scale_mean = 0, err_scale_sd = 0,
  autocor_a, autocor_b, nonstationary, int_scale = 1, include_ints = FALSE,
  num_treated, type = "prior_pred", iter = 1000, quiet = TRUE, ad = 0.98,
  seed = NULL, n_chains = 4
) {
  stopifnot(type %in% c("prior_pred", "posterior"))
  stopifnot(0 < autocor_a)
  stopifnot(0 < autocor_b)  
  stopifnot(K_latent < N_units)
  stopifnot(num_treated >= 0 && num_treated < T_times)
  stopifnot(err_scale > 0 || (err_scale_mean > 0 && err_scale_sd > 0))

  full_size <- iter

  stat_data <- list(
    N_units = N_units,
    T_times = T_times,
    K_latent = K_latent,
    Y = if (is.null(data)) matrix(0, nrow = T_times, ncol = N_units) else data,
    autocor_a = autocor_a,
    autocor_b = autocor_b,
    err_scale_val = err_scale,
    err_scale_mean = err_scale_mean,
    err_scale_sd = err_scale_sd,
    overall_scales_0 =
      if (is.null(overall_scales)) rep(0, N_units) else overall_scales,
    fit_overall_scales = if (type == "prior_pred") 0 else 1,
    nonstationary = nonstationary,
    unit_intercepts = include_ints,
    sample_posterior = (type == "posterior"),
    num_treated = num_treated,
    intercept_scale = int_scale
  )

  model_sample <- ife_mod$sample(
    data = stat_data,
    parallel_chains = n_chains,
    chains = n_chains,
    iter_warmup = iter,
    iter_sampling = iter,
    adapt_delta = ad,
    refresh = if(quiet) 0 else 100,
    show_exceptions = !quiet,
    seed = seed
  )

  sample_index <- sample.int(full_size, size = iter)

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
      extract_variable_array(model_sample$draws(), "treatment_effects")[,1,]
    effect_means <- colMeans(effects)
    effect_sds <- apply(effects, 2, sd)

    err_scale <- as.numeric(model_sample$draws("err_scale"))
    writeLines(paste0("Estimated error scale: ", round(mean(err_scale), 3)))

    mad <- mean(as.numeric(model_sample$draws("mean_abs_diffs")))
    writeLines(paste0("Estimated error MAD: ", round(mean(mad), 3)))

    scale_mult <- extract_variable_array(model_sample$draws(), "overall_scales_param")[,1,]
    writeLines(paste0("Average scale multiplier: ", round(colMeans(scale_mult), 3)))

    abs_cors <- extract_variable_array(model_sample$draws(), "abs_cors")[,1,]
    abs_cors_mean <- colMeans(abs_cors)

    abs_cor_pred <- as.numeric(model_sample$draws("time_cor_pred"))
    abs_cor_data <- abs(cor(data[,1], seq(nrow(data))))
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
    cor_err_mean <- rowMeans(abs(y_cor - t(abs_cors)))

    return(list(
      y_means = y_means_post,
      y_pred = y_pred_post,
      effect_means = effect_means,
      effect_sds = effect_sds,
      mean_abs_diffs = mad,
      abs_cors = abs_cors_mean,
      abs_cors_err = cor_err_mean,
      err_scale = err_scale,
      time_cor_pval = time_cor_pval,
      loc_cor_pval = loc_cor_pval
    ))
  }

}
