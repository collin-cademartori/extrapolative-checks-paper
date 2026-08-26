library(cmdstanr)
library(posterior)
library(ggplot2)
library(forcats)
library(dplyr)

ife_mod <- cmdstan_model(stan_file = "../ife_named.stan")

sample_model <- function(
    N_units = 8, T_times = 20, K_latent = 4,
    data = NULL, overall_scales = NULL, fit_scales = NULL, err_scale = 0.05,
    err_scale_mean = 0, err_scale_sd = 0,
    autocor_a, autocor_b, alpha_diag = 0, nonstationary, int_scale = 1, int_loc = 0, include_ints = FALSE, include_factor_means = FALSE,
    num_treated, type = "prior_pred",
    iter = 1000, iter_warm = NULL, quiet = TRUE, 
    ad = 0.98, max_treedepth = 10, n_chains = 4, parallel_chains = 1,
    seed = NULL, log_file = NULL, log_label = NULL,
    return_draws = NULL, init = NULL,
    pathfinder_init = FALSE) {
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
    # fit_scales overrides the default (estimate sigma in posterior, fix in prior_pred):
    # fit_scales = 0 fixes sigma to the passed overall_scales. The loadings then carry any residual per-unit scale.
    fit_overall_scales = if (!is.null(fit_scales)) fit_scales else if (type == "prior_pred") 0 else 1,
    nonstationary = nonstationary,
    unit_intercepts = include_ints,
    factor_means = include_factor_means,
    sample_posterior = (type == "posterior"),
    num_treated = num_treated,
    gamma_scale = int_scale,
    gamma_loc = int_loc,
    alpha_diag = alpha_diag
  )

  # Optional Pathfinder warm start: seed every chain from a draw in the dominant lp mode, discarding modes which are vanishingly tiny compared to the dominant mode.
  init_arg <- if (isTRUE(pathfinder_init)) {
    pfi <- pathfinder_inits(ife_mod, stat_data, n_chains, seed = seed, quiet = quiet)
    if (is.null(pfi)) (if (is.null(init)) 2 else init) else pfi
  } else if (is.null(init)) 2 else init

  model_sample <- ife_mod$sample(
    data = stat_data,
    parallel_chains = parallel_chains,
    chains = n_chains,
    iter_warmup = ifelse(is.null(iter_warm), iter, iter_warm),
    iter_sampling = iter,
    adapt_delta = ad,
    max_treedepth = max_treedepth,
    init = init_arg,
    refresh = if (quiet) 0 else 100,
    show_exceptions = !quiet,
    show_messages = !quiet,
    seed = seed
  )

  # Sampler diagnostics (posterior fits only): divergences, max-treedepth hits, min
  # E-BFMI, plus a couple model-specific checks
  # Computed once and returned (see `sampler_diag` in the posterior list) for 
  # offline analysis, and  optionally written to a shared log tagged with log_label
  # for cross-correlating warnings with the simulation rep. 
  sampler_diag <- NULL
  if (type == "posterior") {
    ds <- model_sample$diagnostic_summary(quiet = TRUE)
    # rhat is reported three ways, because a single max over all parameters conflates two very
    # different situations. The loadings are only weakly identified: their marginals are bimodal (a
    # loading configuration can re-label or flip sign at essentially equal log density), so Lambda
    # entries mix slowly and rhat on them spikes intermittently -- with n_offmode = 0, no divergences,
    # and healthy ESS everywhere else. Taking |Lambda| removes it (rhat -> ~1.00, ESS x3-6), which is
    # what identifies it as a configuration ambiguity rather than a failure to converge.
    #   rhat_max       : all parameters (kept, so nothing is hidden)
    #   rhat_loadings  : Lambda / Phi_innovations -- expected looser when factors are weakly separated
    #   rhat_estimands : delta, tau, rho, sigma -- rotation-invariant; THESE must be clean
    #   rhat_cor_sq    : the squared loading correlation actually reported from the loadings. Being a
    #                    squared dot product it is invariant to the sign/label ambiguity above, so it
    #                    is the right convergence check for the loadings' scientific content.
    #   rhat_M         : max over the latent-mean matrix M = Lambda_Phi = Phi * Lambda'. THE key
    #                    diagnostic. In ife_named the likelihood depends on (Lambda, Phi) only through
    #                    M, and delta's prior involves only sigma[1]/omega_sq, so
    #                        delta  _||_  (Lambda, Phi)  |  M, sigma, tau.
    #                    Non-mixing confined to a fiber {(L,P): L P' = M} therefore CANNOT affect
    #                    delta, while any missed mode carrying a different delta must carry a
    #                    different M and so shows up here. Empirically it discriminates: at the study's
    #                    K = 4, Lambda R-hat of 1.038 came with rhat_M = 1.004 (fiber-internal, benign),
    #                    whereas an under-specified K = 3 fit gave Lambda 1.567 AND rhat_M 1.24 -- a
    #                    genuine warning that delta's own R-hat (1.002 there) would have missed.
    #                    NOTE this argument requires the error scale NOT to involve ||Lambda||.
    mixing <- tryCatch(
      {
        sv <- model_sample$metadata()$stan_variables
        pick <- function(nms) intersect(nms, sv)
        all_vars <- pick(c("Lambda", "Phi_innovations", "sigma_raw", "rho", "gamma_raw",
          "delta_raw", "omega_sq_param", "Phi_means_param", "tau_param"))
        ps <- model_sample$summary(all_vars, "rhat", "ess_bulk")
        grp <- function(bases) {
          v <- pick(bases); if (!length(v)) return(NA_real_)
          suppressWarnings(max(model_sample$summary(v, "rhat")$rhat, na.rm = TRUE))
        }
        list(
          rhat_max = suppressWarnings(max(ps$rhat, na.rm = TRUE)),
          rhat_loadings = grp(c("Lambda", "Phi_innovations")),
          rhat_estimands = grp(c("delta_raw", "tau_param", "rho", "sigma_raw",
            "gamma_raw", "omega_sq_param", "Phi_means_param")),
          rhat_cor_sq = grp("cor_sq"),
          rhat_M = grp("Lambda_Phi"),
          ess_delta1 = ps$ess_bulk[match("delta_raw[1]", ps$variable)]
        )
      },
      error = function(e) list(rhat_max = NA_real_, rhat_loadings = NA_real_,
        rhat_estimands = NA_real_, rhat_cor_sq = NA_real_, rhat_M = NA_real_, ess_delta1 = NA_real_)
    )
    # Minor-mode flag from the per-chain mean lp__: n_offmode counts every chain more than 5 below
    # the best chain; lp_gap_max is the worst gap.
    lp_gap_max <- NA_real_; n_offmode <- NA_integer_
    lpc <- tryCatch(colMeans(posterior::extract_variable_matrix(model_sample$draws("lp__"), "lp__")),
      error = function(e) NULL)
    if (!is.null(lpc)) { g <- max(lpc) - lpc; lp_gap_max <- max(g); n_offmode <- sum(g > 5) }
    sampler_diag <- list(
      n_div = sum(ds$num_divergent),
      n_tree = sum(ds$num_max_treedepth),
      ebfmi_min = suppressWarnings(min(ds$ebfmi)),
      rhat_max = mixing$rhat_max,
      rhat_loadings = mixing$rhat_loadings,
      rhat_estimands = mixing$rhat_estimands,
      rhat_cor_sq = mixing$rhat_cor_sq,
      rhat_M = mixing$rhat_M,
      ess_delta1 = mixing$ess_delta1,
      n_offmode = n_offmode,
      lp_gap_max = lp_gap_max
    )

    # Log only problematic fits, so a clean run leaves the log to the progress lines.
    if (!is.null(log_file) && (
      sampler_diag$n_div > 0 || sampler_diag$n_tree > 0 ||
        (is.finite(sampler_diag$ebfmi_min) && sampler_diag$ebfmi_min < 0.3) ||
        (is.finite(sampler_diag$ess_delta1) && sampler_diag$ess_delta1 < 100) ||
        (is.finite(sampler_diag$rhat_max) && sampler_diag$rhat_max > 1.01) ||
        (is.finite(sampler_diag$n_offmode) && sampler_diag$n_offmode > 0))) {
      cat(sprintf(
        "[%s] %s  STAN div=%d treedepth=%d ebfmi_min=%.2f ess_delta1=%.0f rhat_max=%.3f rhat_M=%.3f rhat_est=%.3f rhat_load=%.3f rhat_corsq=%.3f lp_gap_max=%.1f\n",
        format(Sys.time(), "%H:%M"), log_label,
        sampler_diag$n_div, sampler_diag$n_tree, sampler_diag$ebfmi_min,
        sampler_diag$ess_delta1, sampler_diag$rhat_max, sampler_diag$rhat_M,
        sampler_diag$rhat_estimands, sampler_diag$rhat_loadings, sampler_diag$rhat_cor_sq,
        sampler_diag$lp_gap_max
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
      loc_cor_pval = loc_cor_pval,
      sampler_diag = sampler_diag,
      draws = if (!is.null(return_draws)) model_sample$draws(return_draws) else NULL,
      sampler_draws = if (!is.null(return_draws)) model_sample$sampler_diagnostics() else NULL
    ))
  }
}
