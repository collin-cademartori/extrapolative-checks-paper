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
    autocor_a, autocor_b, alpha_diag = 0, nonstationary, int_scale = 1, int_loc = 0, include_ints = FALSE, include_factor_means = FALSE,
    num_treated, delta_scale = 0, type = "prior_pred",
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
  stopifnot(num_treated == 0 || delta_scale > 0)

  # If CMDSTAN_OUTPUT_DIR is set, CmdStan writes its CSV output there.
  out_dir <- NULL
  csv_root <- Sys.getenv("CMDSTAN_OUTPUT_DIR", "")
  if (nzchar(csv_root)) {
    dir.create(csv_root, showWarnings = FALSE, recursive = TRUE)
    # tempfile() gives a unique name without advancing the global RNG.
    out_dir <- file.path(csv_root, basename(tempfile("fit_")))
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  }

  sample_index <- if (is.null(seed)) {
    sample.int(iter, size = iter)
  } else {
    have_state <- exists(".Random.seed", envir = globalenv())
    old_state <- if (have_state) get(".Random.seed", envir = globalenv()) else NULL
    set.seed(seed)
    idx <- sample.int(iter, size = iter)
    if (have_state) assign(".Random.seed", old_state, envir = globalenv())
    idx
  }

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
    nonstationary = nonstationary,
    unit_intercepts = include_ints,
    factor_means = include_factor_means,
    sample_posterior = (type == "posterior"),
    num_treated = num_treated,
    gamma_scale = int_scale,
    gamma_loc = int_loc,
    alpha_diag = alpha_diag,
    # Prior scale for the treatment effect, on the data's scale.
    delta_scale = delta_scale
  )

  # Initial value for the noise scale when it has a prior, used when there is no Pathfinder
  # initialization.
  scale_init <- if (err_scale == 0) list(tau_param = array(err_scale_mean, dim = 1)) else NULL

  # Optional Pathfinder initialization
  init_arg <- if (isTRUE(pathfinder_init)) {
    pfi <- pathfinder_inits(ife_mod, stat_data, n_chains, seed = seed, quiet = quiet,
      output_dir = out_dir)
    if (is.null(pfi)) (if (is.null(init)) (if (is.null(scale_init)) 2 else rep(list(scale_init), n_chains)) else init) else pfi
  } else if (is.null(init)) (if (is.null(scale_init)) 2 else rep(list(scale_init), n_chains)) else init

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
    seed = seed,
    output_dir = out_dir
  )

  # Sampler diagnostics for posterior fits.
  sampler_diag <- NULL
  if (type == "posterior") {
    ds <- model_sample$diagnostic_summary(quiet = TRUE)
    # Maximum R-hat over groups of quantities:
    #   rhat_max       : all parameters
    #   rhat_loadings  : Lambda and Phi_innovations
    #   rhat_estimands : delta, the noise scale, rho, gamma, and omega_sq
    #   rhat_cor_sq    : the model's squared correlations between the treated and untreated units
    #   rhat_M         : the latent mean matrix Lambda_Phi. The likelihood depends on Lambda and Phi
    #                    only through this product.
    mixing <- tryCatch(
      {
        sv <- model_sample$metadata()$stan_variables
        pick <- function(nms) intersect(nms, sv)
        all_vars <- pick(c("Lambda", "Phi_innovations", "rho", "gamma_raw",
          "delta_raw", "omega_sq_param", "Phi_means_param", "tau_param"))
        ps <- model_sample$summary(all_vars, "rhat", "ess_bulk")
        grp <- function(bases) {
          v <- pick(bases); if (!length(v)) return(NA_real_)
          suppressWarnings(max(model_sample$summary(v, "rhat")$rhat, na.rm = TRUE))
        }
        list(
          rhat_max = suppressWarnings(max(ps$rhat, na.rm = TRUE)),
          rhat_loadings = grp(c("Lambda", "Phi_innovations")),
          rhat_estimands = grp(c("delta_raw", "tau_param", "rho",
            "gamma_raw", "omega_sq_param")),
          rhat_cor_sq = grp("cor_sq"),
          rhat_M = grp("Lambda_Phi"),
          ess_delta1 = ps$ess_bulk[match("delta_raw[1]", ps$variable)]
        )
      },
      error = function(e) list(rhat_max = NA_real_, rhat_loadings = NA_real_,
        rhat_estimands = NA_real_, rhat_cor_sq = NA_real_, rhat_M = NA_real_, ess_delta1 = NA_real_)
    )

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

    if (!is.null(log_file) && (
      sampler_diag$n_div > 0 || sampler_diag$n_tree > 0 ||
        (is.finite(sampler_diag$ebfmi_min) && sampler_diag$ebfmi_min < 0.3) ||
        (is.finite(sampler_diag$ess_delta1) && sampler_diag$ess_delta1 < 100) ||
        (is.finite(sampler_diag$rhat_max) && sampler_diag$rhat_max > 1.01) ||
        (is.finite(sampler_diag$rhat_M) && sampler_diag$rhat_M > 1.01) ||
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

    out <- list(
      ys = ys_prior,
      ys_latent = ys_latent
    )
    try(unlink(model_sample$output_files(), force = TRUE), silent = TRUE)
    return(out)
  } else if (type == "posterior") {
    y_means_all <-
      extract_variable_array(model_sample$draws("Y_latent"), "Y_latent")
    y_means_post <- y_means_all[sample_index, 1, , ]

    y_pred_all <-
      extract_variable_array(model_sample$draws("Y_pred"), "Y_pred")
    y_pred_post <- y_pred_all[sample_index, 1, , ]

    effects <-
      extract_variable_array(model_sample$draws("delta"), "delta")[, 1, ]
    effect_means <- colMeans(effects)
    effect_sds <- apply(effects, 2, sd)

    err_scale_mat <- posterior::as_draws_matrix(model_sample$draws("tau"))
    err_scale <- as.numeric(err_scale_mat[, 1])
    mad <- mean(as.numeric(model_sample$draws("mean_abs_diffs")))

    cor_sq <- extract_variable_array(model_sample$draws("cor_sq"), "cor_sq")[, 1, ]
    cor_sq_mean <- colMeans(cor_sq)

    abs_cor_pred <- as.numeric(model_sample$draws("time_cor_pred"))
    # Observed S1: the mean over untreated units of |correlation with time|, as in the Stan model.
    abs_cor_data <- mean(abs(cor(data[, -1, drop = FALSE], seq(nrow(data)))))
    time_cor_pval <- mean(abs_cor_pred > abs_cor_data)

    # Observed S2 on the pre-treatment window: the correlation across untreated units between each
    # unit's correlation with the treated unit and its mean.
    loc_cor_pred <- as.numeric(model_sample$draws("loc_cor_pred"))
    pre_times <- seq_len(nrow(data) - num_treated)
    untreated_pre <- data[pre_times, -1, drop = FALSE]
    cor_with_treated <- as.numeric(cor(untreated_pre, data[pre_times, 1]))
    unit_location <- colMeans(untreated_pre)
    loc_cor_data <- cor(cor_with_treated, unit_location)
    loc_cor_pval <- mean(loc_cor_pred > loc_cor_data)

    # Squared sample correlations with the treated unit, compared with the model's cor_sq.
    y_cor_sq <- cor(data)[1, 2:ncol(data)]^2
    cor_err_mean <- rowMeans(abs(y_cor_sq - t(cor_sq)))

    out <- list(
      y_means = y_means_post,
      y_pred = y_pred_post,
      effect_means = effect_means,
      effect_sds = effect_sds,
      mean_abs_diffs = mad,
      cor_sq = cor_sq_mean,
      abs_cors_err = cor_err_mean,
      err_scale = err_scale,
      err_scale_all = err_scale_mat,
      time_cor_pval = time_cor_pval,
      loc_cor_pval = loc_cor_pval,
      sampler_diag = sampler_diag,
      draws = if (!is.null(return_draws)) model_sample$draws(return_draws) else NULL,
      sampler_draws = if (!is.null(return_draws)) model_sample$sampler_diagnostics() else NULL
    )
    try(unlink(model_sample$output_files(), force = TRUE), silent = TRUE)
    return(out)
  }
}


# ---------------------------------------------------------------------------------------------
# Escalation: A fit that fails the convergence criteria is refit with more iterations 
# (for slow mixing) or a higher adapt_delta (for divergences), until it passes 
# we hit the limit of re-runs (max 2 re-runs).
escalation_ladder <- function(iter, warm, ad_floor = 0.95, rhat_M = 1.01,
                              rhat_est = 1.01, ess = 400, div_rate = 0.001) {
  stopifnot(length(iter) == length(warm))
  list(iter = as.integer(iter), warm = as.integer(warm), ad_floor = ad_floor,
       rhat_M = rhat_M, rhat_est = rhat_est, ess = ess, div_rate = div_rate,
       max_rounds = length(iter) + 1L)   # including the initial fit
}

fit_with_escalation <- function(args, seeds, label, progress_log, ladder) {
  iter <- args$iter
  warm <- if (is.null(args$iter_warm)) args$iter else args$iter_warm
  ad <- args$ad
  it_level <- 0L
  fit <- NULL
  for (round in seq_len(ladder$max_rounds)) {
    a <- args
    a$iter <- iter
    a$iter_warm <- warm
    a$ad <- ad
    a$seed <- seeds[round]
    a$log_file <- progress_log
    a$log_label <- if (round == 1L) label else
      sprintf("%s [round %d: iter=%d warm=%d ad=%.3f]", label, round, iter, warm, ad)
    fit <- do.call(sample_model, a)
    sd_ <- fit$sampler_diag
    n_draws <- iter * a$n_chains

    slow <- (is.finite(sd_$rhat_M) && sd_$rhat_M > ladder$rhat_M) ||
      (is.finite(sd_$rhat_estimands) && sd_$rhat_estimands > ladder$rhat_est) ||
      (is.finite(sd_$ess_delta1) && sd_$ess_delta1 < ladder$ess)
    divergent <- is.finite(sd_$n_div) && sd_$n_div > ladder$div_rate * n_draws
    if ((!slow && !divergent) || round == ladder$max_rounds) break

    if (slow && it_level < length(ladder$iter)) {
      it_level <- it_level + 1L
      iter <- ladder$iter[it_level]
      warm <- ladder$warm[it_level]
    }
    ad <- max(if (divergent) min(0.99, 1 - (1 - ad) / 4) else ad, ladder$ad_floor)
  }
  fit$n_rounds <- round
  fit$final_iter <- iter
  fit$final_warm <- warm
  fit$final_ad <- ad
  fit
}
