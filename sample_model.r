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
    autocor_a, autocor_b, alpha_diag = 0, nonstationary, absolute_error = FALSE, int_scale = 1, int_loc = 0, include_ints = FALSE, include_factor_means = FALSE,
    num_treated, delta_scale = 0, type = "prior_pred",
    iter = 1000, iter_warm = NULL, quiet = TRUE, 
    ad = 0.98, max_treedepth = 10, n_chains = 4, parallel_chains = 1,
    seed = NULL, log_file = NULL, log_label = NULL,
    return_draws = NULL, init = NULL,
    pathfinder_init = FALSE) {
  stopifnot(type %in% c("prior_pred", "posterior"))
  stopifnot(0 < autocor_a)
  stopifnot(0 < autocor_b)
  # STRICTLY less than: the saturated case K = M is excluded deliberately. At K = M the loadings
  # are a full invertible lower-triangular matrix, so Phi * Lambda' spans EVERY T x M matrix, the
  # signal is unrestricted, and the likelihood is UNBOUNDED as the error scale goes to zero along
  # the interpolating ridge -- the same pathology as a normal mixture whose component variance
  # collapses onto a data point. Measured, not assumed: at K = M = 8 Stan's optimizer finds that
  # mode in under half a second (eta ~ 2e-4, residual 0 to four decimals) from every random start.
  # The POSTERIOR is still proper and MCMC explores its bulk correctly, so this is not a hard
  # error -- but pathfinder_init = TRUE, which every study fit uses, is optimization-based and can
  # seed chains into that trap. Neither study needs K = M (ex1 fits 4, ex2 fits 3), so the guard
  # costs nothing and removes a real hazard. ife_named.stan's own declaration allows K <= M.
  stopifnot(K_latent < N_units)
  stopifnot(num_treated >= 0 && num_treated < T_times)
  stopifnot(err_scale > 0 || (err_scale_mean > 0 && err_scale_sd > 0))

  # Where CmdStan writes its per-chain CSVs. Left to cmdstanr these go to R's session tempdir, which
  # on the study machine is /tmp -- a TMPFS, so those CSVs are held in RAM against a cap that is
  # typically half of physical memory. That is what killed two long runs: once as an OOM kill, and
  # once (after the $draws() fix cut R's own heap use) as tmpfs hitting its size cap and returning
  # ENOSPC, which surfaced as "disk is full in the temporary directory" from data.table and then
  # "No chains finished successfully" as CmdStan could no longer write at all.
  #
  # The files are unavoidably live for the DURATION of a fit -- roughly 264 MB at the round-2 rung
  # (8000 iterations x 3 chains) and 396 MB at round 3 -- so unlinking them afterwards, which
  # sample_model already does, is not enough on its own when ~19 workers are fitting at once.
  #
  # Set CMDSTAN_OUTPUT_DIR to a path on real disk to move them off tmpfs entirely. Each fit gets its
  # own subdirectory, removed on exit so that a fit which ERRORS does not leak its CSVs either --
  # the ordinary cleanup below only runs on the success path.
  out_dir <- NULL
  csv_root <- Sys.getenv("CMDSTAN_OUTPUT_DIR", "")
  if (nzchar(csv_root)) {
    dir.create(csv_root, showWarnings = FALSE, recursive = TRUE)
    # basename(tempfile()) for a unique name: it uses R's internal counter and, unlike sample(),
    # does NOT advance the global RNG, which would break the study's reproducibility.
    out_dir <- file.path(csv_root, basename(tempfile("fit_")))
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)
  }

  # The returned draws are shuffled, so that consecutive datasets do not come from the same chain
  # (CmdStan writes draws chain-major). Two constraints on where that index comes from:
  #   * It must be drawn *before* $sample(): cmdstanr's $sample() advances R's global RNG, so
  #     drawing it afterward would make the ordering irreproducible. This must precede the fit.
  #   * It is drawn from this call's own `seed`, not from R's global stream. Under the global
  #     stream the ordering depended on how much RNG every earlier line happened to consume, so
  #     inserting any sample()/rnorm() anywhere upstream silently renumbered the datasets -- and a
  #     study that indexes datasets by position (ex1_sim_study.r's study_units) would then be
  #     fitting different data with no visible change. The global stream is saved and restored
  #     around the draw so callers that rely on it are unaffected.
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
    alpha_diag = alpha_diag,
    # 0 = err_sd[n] is tau[n]*sigma[n] (a ratio to each unit's fixed scale); 1 = err_sd[n] is a
    # single absolute eta shared by all units. In absolute mode err_scale / err_scale_mean /
    # err_scale_sd are read on the DATA's scale, not as ratios -- the caller must pass values on the
    # right scale, and nothing can detect the mistake if they do not.
    absolute_error = as.integer(absolute_error),
    # Prior scale for the treatment effect, on the data's own scale. 0 = inherit sigma[1], the
    # historical behaviour. Pass a common positive value across arms being compared so they share a
    # prior on the estimand; see the note in ife_named.stan's data block.
    delta_scale = delta_scale
  )

  # Fallback init for ABSOLUTE mode. eta lives on the data's own scale (order 2 here), but Stan's
  # default init draws a positive parameter from exp(U(-2, 2)) = (0.14, 7.4); the low end makes the
  # likelihood extremely sharp at a random Lambda/Phi and every proposal is rejected. In RATIO mode
  # the same draw is multiplied by sigma[n] (order 8), which rescues it. Pathfinder finds a sensible
  # eta on its own (~1.2 against a truth of 2 in testing), so this only matters when Pathfinder is
  # off or has failed. Only absolute mode is touched; the ratio configuration is unchanged.
  scale_init <- if (isTRUE(absolute_error) && err_scale == 0)
    list(tau_param = array(err_scale_mean, dim = 1)) else NULL

  # Optional Pathfinder warm start: seed every chain from a draw in the dominant lp mode, discarding modes which are vanishingly tiny compared to the dominant mode.
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
        # rhat_M is NOT a subset of rhat_max (that is computed over the parameters, while
        # Lambda_Phi is a transformed quantity), so M can in principle be poor while the raw
        # parameters look fine. Trip on it in its own right.
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
    # Drop the CmdStan CSVs now that everything is materialized in R. See the note at the
    # matching call in the posterior branch below for why this cannot be left to the GC.
    try(unlink(model_sample$output_files(), force = TRUE), silent = TRUE)
    return(out)
  } else if (type == "posterior") {
    y_means_all <-
      extract_variable_array(model_sample$draws("Y_latent"), "Y_latent")
    y_means_post <- y_means_all[sample_index, 1, , ]

    y_pred_all <-
      extract_variable_array(model_sample$draws("Y_pred"), "Y_pred")
    y_pred_post <- y_pred_all[sample_index, 1, , ]

    # $draws("delta"), NOT $draws(): the bare call materializes every one of the model's ~1270
    # columns and cmdstanr then caches that array for the life of the fit, so it stays resident.
    # Measured at 10.1 KB per stored draw, which is 61 MB at iter = 2000 x 3 chains but 607 MB at
    # iter = 20000 -- and the escalation ladder reaches iter = 20000. Naming the variable brings it
    # down to ~4 MB there. This is what put the study into the OOM killer at round 3.
    effects <-
      extract_variable_array(model_sample$draws("delta"), "delta")[, 1, ]
    effect_means <- colMeans(effects)
    effect_sds <- apply(effects, 2, sd)

    # tau is now a VECTOR over units, so draws("tau") is [draws x M]. err_scale keeps its scalar
    # per-draw meaning by reporting the TREATED unit's tau -- unit 1, the only one delta depends on,
    # and the one every tau summary downstream is about. The full matrix is returned alongside as
    # err_scale_all so nothing is lost. (Which unit to headline is a reporting choice, not a
    # modelling one; change it here if the treated unit is not the right summary.)
    err_scale_mat <- posterior::as_draws_matrix(model_sample$draws("tau"))
    err_scale <- as.numeric(err_scale_mat[, 1])
    mad <- mean(as.numeric(model_sample$draws("mean_abs_diffs")))

    cor_sq <- extract_variable_array(model_sample$draws("cor_sq"), "cor_sq")[, 1, ]
    cor_sq_mean <- colMeans(cor_sq)

    abs_cor_pred <- as.numeric(model_sample$draws("time_cor_pred"))
    # Observed S1, matching the Stan definition exactly: the MEAN over untreated units of each
    # unit's |correlation with time|. Both sides must average over the same set (columns 2..M of
    # whatever matrix was fitted), or the p-value compares different quantities.
    abs_cor_data <- mean(abs(cor(data[, -1, drop = FALSE], seq(nrow(data)))))
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

    # SQUARED sample correlation, to match cor_sq, which ife_named.stan defines as a squared
    # correlation (squared factor covariance over the product of total variances). This used to be
    # abs(cor(...)) -- an unsquared |r| differenced against an r^2 -- so the "error" was dominated by
    # the gap between the two quantities rather than by any model misfit. Measured on ex2: for units
    # whose r^2 the model put at 0.887, the reported error was 0.052, essentially exactly the
    # |r| - r^2 = 0.942 - 0.887 that the mismatch alone produces.
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
    # Delete the CmdStan output CSVs explicitly. cmdstanr removes them only when the CmdStanFit
    # object is garbage collected, and R's GC is driven by *heap* pressure -- which for one fit is a
    # few MB, while the CSVs it leaves behind are 66 MB at iter = 2000 x 3 chains and 264 MB at
    # iter = 8000 (this model writes 1272 columns, ~11 KB per stored draw). Nothing in R's heap
    # therefore prompts the finalizers to run, so across a few hundred fits the temp filesystem
    # accumulates tens of GB and the study dies on a disk quota. By this point every quantity we
    # return has already been read into R (cmdstanr caches draws after the first read), so the files
    # are safe to drop. Must come AFTER `out` is built: the $draws()/$sampler_diagnostics() calls
    # above are the last readers.
    try(unlink(model_sample$output_files(), force = TRUE), silent = TRUE)
    return(out)
  }
}


# ---------------------------------------------------------------------------------------------
# Convergence escalation, shared by both simulation studies.
#
# A minority of fits mix slowly in the loadings, and a few in the latent means M as well; longer
# chains fix most of them, and a higher adapt_delta fixes the divergent ones. Rather than
# post-processing or discarding those fits, refit them in place with progressively more computation
# until they meet the criterion or the ladder is spent.
#
# This is adaptive computation, not selection: the stopping rule is a function of convergence
# diagnostics only. Under a correct sampler those are independent of the estimand, so escalating does
# not bias delta. (Discarding non-converged fits instead WOULD be a selection procedure, because the
# fits that struggle are the harder datasets, which are plausibly not exchangeable with the rest.)
# Note the flip side, worth stating wherever the diagnostic distribution is reported: because a fit
# stops as soon as it passes, the RECORDED R-hat distribution is stopped-on-success and so reads low.
# Report the ladder (n_rounds) alongside it rather than the final R-hats alone.
#
# The criterion is rhat_M, not rhat_max: the likelihood depends on (Lambda, Phi) only through
# M = Lambda_Phi, so slow mixing confined to the loadings cannot affect delta, while anything that
# could must show up in M. This holds for BOTH examples -- with unit intercepts and factor means the
# intercept enters Y_means additively (gamma[n] + sigma[n] * Lambda_Phi[:,n]) and the factor means
# fold into Phi, so M remains the sufficient statistic; the conditioning set merely grows to include
# gamma and omega_sq, which rhat_estimands already covers. ess_delta1 guards the estimand directly.
#
# The ladder itself is per-example, since the two studies have different base configurations. Build
# one with escalation_ladder() and pass it in.
# A ZERO-LENGTH ladder disables escalation: max_rounds becomes 1, so fit_with_escalation() breaks
# after the first fit by construction. That is the honest way to express "no escalation" -- the
# alternative, setting rhat_M = Inf / ess = 0 / div_rate = 1 so no criterion can ever fire, works
# but hides the intent behind three unreachable numbers.
escalation_ladder <- function(iter, warm, ad_floor = 0.95, rhat_M = 1.01,
                              ess = 400, div_rate = 0.001) {
  stopifnot(length(iter) == length(warm))
  list(iter = as.integer(iter), warm = as.integer(warm), ad_floor = ad_floor,
       rhat_M = rhat_M, ess = ess, div_rate = div_rate,
       max_rounds = length(iter) + 1L)   # rounds INCLUDING the unescalated first attempt
}

fit_with_escalation <- function(args, seeds, label, progress_log, ladder) {
  iter <- args$iter
  warm <- if (is.null(args$iter_warm)) args$iter else args$iter_warm
  ad <- args$ad
  # Rungs of ladder$iter already spent. Tracked separately from `round` so that a round bought
  # purely by divergences does not consume a rung of the iteration ladder.
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
      (is.finite(sd_$ess_delta1) && sd_$ess_delta1 < ladder$ess)
    divergent <- is.finite(sd_$n_div) && sd_$n_div > ladder$div_rate * n_draws
    if ((!slow && !divergent) || round == ladder$max_rounds) break
    # Slow mixing buys iterations and matching warmup; divergences buy adapt_delta on top of the
    # floor that every escalated round gets regardless.
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
