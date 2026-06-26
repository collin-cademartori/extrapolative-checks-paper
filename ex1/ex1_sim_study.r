## This file runs the simulation study for the nonstationary example, fitting the
## nonstationary model and two stationary models (with weak and strong priors on the 
## iid error scale) to samples from the prior predictive distribution of the
## nonstationary model.

source(".,/sample_model.r")

run_sim_stat <- function(test_data, i, K_latent, post_check = FALSE) {
  test_ys <- test_data$ys[i, , ]
  N_units <- ncol(test_ys)
  T_times <- nrow(test_ys)
  overall_scales_stat <- apply(test_ys, 2, sd)
  overall_scales_nonstat <- apply(test_ys, 2, function(y) sd(diff(y)))

  fits <- list()

  fits$nonstat <- sample_model(
    overall_scales = overall_scales_nonstat, err_scale = 0,
    err_scale_mean = 3,
    err_scale_sd = 2,
    data = test_ys,
    autocor_a = 8, autocor_b = 2,
    nonstationary = TRUE, num_treated = 5,
    type = "posterior", K_latent = K_latent,
    iter = 500
  )
  fits$nonstat$name <- "nonstat"

  fits$stat2 <- sample_model(
    overall_scales = overall_scales_stat, err_scale = 0,
    err_scale_mean = 0.2, # 0.1
    err_scale_sd = 0.2, # 0.5
    data = test_ys,
    autocor_a = 99, autocor_b = 1, # autocor_a = 99
    nonstationary = FALSE, num_treated = 5,
    type = "posterior", K_latent = K_latent,
    iter = 500
  )
  fits$stat2$name <- "stat2"

  fits$stat1 <- sample_model(
    overall_scales = overall_scales_stat, err_scale = 0.05,
    # err_scale_mean = 0.1, # 0.1
    # err_scale_sd = 0.01, # 0.5
    data = test_ys,
    autocor_a = 99, autocor_b = 1, # autocor_a = 99
    nonstationary = FALSE, num_treated = 5,
    type = "posterior", K_latent = K_latent,
    iter = 500
  )
  fits$stat1$name <- "stat1"

  res <- fits |> purrr::map(function(pfit) {

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

    absz <- abs(pfit$effect_means / pfit$effect_sds)
    res[paste0("absz_", seq_along(absz))] <- absz

    pmean <- pfit$effect_means
    res[paste0("mean_", seq_along(pmean))] <- pmean

    psds <- pfit$effect_sds
    res[paste0("sd_", seq_along(psds))] <- psds

    pred_mad <- pfit$mean_abs_diffs
    res$pred_mad <- pred_mad

    return(res)

  }) |> purrr::list_flatten()

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

run_sim_study_stat <- function(K_latent = 3, reps, post_check = FALSE) {
  test_data <- sample_model(overall_scales = rep(1, 8), err_scale = 3, # err_scale = 5
                            autocor_a = 8, autocor_b = 2,
                            nonstationary = TRUE, num_treated = 0,
                            type = "prior_pred", K_latent = K_latent)

  study_res <- data.frame()
  study_units <- sample.int(1000, size = reps, replace = FALSE)
  iter <- 1
  for(s in study_units) {
    unit_res <- as.data.frame(run_sim_stat(test_data, s, K_latent, post_check))
    study_res <- rbind(study_res, unit_res)

    writeLines(paste0("Iteration ", iter, ":"))
    absz_res <- study_res |> dplyr::select(contains("absz"))
    print(colMeans(absz_res))
    writeLines("-----------------------------")

    iter <- iter + 1
  }

  return(study_res)
}

study_reps <- 300
sim_study_stat <- run_sim_study_stat(K_latent = 4, study_reps)

save(sim_study_stat, file="sim_study_ns.RData")

study_reps <- 1
sim_study_stat <- run_sim_study_stat(K_latent = 4, study_reps, post_check = TRUE)
