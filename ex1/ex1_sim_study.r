## This file runs the simulation study for the nonstationary example, fitting the
## nonstationary model and two stationary models (with weak and strong priors on the 
## iid error scale) to samples from the prior predictive distribution of the
## nonstationary model.

library(foreach)
library(doParallel)

cl <- makeCluster(detectCores() - 1, outfile = "")
registerDoParallel(cl)

source("../sample_model.r")
source("../plotting.r")

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
    iter = 400,
    n_chains = 1
  )
  fits$nonstat$name <- "nonstat"

  fits$stat2 <- sample_model(
    overall_scales = overall_scales_stat, err_scale = 0,
    err_scale_mean = 0.2, # 0.1
    err_scale_sd = 0.2, # 0.5
    data = test_ys,
    autocor_a = 97, autocor_b = 3,
    nonstationary = FALSE, num_treated = 5,
    type = "posterior", K_latent = K_latent,
    iter = 400,
    n_chains = 1
  )
  fits$stat2$name <- "stat2"

  fits$stat1 <- sample_model(
    overall_scales = overall_scales_stat, err_scale = 0.05,
    # err_scale_mean = 0.1, # 0.1
    # err_scale_sd = 0.01, # 0.5
    data = test_ys,
    autocor_a = 97, autocor_b = 3,
    nonstationary = FALSE, num_treated = 5,
    type = "posterior", K_latent = K_latent,
    iter = 400,
    n_chains = 1
  )
  fits$stat1$name <- "stat1"

  res <- fits |> purrr::map(function(pfit) {

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
  test_data <- sample_model(overall_scales = rep(1, 8), err_scale = 3,
                            autocor_a = 8, autocor_b = 2,
                            nonstationary = TRUE, num_treated = 0,
                            type = "prior_pred", K_latent = K_latent,
                            iter = 2 * reps)

  #study_res <- data.frame()
  study_units <- sample.int(2 * reps, size = reps, replace = FALSE)
  iter <- 1
  exp_vars <- c('run_sim_stat', 'sample_model', 'ife_mod', 'plot_post_fits_all', 'plot_data_matrix_post')
  exp_packages <- c('cmdstanr', 'posterior', 'forcats', 'dplyr', 'ggplot2')
  study_res <- 
    foreach(
      s = study_units, iter = seq(reps),
      .combine='rbind', .export = exp_vars, .packages = exp_packages
    ) %dopar% {
      print(paste0(">>>> Beginning iteration ", iter, "."))
      unit_res <- as.data.frame(run_sim_stat(test_data, s, K_latent, post_check))
    }

  return(study_res)
}

study_reps <- 100 #2000
sim_study_stat <- run_sim_study_stat(K_latent = 4, study_reps)

stopCluster(cl)

absz_res <- sim_study_stat |> dplyr::select(contains("absz"))
print(colMeans(absz_res))

ci_width_res <- sim_study_stat |> dplyr::select(contains("width"))
print(colMeans(ci_width_res))

ci_cov_res <- sim_study_stat |> dplyr::select(contains("perc"))
print(colMeans(ci_cov_res))

time_pval <- sim_study_stat |> dplyr::select(contains("time_cor"))
print(colMeans(time_pval))

save(sim_study_stat, file="sim_study_ns.RData")

# study_reps <- 1
# sim_study_stat <- run_sim_study_stat(K_latent = 4, study_reps, post_check = TRUE)
