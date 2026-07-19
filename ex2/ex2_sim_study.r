## This file runs the simulation study for the intercepts example, fitting the
## fitting the models with and without intercepts to data generated from an
## adversarial process whereby some units are simulated to have different
## long run means than the treated but spuriously correlation in the pre-
## treatment period only.

library(foreach)
library(doParallel)

cl <- makeCluster(detectCores() - 1, outfile = "")
registerDoParallel(cl)

source("../sample_model.r")
source("../plotting.r")

ruv <- function(d) {
  v <- rnorm(d)
  uv <- v / sqrt(sum(v * v))
  return(uv)
}

sim_model_intercepts <- function(
  N_unc = 5, N_comp_true = 2, N_comp_spur = 2, T_times = 20, T_treated = 5, K_unc = 3, sim = 0.9) {

    treat_time <- T_times - T_treated + 1
    
    f_treat <- arima.sim(model = list(ar = 0.96), n = T_times)

    f_alt <- f_treat + 
      c(rep(0, T_times - T_treated), (-1 / 2) * seq(T_treated))
    
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
    K_latent <- 2 + K_unc
    loads <- matrix(nrow = N_units, ncol = K_latent)

    loads[1, ] <- c(1, rep(0, K_latent - 1))
    
    for(n in 1:N_comp_true) {
      loads[1 + n, ] <- c(sqrt(0.9), 0, sqrt(1 - 0.9) * ruv(K_latent - 2))
    }

    for(n in 1:N_comp_spur) {
      loads[1 + N_comp_true + n, ] <- c(0, sqrt(sim), sqrt(1 - sim) * ruv(K_latent - 2))
    }

    for(n in 1:N_unc) {
      loads[1 + N_comp_true + N_comp_spur + n, ] <- c(0, 0 , ruv(K_latent - 2))
    }

    lat <- loads %*% facs
    intercepts <- c(
      1, rnorm(N_comp_true, mean = 1, sd = 0), 
      rnorm(N_comp_spur, mean = -5, sd = 0),
      rnorm(N_unc, mean = 0, sd = 1)
    )
    
    Y <- lat + intercepts + rnorm(N_units * T_times, sd = 0.02)

    return(t(Y))

}

run_sim_intercepts <- function(N_comp, sim, K_latent, post_check = FALSE) {
  test_ys <- sim_model_intercepts(N_unc = 7 - N_comp, N_comp_spur = N_comp, K_unc = 3, sim = sim)

  N_units <- ncol(test_ys)
  T_times <- nrow(test_ys)
  overall_scales <- apply(test_ys, 2, sd)

  fits <- list()

  fits$nint <- sample_model(N_units = 10, T_times = 20, K_latent = 6,
                            overall_scales = overall_scales, err_scale = 0.2,
                            data = test_ys,
                            autocor_a = 99, autocor_b = 1,
                            nonstationary = FALSE, num_treated = 5,
                            type = "posterior", quiet = TRUE, ad = 0.8, iter = 500,
                            n_chains = 1)
  fits$nint$name <- "nint"

  fits$ints <- sample_model(N_units = 10, T_times = 20, K_latent = 6,
                            overall_scales = overall_scales, err_scale = 0.2,
                            data = test_ys,
                            autocor_a = 99, autocor_b = 1,
                            nonstationary = FALSE, num_treated = 5,
                            include_ints = TRUE, int_scale = 1, #max(overall_scales)
                            type = "posterior", quiet = TRUE, iter = 500,
                            n_chains = 1)
  fits$ints$name <- "ints"

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

  }) |> purrr::list_flatten()

  if (post_check) {
    # TBD
    # check_plot <- plot_data_matrix_post(test_ys, fits$stat2$y_pred)
    # ggsave(
    #     check_plot, file=paste0("sim_stat_figs/check_plot_", i, ".pdf"),
    #     device = "pdf", height = 4, width = 8
    # )
  }

  return(res)
}

run_sim_study_intercepts <- function(K_latent = 3, reps, N_comps, sims, post_check = FALSE) {
  # Functions called inside the worker must be exported explicitly; foreach does
  # not auto-export them across the %:% nesting (K_latent, a local variable, is
  # auto-exported). posterior is attached because sample_model() calls
  # extract_variable_array() unqualified.
  exp_vars <- c('run_sim_intercepts', 'sim_model_intercepts', 'ruv',
                'sample_model', 'ife_mod')
  exp_packages <- c('cmdstanr', 'posterior')

  # The %:% operator flattens the three nested loops into a single stream of
  # sim x N_comp x rep tasks, distributed across workers by %dopar%.
  study_res <-
    foreach(sim = sims, .combine = 'rbind') %:%
    foreach(N_comp = N_comps, .combine = 'rbind') %:%
    foreach(
      r = seq_len(reps), .combine = 'rbind',
      .export = exp_vars, .packages = exp_packages
    ) %dopar% {
      print(paste0(">>>> Condition (sim = ", sim, ", num_comp = ", N_comp, "), rep ", r, "."))
      unit_res <- run_sim_intercepts(N_comp = N_comp, sim = sim, K_latent = K_latent)
      unit_res <- c(unit_res, list(sim = sim, num_comp = N_comp))
      as.data.frame(unit_res)
    }

  return(study_res)
}

sim_study_ints <- run_sim_study_intercepts(
  reps = 1000, # 50
  N_comps = c(1, 3),
  sims = c(0.7, 0.9)
)

stopCluster(cl)

absz_res <- sim_study_ints |> dplyr::select(contains("absz"))
print(colMeans(absz_res))

save(sim_study_ints, file="sim_study_ints.RData")
