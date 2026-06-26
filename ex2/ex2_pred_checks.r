#################################
## Per-Unit Intercepts Example ##
#################################

source("../sample_model.r")

# Step 1. Generate a sample (with unit-level intercepts) with three types of untreated units:
#  1. Units which are uncorrelated with the treated unit
#  2. Units which are correlated with the treated unit but with a lower intercept
#  3. Units which are correlated with the treated unit and have a comparable intercept
#  In this sample, we will also make it so that the second and third groups of units have
#  distinctly different behavior in unit 1's post-treatment period
 
ruv <- function(d) {
  v <- rnorm(d)
  uv <- v / sqrt(sum(v * v))
  return(uv)
}

sim_model_intercepts <- function(
  N_unc = 5, N_comp = 2, T_times = 30, T_treated = 5, K_unc = 3, sim = 0.9) {

    treat_time <- T_times - T_treated + 1
    
    f_treat <- arima.sim(model = list(ar = 0.96), n = T_times)

    # f_alt <- sqrt(0.95) * f_treat + sqrt(0.05) * rnorm(T_times) + 
    #   c(rep(0, T_times - T_treated), -1 * seq(T_treated))

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

    N_units <- 1 + 2 * N_comp + N_unc
    K_latent <- 2 + K_unc
    loads <- matrix(nrow = N_units, ncol = K_latent)

    loads[1, ] <- c(1, rep(0, K_latent - 1))
    
    for(n in 1:N_comp) {
      loads[1 + n, ] <- c(sqrt(sim), 0, sqrt(1 - sim) * ruv(K_latent - 2))
    }

    for(n in 1:N_comp) {
      loads[1 + N_comp + n, ] <- c(0, sqrt(sim), sqrt(1 - sim) * ruv(K_latent - 2))
    }

    for(n in 1:N_unc) {
      loads[1 + 2*N_comp + n, ] <- c(0, 0 , ruv(K_latent - 2))
    }

    lat <- loads %*% facs
    intercepts <- c(
      1, rnorm(N_comp, mean = 1, sd = 0), 
      rnorm(N_comp, mean = -5, sd = 0),
      rnorm(N_unc, mean = 0, sd = 1)
    )
    
    Y <- lat + intercepts + rnorm(N_units * T_times, sd = 0.02)

    return(t(Y))

}

y_example <- sim_model_intercepts(sim = 0.75)

# Step 2. Generate prior predictive simulations from the model with unit-level intercepts.
# Plot these to demonstrate (a) compatibility with the observed data, and (b) the difficulty
# of assessing the decoupling between location and correlation.

test_data <- sample_model(N_units = 200, T_times = 100, K_latent = 6,
                          overall_scales = rep(1, 200), err_scale = 0.2,
                          autocor_a = 99, autocor_b = 1,
                          include_ints = TRUE,
                          nonstationary = FALSE, num_treated = 0,
                          int_scale = 1.5,
                          type = "prior_pred", quiet = FALSE)

plot_ppd <- plot_data_highlight(test_data, use_exp = FALSE, cor_perc = 0.95, num_samples = 10)
ggsave(plot_ppd, file="ppd_intercepts.pdf", device = "pdf", width = 7, height = 4)

test_data_small <- sample_model(N_units = 6, T_times = 30, K_latent = 6,
                          overall_scales = rep(1, 6), err_scale = 0.2,
                          autocor_a = 99, autocor_b = 1,
                          include_ints = TRUE,
                          nonstationary = FALSE, num_treated = 0,
                          int_scale = 1,
                          type = "prior_pred")

plot_ppd_hsmall <- plot_data_highlight(test_data_small, use_exp = FALSE, cor_perc = 0.66, num_samples = 10)
ggsave(plot_ppd_hsmall, file="ppd_intercepts_hsmall.pdf", device = "pdf", width = 7, height = 4)
