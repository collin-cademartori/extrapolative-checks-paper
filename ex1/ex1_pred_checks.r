## Reproduce classic and extrapolated prior checks in the stationarity example

source("../sample_model.r")
source("../plotting.r")

# Use the same model as ex1_sim_study: the Cartesian variant, whose error scale is anchored to the
# ESTIMATED signal scale sigma_data[n]*||Lambda[n,:]|| rather than to sigma_data alone. That changes
# the prior predictive (the error co-scales with the loadings), so these figures must be generated
# from the same model the study fits, not from sample_model.r's default ../ife_named.stan.
ife_mod <- cmdstan_model(stan_file = "../ife_named_cartesian.stan")
stopifnot(basename(ife_mod$stan_file()) == "ife_named_cartesian.stan")


# This seed reproduces the paper's figures; change it to vary the prior
# predictive sample.
seed <- 72385614

N_units <- 8
T_times <- 20
K_latent <- 7

stat_prior_data <- sample_model(
  N_units = N_units, T_times = T_times,
  overall_scales = rep(1, N_units), alpha_diag = 10,
  err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.1,
  data = NULL,
  autocor_a = 97, autocor_b = 3,
  nonstationary = FALSE, num_treated = 0,
  type = "prior_pred", K_latent = K_latent,
  seed = seed, quiet = FALSE
)

stat_plot_ppc <- plot_data_units(
  stat_prior_data,
  unit = 1, samples = 14, hide_y = TRUE
)
ggsave(stat_plot_ppc, file = "../figs/ppc_stat_ns.pdf", device = "pdf", width = 7, height = 4, create.dir = TRUE)


stat_prior_data_long <- sample_model(
  N_units = N_units, T_times = 500,
  overall_scales = rep(1, N_units), alpha_diag = 10,
  err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.1,
  data = NULL,
  autocor_a = 97, autocor_b = 3,
  nonstationary = FALSE, num_treated = 0,
  type = "prior_pred", iter = 1000, quiet = FALSE,
  seed = seed
)

stat_plot_ppc_long <- plot_data_units(
  stat_prior_data_long,
  unit = 1, samples = 14, n_x_breaks = 3
)
ggsave(stat_plot_ppc_long, file = "../figs/ppc_stat_long.pdf", device = "pdf", width = 7, height = 4, create.dir = TRUE)

nonstat_prior_data <- sample_model(
  N_units = N_units, T_times = T_times,
  overall_scales = rep(1, N_units), alpha_diag = 10,
  err_scale = 0, err_scale_mean = 2, err_scale_sd = 2,
  data = NULL,
  autocor_a = 8, autocor_b = 2,
  nonstationary = TRUE, num_treated = 0,
  type = "prior_pred", iter = 1000, quiet = FALSE,
  seed = seed
)

nonstat_plot_ppc <- plot_data_units(
  nonstat_prior_data,
  unit = 1, samples = 14
)
ggsave(nonstat_plot_ppc, file = "../figs/ppc_nonstat.pdf", device = "pdf", width = 7, height = 4, create.dir = TRUE)

nonstat_prior_data_long <- sample_model(
  N_units = N_units, T_times = 500,
  overall_scales = rep(1, N_units), alpha_diag = 10,
  err_scale = 0, err_scale_mean = 2, err_scale_sd = 2,
  data = NULL,
  autocor_a = 8, autocor_b = 2,
  nonstationary = TRUE, num_treated = 0,
  type = "prior_pred", iter = 1000, quiet = FALSE,
  seed = seed
)

nonstat_plot_ppc_long <- plot_data_units(
  nonstat_prior_data_long,
  unit = 1, samples = 14, n_x_breaks = 3
)
ggsave(nonstat_plot_ppc_long, file = "../figs/ppc_nonstat_long.pdf", device = "pdf", width = 7, height = 4, create.dir = TRUE)
