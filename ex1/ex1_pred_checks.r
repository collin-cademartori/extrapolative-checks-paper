## Reproduce classic and extrapolated prior checks in the stationarity example

source("../sample_model.r")
source("../plotting.r")

# Same model as ex1_sim_study: sample_model.r's default ../ife_named.stan, whose error scale is
# tau*sigma with no ||Lambda|| term. Stated explicitly so these figures cannot silently diverge from
# the model the study fits.


# This seed reproduces the paper's figures; change it to vary the prior
# predictive sample.
seed <- 72385614

N_units <- 8
T_times <- 20
K_latent <- 4

# SCALES AND ERROR MATCH ex1_sim_study EXACTLY. The study sets sigma from each dataset, so the
# expected values are used here (measured over 200 prior-predictive datasets):
#   stat    sigma = 2 x RMS(y)       -> E = 13.33 ;  err_sd = 1.0 fixed, absolute
#   nonstat sigma = 1.2 x sd(diff y) -> E =  3.56 ;  err_sd = 2.0 fixed, absolute
# absolute_error = TRUE means err_scale is the observation-error sd on the DATA's own scale, not a
# ratio to sigma. Both studies moved to this; leaving these figures on the old estimated-ratio
# parameterization would have them describe a model nobody fits.
SIGMA_STAT    <- 13.33
SIGMA_NONSTAT <-  3.56
ERR_STAT      <-  1.0
ERR_NONSTAT   <-  2.0

overall_scales_stat <- rep(1, 8)
overall_scales_nonstat <- rep(1, 8)

stat_prior_data <- sample_model(
  N_units = N_units, T_times = T_times,
  overall_scales = 2 * overall_scales_stat, alpha_diag = 20,
  err_scale = 0.1 * mean(overall_scales_stat), absolute_error = TRUE,
  data = NULL,
  autocor_a = 98, autocor_b = 2,
  nonstationary = FALSE, num_treated = 0,
  type = "prior_pred", K_latent = K_latent,
  seed = seed, quiet = FALSE,
  parallel_chains = 4
)

cat(paste0("E(SD): ", mean(apply(stat_prior_data$ys[, , 1], 1, sd)), "\n"))

stat_plot_ppc <- plot_data_units(
  stat_prior_data,
  unit = 1, samples = 14, hide_y = TRUE
)
ggsave(stat_plot_ppc, file = "../figs/ppc_stat_ns.pdf", device = "pdf", width = 7, height = 4, create.dir = TRUE)

nonstat_prior_data <- sample_model(
  N_units = N_units, T_times = T_times,
  overall_scales = (1/7) * overall_scales_nonstat, alpha_diag = 20,
  err_scale = 2 * mean(overall_scales_nonstat), absolute_error = TRUE,
  data = NULL,
  autocor_a = 8, autocor_b = 2,
  nonstationary = TRUE, num_treated = 0,
  type = "prior_pred", iter = 1000, quiet = FALSE,
  seed = seed,
  parallel_chains = 4
)

cat(paste0("E(SD): ", mean(apply(nonstat_prior_data$ys[, , 1], 1, sd)), "\n"))

nonstat_plot_ppc <- plot_data_units(
  nonstat_prior_data,
  unit = 1, samples = 14
)
ggsave(nonstat_plot_ppc, file = "../figs/ppc_nonstat.pdf", device = "pdf", width = 7, height = 4, create.dir = TRUE)

stat_prior_data_long <- sample_model(
  N_units = N_units, T_times = 500,
 overall_scales = 2 * overall_scales_stat, alpha_diag = 20,
  err_scale = 0.1 * mean(overall_scales_stat), absolute_error = TRUE,
  data = NULL,
  autocor_a = 98, autocor_b = 2,
  nonstationary = FALSE, num_treated = 0,
  type = "prior_pred", iter = 1000, quiet = FALSE,
  seed = seed,
  parallel_chains = 4
)

stat_plot_ppc_long <- plot_data_units(
  stat_prior_data_long,
  unit = 1, samples = 14, n_x_breaks = 3
)
ggsave(stat_plot_ppc_long, file = "../figs/ppc_stat_long.pdf", device = "pdf", width = 7, height = 4, create.dir = TRUE)

nonstat_prior_data_long <- sample_model(
  N_units = N_units, T_times = 500,
 overall_scales = overall_scales_nonstat, alpha_diag = 20,
  err_scale = 2 * mean(overall_scales_nonstat), absolute_error = TRUE,
  data = NULL,
  autocor_a = 8, autocor_b = 2,
  nonstationary = TRUE, num_treated = 0,
  type = "prior_pred", iter = 1000, quiet = FALSE,
  seed = seed,
  parallel_chains = 4
)

nonstat_plot_ppc_long <- plot_data_units(
  nonstat_prior_data_long,
  unit = 1, samples = 14, n_x_breaks = 3
)
ggsave(nonstat_plot_ppc_long, file = "../figs/ppc_nonstat_long.pdf", device = "pdf", width = 7, height = 4, create.dir = TRUE)
