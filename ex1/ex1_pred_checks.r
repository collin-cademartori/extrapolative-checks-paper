## Prior predictive figures for the stationarity example: the classic check over the study's
## T = 20 window, and the extrapolative check over a long one.

source("../sample_model.r")
source("../plotting.r")

seed <- 72385614

# Shape comes from ex1_config.r; only T_LONG is local, for the extrapolative check.
N_units <- N_UNITS
T_times <- T_TIMES
K_latent <- K_LATENT
T_LONG <- 500

# Scale constants come from ex1_config.r.
source("ex1_config.r")

rms_y <- rep(1, N_UNITS)
overall_scales_stat    <- SIGMA_MULT_STAT * rms_y
overall_scales_nonstat <- SIGMA_MULT_NONSTAT * rms_y
eta_anchor <- mean(rms_y)

stat_prior_data <- sample_model(
  N_units = N_units, T_times = T_times,
  overall_scales = overall_scales_stat, alpha_diag = ALPHA_DIAG,
  err_scale = ETA_FRAC_STAT * eta_anchor,
  data = NULL,
  autocor_a = RHO_STAT[1], autocor_b = RHO_STAT[2],
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
  overall_scales = overall_scales_nonstat, alpha_diag = ALPHA_DIAG,
  err_scale = ETA_FRAC_NONSTAT * eta_anchor,
  data = NULL,
  autocor_a = RHO_NONSTAT[1], autocor_b = RHO_NONSTAT[2],
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
  N_units = N_units, T_times = T_LONG,
  overall_scales = overall_scales_stat, alpha_diag = ALPHA_DIAG,
  err_scale = ETA_FRAC_STAT * eta_anchor,
  data = NULL,
  autocor_a = RHO_STAT[1], autocor_b = RHO_STAT[2],
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
  N_units = N_units, T_times = T_LONG,
  overall_scales = overall_scales_nonstat, alpha_diag = ALPHA_DIAG,
  err_scale = ETA_FRAC_NONSTAT * eta_anchor,
  data = NULL,
  autocor_a = RHO_NONSTAT[1], autocor_b = RHO_NONSTAT[2],
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
