#################################
## Per-Unit Intercepts Example ##
#################################

source("../sample_model.r")
source("../plotting.r")
source("ex2_config.r")

seed <- 60412
set.seed(seed)
pp_seed_big <- sample.int(.Machine$integer.max, 1)
pp_seed_small <- sample.int(.Machine$integer.max, 1)

# Scale and prior constants come from ex2_config.r. Some constants are relative to other scales.
# Whereas the simulation study calculates sd(y_n) from each dataset to define these scales,
# the predictive checks here set these scales to 1.

# The large panel uses more units and a longer window than the study, yielding an
# extrapolative check where the association between location and correlation is visible.
test_data <- sample_model(
  N_units = 200, T_times = 100, K_latent = K_LATENT,
  overall_scales = rep(1, 200),
  err_scale = ETA_FRAC_EX2,
  alpha_diag = ALPHA_DIAG,
  autocor_a = RHO_EX2[1], autocor_b = RHO_EX2[2],
  include_ints = TRUE,
  nonstationary = FALSE, num_treated = 0,
  int_scale = INT_FRAC * LEVEL_SPREAD_FRAC, int_loc = 0,
  type = "prior_pred", quiet = FALSE, seed = pp_seed_big
)

plot_ppd <- plot_data_highlight(test_data, use_exp = FALSE, cor_perc = 0.95, num_samples = 10)
ggsave(plot_ppd, file = "../figs/ppd_intercepts.pdf", device = "pdf", width = 7, height = 4, create.dir = TRUE)

# The small panel matches the dimensions of the data in the study: DGP_N_UNITS units over DGP_T_TIMES periods.
test_data_small <- sample_model(
  N_units = DGP_N_UNITS, T_times = DGP_T_TIMES, K_latent = K_LATENT,
  overall_scales = rep(1, DGP_N_UNITS),
  err_scale = ETA_FRAC_EX2,
  alpha_diag = ALPHA_DIAG,
  autocor_a = RHO_EX2[1], autocor_b = RHO_EX2[2],
  include_ints = TRUE,
  nonstationary = FALSE, num_treated = 0,
  int_scale = INT_FRAC * LEVEL_SPREAD_FRAC, int_loc = 0,
  type = "prior_pred", seed = pp_seed_small
)

plot_ppd_hsmall <- plot_data_highlight(test_data_small, use_exp = FALSE, cor_perc = 0.66, num_samples = 10)
ggsave(plot_ppd_hsmall, file = "../figs/ppd_intercepts_hsmall.pdf", device = "pdf", width = 7, height = 4, create.dir = TRUE)
