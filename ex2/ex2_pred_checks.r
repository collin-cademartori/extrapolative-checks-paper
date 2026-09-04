#################################
## Per-Unit Intercepts Example ##
#################################

source("../sample_model.r")
source("../plotting.r")
source("ex2_config.r")

# Stan seeds are drawn up front, before any sample_model() call, because cmdstanr's $sample()
# advances R's RNG.
seed <- 60412
set.seed(seed)
pp_seed_big <- sample.int(.Machine$integer.max, 1)
pp_seed_small <- sample.int(.Machine$integer.max, 1)

# Scale and prior constants come from ex2_config.r, the same file the study and
# ex2_derive_scales.r read, so this figure cannot describe a model nobody fits. Only the ANCHOR is
# local: the study measures sd(y_n) from each dataset and this file has no dataset, so it works in
# units where that anchor is 1 and carries the study's multiples unchanged.
#
# The figure contains no data, by design: the check asks whether the configuration implies what the
# analyst expects, before any data exists. eta is fixed here at its prior location, where the study
# estimates it, so each panel is a statement about one number.

# Generate prior predictive simulations from the model with unit-level intercepts.
# Plot these to demonstrate (a) compatibility with the observed data, and (b) the difficulty
# of assessing the decoupling between location and correlation.

# The large panel uses more units and a longer window than the study, at a scale where the
# location/correlation entanglement is visible. Everything else is the study's configuration.
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

# The small panel matches the study exactly: DGP_N_UNITS units over DGP_T_TIMES periods.
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
