#################################
## Per-Unit Intercepts Example ##
#################################

source("../sample_model.r")
source("../plotting.r")
source("ex2_config.r")

# Reproducibility: seed the base-R RNG (used to pick the highlighted samples) and
# derive fixed Stan seeds for the two prior predictive draws up front, before any
# sample_model() call (cmdstanr's $sample() advances R's RNG).
seed <- 60412
set.seed(seed)
pp_seed_big <- sample.int(.Machine$integer.max, 1)
pp_seed_small <- sample.int(.Machine$integer.max, 1)

# SCALES AND PRIORS MATCH ex2_sim_study. Every constant comes from ex2_config.r, the same file the
# study and ex2_derive_scales.r read, so this figure cannot describe a model nobody fits. It did:
# it drew intercepts from N(0, 5) where the study fits an anchored prior, ran at
# alpha_diag = 0 against the study's 20, and used an error scale of 0.2 where these units call for
# ETA_FRAC_EX2 = 0.12.
#
# UNITS. The study sets the `ints` arm's sigma to sd(y_n), measured per dataset. This file has no
# dataset, so it works in units where that anchor is 1 -- overall_scales = 1 -- and carries the
# study's multiples unchanged. ETA_FRAC_EX2 is then read directly as the error sd, and it matches
# the DGP's own noise-to-dispersion ratio (0.201 / 1.66 = 0.121) as it should.
#
# The intercept prior follows the same one rule as the study -- location = the grand mean, scale =
# INT_FRAC * sd(colMeans(y)) -- expressed in these units: the grand mean is 0 because the units are
# centred, and the level spread is LEVEL_SPREAD_FRAC times the dispersion anchor.
#
# Like ex1's prior predictive figures this one contains NO DATA, by design: the check asks whether
# the configuration implies what the analyst expects, before any data exists.
#
# One deliberate difference from the study: eta is FIXED here at the prior location where the study
# ESTIMATES it under a truncated normal. Fixing keeps each panel a statement about one number; on
# ex1 the same simplification moved the plotted statistic by under 2%.

# Generate prior predictive simulations from the model with unit-level intercepts.
# Plot these to demonstrate (a) compatibility with the observed data, and (b) the difficulty
# of assessing the decoupling between location and correlation.

# The large panel deliberately uses more units and a longer window than the study, to show the
# prior's behaviour at a scale where the location/correlation entanglement is visible. Everything
# else is the study's configuration.
test_data <- sample_model(
  N_units = 200, T_times = 100, K_latent = K_LATENT,
  overall_scales = rep(1, 200),
  err_scale = ETA_FRAC_EX2, absolute_error = TRUE,
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
  err_scale = ETA_FRAC_EX2, absolute_error = TRUE,
  alpha_diag = ALPHA_DIAG,
  autocor_a = RHO_EX2[1], autocor_b = RHO_EX2[2],
  include_ints = TRUE,
  nonstationary = FALSE, num_treated = 0,
  int_scale = INT_FRAC * LEVEL_SPREAD_FRAC, int_loc = 0,
  type = "prior_pred", seed = pp_seed_small
)

plot_ppd_hsmall <- plot_data_highlight(test_data_small, use_exp = FALSE, cor_perc = 0.66, num_samples = 10)
ggsave(plot_ppd_hsmall, file = "../figs/ppd_intercepts_hsmall.pdf", device = "pdf", width = 7, height = 4, create.dir = TRUE)
