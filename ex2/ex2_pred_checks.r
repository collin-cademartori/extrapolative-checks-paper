#################################
## Per-Unit Intercepts Example ##
#################################

source("../sample_model.r")
source("../plotting.r")

# Reproducibility: seed the base-R RNG (used to pick the highlighted samples) and
# derive fixed Stan seeds for the two prior predictive draws up front, before any
# sample_model() call (cmdstanr's $sample() advances R's RNG).
seed <- 60412
set.seed(seed)
pp_seed_big <- sample.int(.Machine$integer.max, 1)
pp_seed_small <- sample.int(.Machine$integer.max, 1)

# Generate prior predictive simulations from the model with unit-level intercepts.
# Plot these to demonstrate (a) compatibility with the observed data, and (b) the difficulty
# of assessing the decoupling between location and correlation.

test_data <- sample_model(
  N_units = 200, T_times = 100, K_latent = 3,
  overall_scales = rep(1, 200), err_scale = 0.2,
  autocor_a = 90, autocor_b = 10,
  include_ints = TRUE,
  nonstationary = FALSE, num_treated = 0,
  int_scale = 5,
  type = "prior_pred", quiet = FALSE, seed = pp_seed_big
)

plot_ppd <- plot_data_highlight(test_data, use_exp = FALSE, cor_perc = 0.95, num_samples = 10)
ggsave(plot_ppd, file = "../figs/ppd_intercepts.pdf", device = "pdf", width = 7, height = 4, create.dir = TRUE)

test_data_small <- sample_model(
  N_units = 6, T_times = 30, K_latent = 3,
  overall_scales = rep(1, 6), err_scale = 0.2,
  autocor_a = 90, autocor_b = 10,
  include_ints = TRUE,
  nonstationary = FALSE, num_treated = 0,
  int_scale = 5,
  type = "prior_pred", seed = pp_seed_small
)

plot_ppd_hsmall <- plot_data_highlight(test_data_small, use_exp = FALSE, cor_perc = 0.66, num_samples = 10)
ggsave(plot_ppd_hsmall, file = "../figs/ppd_intercepts_hsmall.pdf", device = "pdf", width = 7, height = 4, create.dir = TRUE)
