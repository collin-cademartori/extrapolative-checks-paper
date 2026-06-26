## Reproduce classic and extrapolated prior checks in the stationarity example

source("../sample_model.r")

# Seed controls the prior predictive sample produced by Stan
# This seed reproduces exact figures from paper. Change to 
# vary the prior predictive sample.
seed <- 72385614

nonstat_prior_data <- sample_model(overall_scales = 3, 
                                   err_scale = 3,
                                   data = test_ys,
                                   autocor_a = 8, autocor_b = 2,
                                   nonstationary = TRUE, num_treated = 0,
                                   type = "prior_pred", K_latent = 7,
                                   seed = seed)

nonstat_plot_ppc <- plot_data_units(
  nonstat_prior_data, 
  unit = 1, samples = 14, hide_y = FALSE
)
ggsave(nonstat_plot_ppc, file="../figs/ppc_nonstat_ns.pdf", device = "pdf", width = 7, height = 4, create.dir = TRUE)

stat_plot_ppc <- plot_data_units(
  stat_prior_data,
  unit = 1, samples = 14, hide_y = FALSE
)
ggsave(stat_plot_ppc, file="../figs/ppc_stat_ns.pdf", device = "pdf", width = 7, height = 4, create.dir = TRUE)