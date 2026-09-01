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

# SCALES AND ERROR MATCH ex1_sim_study. Same names, same multiples -- see the long comment in
# ex1_sim_study.r for what each constant is and why it takes the value it does.
#
# The study anchors every scale on RMS(y_n), measured from each dataset. These figures have no
# dataset to measure, so they work in units where RMS(y) = 1 and carry the study's multiples
# unchanged. Both things these figures show -- the shape of a series, and |cor(y, t)| -- are scale
# invariant, so the choice of unit is a y-axis relabeling and nothing more.
#
# absolute_error = TRUE means err_scale is the observation-error sd on the DATA's own scale, not a
# ratio to sigma, and it is the LEVEL error sd in both the stationary and the differenced branch
# (ife_named.stan applies the sqrt(2) differencing inflation itself). So the ETA_FRAC_* constants
# below take no differencing correction, and a caller-side one would double-count.
#
# One consequence of carrying the study's multiples is worth knowing before reading the panels
# against each other: the arms do NOT come out on a common scale. Measured at these settings,
# E[RMS(y)] is 0.93 for nonstat but 1.62 for the stationary arms, because SIGMA_MULT_STAT = 2 is
# deliberately ~1.6x generous (the stationary self-consistency fixed point is 1.22). Their
# dispersion does match closely -- E[sd(y)] is 0.56 against 0.62 -- since the gap is level, not
# spread. Inherited from the study by design, not a defect here.
# Scale constants come from ex1_config.r, the same file the study and ex1_derive_scales.r read, so
# these figures cannot describe a model nobody fits. Only the ANCHOR is local: the study measures
# RMS(y_n) from each dataset, and this file has no dataset, so it works in units where the anchor is
# 1 and carries the study's multiples unchanged. Both things these figures show -- the shape of a
# series and |cor(y, t)| -- are scale invariant, so the unit is a y-axis relabeling and nothing more.
source("ex1_config.r")

rms_y <- rep(1, N_UNITS)
overall_scales_stat    <- SIGMA_MULT_STAT * rms_y
overall_scales_nonstat <- SIGMA_MULT_NONSTAT * rms_y
eta_anchor <- mean(rms_y)

stat_prior_data <- sample_model(
  N_units = N_units, T_times = T_times,
  overall_scales = overall_scales_stat, alpha_diag = 20,
  err_scale = ETA_FRAC_STAT * eta_anchor, absolute_error = TRUE,
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
  overall_scales = overall_scales_nonstat, alpha_diag = 20,
  err_scale = ETA_FRAC_NONSTAT * eta_anchor, absolute_error = TRUE,
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
  overall_scales = overall_scales_stat, alpha_diag = 20,
  err_scale = ETA_FRAC_STAT * eta_anchor, absolute_error = TRUE,
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
  err_scale = ETA_FRAC_NONSTAT * eta_anchor, absolute_error = TRUE,
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
