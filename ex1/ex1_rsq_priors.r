## This file generates samples from the prior predictive distribution of the nonstationary
## model and the stationary model with several different priors on the iid error scale.
## Produces a plot of histograms comparing the implied predictive distribution of the
## absolute correlation between the observed series and time. Distributions favoring
## larger absolute correlations more frequently produce series with approximately
## linear trends.

source("../sample_model.r")
source("../plotting.r")

# Same model as ex1_sim_study: sample_model.r's default ../ife_named.stan, whose error scale is
# tau*sigma with no ||Lambda|| term. Stated explicitly so these figures cannot silently diverge from
# the model the study fits.


plot_prior_absr <- function(absrs) {
  absr_df <- stack(absrs)
  names(absr_df) <- c("absr", "model")

  prior_plot <- ggplot(absr_df) +
    facet_wrap(vars(model), ncol = 1, scales = "free_y", strip.position = "left") +
    geom_histogram(aes(x = absr), fill = "#aaaaaa", color = "#aaaaaa") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(strip.background = element_rect(fill = "white", color = "black")) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    xlab("Absolute Correlation\n (Outcome vs Time)") +
    ylab("") +
    ggtitle("(A)")

  return(prior_plot)
}

# SCALES AND ERROR MATCH ex1_sim_study. Same names, same multiples -- see the long comment in
# ex1_sim_study.r for what each constant is and why it takes the value it does.
#
# The study anchors every scale on RMS(y_n), measured from each dataset. This file has no dataset to
# measure, so it works in units where RMS(y) = 1 and carries the study's multiples unchanged. What
# this figure shows -- |cor(y, t)| -- is scale invariant, so the choice of unit does not enter it.
#
# absolute_error = TRUE means err_scale is the observation-error sd on the DATA's own scale, not a
# ratio to sigma, and it is the LEVEL error sd in both the stationary and the differenced branch
# (ife_named.stan applies the sqrt(2) differencing inflation itself). So the ETA_FRAC_* constants
# below take no differencing correction, and a caller-side one would double-count.
#
# The four panels are a ladder in ETA_FRAC alone -- ETA_FRAC_NONSTAT (0.323) / 0.100 / 0.050 /
# 0.500 as a fraction of
# RMS(y) -- which is the point of the figure: a large error attenuates a linear trend, so only a
# small one lets a stationary model put enough prior mass on high |cor(y, t)| to be plausible for
# these data. The stationary arms' sigma is held fixed across the ladder so the error scale is the
# only thing that moves.
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
# The vague rung, this file only: a stationary configuration whose error scale is large enough to
# flatten any trend. It is NOT a study arm. Its role is to show what the check rules OUT -- the
# implied prior on |cor(y, t)| collapses toward zero, which is incompatible with expecting series
# that trend -- and thereby what pins ETA_FRAC_STAT at the study's much smaller value.
ETA_FRAC_VAGUE <- 0.5

nonstat_prior_data <- sample_model(
  overall_scales = overall_scales_nonstat, alpha_diag = ALPHA_DIAG,
  err_scale = ETA_FRAC_NONSTAT * eta_anchor, absolute_error = TRUE,
  data = NULL,
  autocor_a = RHO_NONSTAT[1], autocor_b = RHO_NONSTAT[2],
  nonstationary = TRUE, num_treated = 0,
  type = "prior_pred", K_latent = K_LATENT, iter = 6000,
  n_chains = 1
)

nonstat_absr <- apply(nonstat_prior_data$ys[, , 1], 1, \(x) sqrt(summary(lm(x ~ seq_len(T_TIMES)))$r.squared))

stat_prior_data <- sample_model(
  overall_scales = overall_scales_stat, alpha_diag = ALPHA_DIAG,
  err_scale = ETA_FRAC_STAT * eta_anchor, absolute_error = TRUE,
  data = NULL,
  autocor_a = RHO_STAT[1], autocor_b = RHO_STAT[2],
  nonstationary = FALSE, num_treated = 0,
  type = "prior_pred", K_latent = K_LATENT, iter = 6000,
  n_chains = 1
)

stat_absr <- apply(stat_prior_data$ys[, , 1], 1, \(x) sqrt(summary(lm(x ~ seq_len(T_TIMES)))$r.squared))


vague_prior_data <- sample_model(
  overall_scales = overall_scales_stat, alpha_diag = ALPHA_DIAG,
  err_scale = ETA_FRAC_VAGUE * eta_anchor, absolute_error = TRUE,
  data = NULL,
  autocor_a = RHO_STAT[1], autocor_b = RHO_STAT[2],
  nonstationary = FALSE, num_treated = 0,
  type = "prior_pred", K_latent = K_LATENT, iter = 6000,
  n_chains = 1
)

vague_absr <- apply(vague_prior_data$ys[, , 1], 1, \(x) sqrt(summary(lm(x ~ seq_len(T_TIMES)))$r.squared))

# Three panels: the two configurations the study actually fits, plus the vague stationary variant
# that the check rules out. Only Nonstationary and Stationary are study arms.
#
# One difference from the study remains, and it is deliberate: the study ESTIMATES eta under
# TruncNormal(ETA_FRAC * eta_anchor, ETA_FRAC * eta_anchor) while these panels FIX it at the prior
# location. Measured, the difference is negligible -- mean|r| 0.565 fixed against 0.555 estimated,
# P(|r| > 0.8) 0.269 against 0.252 -- because at eta this small against sigma = 2 the error barely
# moves the statistic either way. Fixing keeps each panel a statement about one number.
absrs <- list(
  `Nonstationary` = nonstat_absr,
  `Stationary` = stat_absr,
  `Vague Error` = vague_absr
)

absr_hists <- plot_prior_absr(absrs)
ggsave(absr_hists, file = "../figs/stat_rsq_hists.pdf", device = "pdf", width = 2.5, height = 5, create.dir = TRUE)
