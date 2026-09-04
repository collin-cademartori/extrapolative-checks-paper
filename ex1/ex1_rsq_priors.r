## This file generates samples from the prior predictive distribution of the nonstationary
## model and the stationary model with several different priors on the iid error scale.
## Produces a plot of histograms comparing the implied predictive distribution of the
## absolute correlation between the observed series and time. Distributions favoring
## larger absolute correlations more frequently produce series with approximately
## linear trends.

source("../sample_model.r")
source("../plotting.r")

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

# Scale constants come from ex1_config.r.
source("ex1_config.r")

rms_y <- rep(1, N_UNITS)
overall_scales_stat    <- SIGMA_MULT_STAT * rms_y
overall_scales_nonstat <- SIGMA_MULT_NONSTAT * rms_y
eta_anchor <- mean(rms_y)
# Constant defining the iid error scale for a stationary model with wide prior on the error.
ETA_FRAC_VAGUE <- 0.5

nonstat_prior_data <- sample_model(
  overall_scales = overall_scales_nonstat, alpha_diag = ALPHA_DIAG,
  err_scale = ETA_FRAC_NONSTAT * eta_anchor,
  data = NULL,
  autocor_a = RHO_NONSTAT[1], autocor_b = RHO_NONSTAT[2],
  nonstationary = TRUE, num_treated = 0,
  type = "prior_pred", K_latent = K_LATENT, iter = 6000,
  n_chains = 1
)

nonstat_absr <- apply(nonstat_prior_data$ys[, , 1], 1, \(x) sqrt(summary(lm(x ~ seq_len(T_TIMES)))$r.squared))

stat_prior_data <- sample_model(
  overall_scales = overall_scales_stat, alpha_diag = ALPHA_DIAG,
  err_scale = ETA_FRAC_STAT * eta_anchor,
  data = NULL,
  autocor_a = RHO_STAT[1], autocor_b = RHO_STAT[2],
  nonstationary = FALSE, num_treated = 0,
  type = "prior_pred", K_latent = K_LATENT, iter = 6000,
  n_chains = 1
)

stat_absr <- apply(stat_prior_data$ys[, , 1], 1, \(x) sqrt(summary(lm(x ~ seq_len(T_TIMES)))$r.squared))


vague_prior_data <- sample_model(
  overall_scales = overall_scales_stat, alpha_diag = ALPHA_DIAG,
  err_scale = ETA_FRAC_VAGUE * eta_anchor,
  data = NULL,
  autocor_a = RHO_STAT[1], autocor_b = RHO_STAT[2],
  nonstationary = FALSE, num_treated = 0,
  type = "prior_pred", K_latent = K_LATENT, iter = 6000,
  n_chains = 1
)

vague_absr <- apply(vague_prior_data$ys[, , 1], 1, \(x) sqrt(summary(lm(x ~ seq_len(T_TIMES)))$r.squared))

# Three panels: the two configurations used in simulation stidy, plus the stationary variant with vague
# error prior.
absrs <- list(
  `Nonstationary` = nonstat_absr,
  `Stationary` = stat_absr,
  `Vague Error` = vague_absr
)

absr_hists <- plot_prior_absr(absrs)
ggsave(absr_hists, file = "../figs/stat_rsq_hists.pdf", device = "pdf", width = 2.5, height = 5, create.dir = TRUE)
