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

nonstat_prior_data <- sample_model(
  overall_scales = SIGMA_NONSTAT * overall_scales_nonstat, alpha_diag = 20,
  err_scale = ERR_NONSTAT, absolute_error = TRUE,
  data = NULL,
  autocor_a = 8, autocor_b = 2,
  nonstationary = TRUE, num_treated = 0,
  type = "prior_pred", K_latent = 4, iter = 6000,
  n_chains = 1
)

nonstat_absr <- apply(nonstat_prior_data$ys[, , 1], 1, \(x) sqrt(summary(lm(x ~ seq(1, 20)))$r.squared))

stat_prior_data <- sample_model(
  overall_scales = SIGMA_STAT * overall_scales_stat, alpha_diag = 20,
  err_scale = ERR_STAT, absolute_error = TRUE,
  data = NULL,
  autocor_a = 97, autocor_b = 3,
  nonstationary = FALSE, num_treated = 0,
  type = "prior_pred", K_latent = 4, iter = 6000,
  n_chains = 1
)

stat_absr <- apply(stat_prior_data$ys[, , 1], 1, \(x) sqrt(summary(lm(x ~ seq(1, 20)))$r.squared))

stat_prior_data1 <- sample_model(
  overall_scales = SIGMA_STAT * overall_scales_stat, alpha_diag = 20,
  err_scale = ERR_STAT, absolute_error = TRUE,
  data = NULL,
  autocor_a = 97, autocor_b = 3,
  nonstationary = FALSE, num_treated = 0,
  type = "prior_pred", K_latent = 4, iter = 6000,
  n_chains = 1
)

stat_strong_absr <- apply(stat_prior_data1$ys[, , 1], 1, \(x) sqrt(summary(lm(x ~ seq(1, 20)))$r.squared))

weak_prior_data <- sample_model(
  overall_scales = SIGMA_STAT * overall_scales_stat, alpha_diag = 20,
  err_scale = ERR_STAT, absolute_error = TRUE,
  data = NULL,
  autocor_a = 97, autocor_b = 3,
  nonstationary = FALSE, num_treated = 0,
  type = "prior_pred", K_latent = 4, iter = 6000,
  n_chains = 1
)

weak_absr <- apply(weak_prior_data$ys[, , 1], 1, \(x) sqrt(summary(lm(x ~ seq(1, 20)))$r.squared))

# WARNING -- this figure's premise no longer matches ex1_sim_study. It contrasts three STRENGTHS of
# prior on the error scale, but the study now FIXES the error scale (err_sd = 1.0, absolute) for both
# stationary arms, which differ only in K_latent (K vs K + 1). Matching the study therefore makes the
# three stationary arms below identical, and the figure degenerate. Three ways out, all decisions for
# the paper rather than the code:
#   (a) contrast K_latent 4 vs 5, mirroring the study -- but then it is not a prior-strength figure;
#   (b) keep it a prior-strength figure by varying err_sd (0.5 / 1.0 / 2.0), illustrating the dial
#       even though the study fixes it at 1.0 -- honest only if labelled as an illustration;
#   (c) drop it, if the paper no longer makes a prior-strength claim for ex1.
absrs <- list(
  `Nonstationary` = nonstat_absr,
  `Stronger Prior` = stat_strong_absr,
  `Weaker Prior` = stat_absr,
  `Vague Prior` = weak_absr
)

absr_hists <- plot_prior_absr(absrs)
ggsave(absr_hists, file = "../figs/stat_rsq_hists.pdf", device = "pdf", width = 2.5, height = 5, create.dir = TRUE)
