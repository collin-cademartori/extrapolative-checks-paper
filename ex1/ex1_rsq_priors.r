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

overall_scales_stat <- rep(1, 8)
overall_scales_nonstat <- rep(1, 8)

nonstat_prior_data <- sample_model(
  overall_scales = 0.5 * overall_scales_nonstat,
  err_scale = 0, err_scale_mean = 2, err_scale_sd = 2,
  data = NULL,
  autocor_a = 8, autocor_b = 2,
  nonstationary = TRUE, num_treated = 0,
  type = "prior_pred", K_latent = 7, iter = 6000
)

nonstat_absr <- apply(nonstat_prior_data$ys[, , 1], 1, \(x) sqrt(summary(lm(x ~ seq(1, 20)))$r.squared))

stat_prior_data <- sample_model(
  overall_scales = 2 * overall_scales_stat,
  err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.1,
  data = NULL,
  autocor_a = 97, autocor_b = 3,
  nonstationary = FALSE, num_treated = 0,
  type = "prior_pred", K_latent = 7, iter = 6000
)

stat_absr <- apply(stat_prior_data$ys[, , 1], 1, \(x) sqrt(summary(lm(x ~ seq(1, 20)))$r.squared))

stat_prior_data1 <- sample_model(
  overall_scales = 2 * overall_scales_stat,
  err_scale = 0.05,
  data = NULL,
  autocor_a = 97, autocor_b = 3,
  nonstationary = FALSE, num_treated = 0,
  type = "prior_pred", K_latent = 7, iter = 6000
)

stat_strong_absr <- apply(stat_prior_data1$ys[, , 1], 1, \(x) sqrt(summary(lm(x ~ seq(1, 20)))$r.squared))

weak_prior_data <- sample_model(
  overall_scales = 2 * overall_scales_stat,
  err_scale = 0, err_scale_mean = 0.5, err_scale_sd = 0.4,
  data = NULL,
  autocor_a = 97, autocor_b = 3,
  nonstationary = FALSE, num_treated = 0,
  type = "prior_pred", K_latent = 7, iter = 6000
)

weak_absr <- apply(weak_prior_data$ys[, , 1], 1, \(x) sqrt(summary(lm(x ~ seq(1, 20)))$r.squared))

absrs <- list(
  `Nonstationary` = nonstat_absr,
  `Stronger Prior` = stat_strong_absr,
  `Weaker Prior` = stat_absr,
  `Vague Prior` = weak_absr
)

absr_hists <- plot_prior_absr(absrs)
ggsave(absr_hists, file = "../figs/stat_rsq_hists.pdf", device = "pdf", width = 2.5, height = 5, create.dir = TRUE)
