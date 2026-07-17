## This file generates samples from the prior predictive distribution of the nonstationary
## model and the stationary model with several different priors on the iid error scale.
## Produces plot of histograms comparing implied predictive distribution on R squared
## for predicting observed series from time. Distributions favoring larger R squared
## values more frequently produce series with approximately linear trends.

source("../sample_model.r")
source("../plotting.r")

plot_prior_rsq <- function(rsqs) {
  rsq_df <- stack(rsqs)
  names(rsq_df) <- c("rsq", "model")
  
  prior_plot <- ggplot(rsq_df) +
    facet_wrap(vars(model), ncol = 1, scales = "free_y", strip.position = "left") +
    geom_histogram(aes(x=rsq), fill = "#aaaaaa", color = "#aaaaaa") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(strip.background = element_rect(fill = "white", color = "black")) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    xlab("Absolute Correlation") +
    ylab("") +
    ggtitle("(A)")

  return(prior_plot)
    
}

overall_scales_stat <- rep(1,8)
overall_scales_nonstat <- rep(1,8)

nonstat_prior_data <- sample_model(overall_scales = 0.5 * overall_scales_nonstat, 
                                   #err_scale = 3,
                                   err_scale = 0, err_scale_mean = 2, err_scale_sd = 2,
                                   data = NULL,
                                   autocor_a = 8, autocor_b = 2,
                                   nonstationary = TRUE, num_treated = 0,
                                   type = "prior_pred", K_latent = 7, iter = 4000)

nonstat_rsq <- apply(nonstat_prior_data$ys[,,1], 1, \(x) summary(lm(x ~ seq(1, 20)))$r.squared)

stat_prior_data <- sample_model(overall_scales = 2 * overall_scales_stat, 
                                err_scale = 0, err_scale_mean = 0.2, err_scale_sd = 0.2,
                                data = NULL,
                                autocor_a = 99, autocor_b = 1,
                                nonstationary = FALSE, num_treated = 0,
                                type = "prior_pred", K_latent = 7, iter = 4000)

stat_rsq <- apply(stat_prior_data$ys[,,1], 1, \(x) summary(lm(x ~ seq(1, 20)))$r.squared)

stat_prior_data1 <- sample_model(overall_scales = 2 * overall_scales_stat, 
                                err_scale = 0.05,
                                data = NULL,
                                autocor_a = 99, autocor_b = 1,
                                nonstationary = FALSE, num_treated = 0,
                                type = "prior_pred", K_latent = 7, iter = 4000)

stat1_rsq <- apply(stat_prior_data1$ys[,,1], 1, \(x) summary(lm(x ~ seq(1, 20)))$r.squared)

weak_prior_data <- sample_model(overall_scales = 2 * overall_scales_stat, 
                                err_scale = 0, err_scale_mean = 0.5, err_scale_sd = 0.4,
                                data = NULL,
                                autocor_a = 99, autocor_b = 1,
                                nonstationary = FALSE, num_treated = 0,
                                type = "prior_pred", K_latent = 7, iter = 4000)

weak_rsq <- apply(weak_prior_data$ys[,,1], 1, \(x) summary(lm(x ~ seq(1, 20)))$r.squared)

rsqs <- list(
  `Nonstationary` = nonstat_rsq,
  `Stronger Prior` = stat1_rsq,
  `Weaker Prior` = stat_rsq,
  `Vague Prior` = weak_rsq
)
rsqs <- lapply(rsqs, sqrt)

rsq_hists <- plot_prior_rsq(rsqs)
ggsave(rsq_hists, file = "../figs/stat_rsq_hists.pdf", device = "pdf", width = 4, height = 5, create.dir = TRUE)
