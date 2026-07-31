## This file loads the result of the simulation study for the nonstationary example
## and generates summary plots.
## Plots summarize (i) absolute and (ii) standardized error against time, with
## bands to show standard errors. Plot (iii) shows standardized error against
## predictive error for posterior predictive mean, with models labelled.

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

source("../plotting.r")
load("sim_study_ns.RData")

# Numeric-only results (the standardized-error curves and overfit trade-off are
# shown in the plots below): 99% posterior-predictive interval coverage -- at
# least 99% for every model, so the check cannot rule out the misspecified model
# (paper Section 5) -- and the S1 time-correlation predictive p-value, averaged
# over the study for each model.
models <- c("nonstat", "stat_weak", "stat_strong")
num_summary <- cbind(
  coverage   = colMeans(sim_study_stat[paste0(models, "_pred_perc")]),
  `S1 p-val 0.05Q` = apply(sim_study_stat[paste0(models, "_time_cor_pval")], 2, \(x) quantile(x, 0.05)),
  `S1 p-val 0.95Q` = apply(sim_study_stat[paste0(models, "_time_cor_pval")], 2, \(x) quantile(x, 0.95))
)
rownames(num_summary) <- models

cat("\nPer-model 99% interval coverage and S1 time-correlation predictive p-value:\n")
print(round(num_summary, 3))

# Build the per-time mean and +/-2 SE bands for the nonstationary ("ns") and
# weaker-prior stationary ("st", i.e. stat_weak) models, for a given per-time
# statistic: the signed error "mean" or the standardized error "absz".
summarize_error <- function(stat) {
  abs(sim_study_stat) |>
    select(contains(paste0(stat, "_"))) |>
    pivot_longer(
      cols = everything(),
      names_to = c(".value", "time"),
      names_transform = list(time = as.integer),
      names_pattern = "(.*)_(\\d+)$"
    ) |>
    group_by(time) |>
    summarize(
      ns_mean = mean(.data[[paste0("nonstat_", stat)]]),
      ns_se   = sd(.data[[paste0("nonstat_", stat)]]) / sqrt(n()),
      st_mean = mean(.data[[paste0("stat_weak_", stat)]]),
      st_se   = sd(.data[[paste0("stat_weak_", stat)]]) / sqrt(n()),
      .groups = "drop"
    ) |>
    mutate(
      ns_lower = ns_mean - 2 * ns_se,
      ns_upper = ns_mean + 2 * ns_se,
      st_lower = st_mean - 2 * st_se,
      st_upper = st_mean + 2 * st_se
    )
}

# Time series of the two models' means with shaded +/-2 SE bands
# (nonstationary solid, stationary dashed).
plot_error_bands <- function(df, y_label, label) {
  ggplot(data = df) +
    geom_ribbon(
      aes(x = time, ymin = ns_lower, ymax = ns_upper),
      alpha = 0.2, fill = "#858585"
    ) +
    geom_ribbon(
      aes(x = time, ymin = st_lower, ymax = st_upper),
      alpha = 0.2, fill = "#858585"
    ) +
    geom_line(aes(x = time, y = ns_mean)) +
    geom_line(aes(x = time, y = st_mean), linetype = "dashed") +
    xlab("Post-Treatment Time") +
    ylab(y_label) +
    ggtitle(label) +
    theme_bw()
}

sim_study_abs_err <- summarize_error("mean")
abs_mad_plot <- plot_error_bands(
  sim_study_abs_err, "Average Absolute Error of Posterior\n Expected Treatment Effect", "(A)"
)
ggsave(abs_mad_plot, device = "pdf", width = 5, height = 4, file = "../figs/stat_abs_err.pdf", create.dir = TRUE)

sim_study_std_err <- summarize_error("absz")
std_err_plot <- plot_error_bands(
  sim_study_std_err, "Average Standardized Error of Posterior\n Expected Treatment Effect", "(B)"
)
ggsave(std_err_plot, device = "pdf", width = 5, height = 4, file = "../figs/stat_std_err.pdf", create.dir = TRUE)

sim_study_overfit <- abs(sim_study_stat) |>
  pivot_longer(
    cols = everything(),
    names_to = c("model", ".value"),
    names_transform = list(time = as.integer),
    names_pattern = "^(nonstat|stat_weak|stat_strong)_(.*)$"
  ) |>
  # Second pivot only over the time-indexed columns (name ends in _<digit>);
  # the scalar per-model stats (pred_mad, pred_perc, ...) are left untouched, so
  # new scalar stats can be added without breaking this pivot.
  pivot_longer(
    cols = matches("_\\d+$"),
    names_to = c(".value", "time"),
    names_transform = list(time = as.integer),
    names_pattern = "(.*)_(\\d+)$"
  ) |>
  group_by(time, model) |>
  summarize(
    mean_absz = mean(absz),
    mean_mad = mean(pred_mad),
    .groups = "drop"
  ) |>
  mutate(
    model = fct_recode(as.factor(model),
      "Nonstat" = "nonstat",
      "Weaker" = "stat_weak",
      "Stronger" = "stat_strong"
    ),
    time = paste0("Time ", time)
  )

overfit_plot <- ggplot(data = sim_study_overfit) +
  geom_line(aes(x = mean_mad, y = mean_absz, group = time), linewidth = 0.8) +
  geom_label(aes(label = model, x = mean_mad, y = mean_absz), size = 3) +
  facet_wrap(vars(time), ncol = 1, strip.position = "right") +
  xlab("Average MAD of Posterior\n Predictive Expectation") +
  ylab("Average Standardized Error of Posterior Expected Treatment Effect") +
  scale_x_continuous(expand = expansion(mult = 0.3)) +
  scale_y_continuous(expand = expansion(mult = 0.2)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(strip.background = element_rect(fill = "white", color = "black")) +
  ggtitle("(B)")

ggsave(overfit_plot, device = "pdf", width = 2.5, height = 5, file = "../figs/stat_overfit.pdf", create.dir = TRUE)
