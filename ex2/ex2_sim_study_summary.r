library(dplyr)
library(forcats)
library(ggplot2)
library(tidyr)

source("../plotting.r")
load("sim_study_ints.RData")

# 99% posterior-predictive interval coverage, averaged per condition. Both models
# cover in excess of 99%, confirming the check cannot rule out the intercepts
# model (paper Section 5). This is a purely numeric result -- the standardized
# error and correlation summaries are instead conveyed by the plots below.
perc_summary <- sim_study_ints |>
  group_by(num_comp, sim) |>
  summarize(
    q5_perc_no_ints = quantile(no_ints_pred_perc, 0.05),
    q95_perc_no_ints = quantile(no_ints_pred_perc, 0.95),
    q5_perc_ints = quantile(ints_pred_perc, 0.05),
    q95_perc_ints = quantile(ints_pred_perc, 0.95),
    mean_perc_ints = mean(ints_pred_perc),
    mean_perc_no_ints = mean(no_ints_pred_perc),
    .groups = "drop"    
  )

cat("\n99% posterior-predictive interval coverage by condition (no-int vs with-int):\n")
print(perc_summary, width = Inf)

## Statistic S2 predictive p-value (see paper Section 5), averaged separately for
## each data generating condition (num_comp x sim).
loc_cor_summary <- sim_study_ints |>
  group_by(num_comp, sim) |>
  summarize(
    q5_loc_cor_no_ints = quantile(no_ints_loc_cor_pval, 0.05),
    q95_loc_cor_no_ints = quantile(no_ints_loc_cor_pval, 0.95),
    q5_loc_cor_ints = mean(ints_loc_cor_pval, 0.05),
    q95_loc_cor_ints = mean(ints_loc_cor_pval, 0.95),
    .groups = "drop"
  )

cat("\nS2 location-correlation predictive p-value by condition (no-int vs with-int):\n")
print(loc_cor_summary, width = Inf)

# Per-condition mean and +/-2 SE bands over post-treatment time for a per-time
# statistic, for the no-intercepts and with-intercepts models. stat is "absz"
# (standardized error) or "mean" (the posterior mean, whose absolute value is the
# error since the true effect is 0).
summarize_error <- function(stat) {
  abs(sim_study_ints) |>
    pivot_longer(
      cols = contains(paste0(stat, "_")),
      names_to = c(".value", "time"),
      names_transform = list(time = as.integer),
      names_pattern = "(.*)_(\\d+)$"
    ) |>
    mutate(
      sim_f = as.factor(paste0("l == ", sim)),
      num_f = as.factor(paste0("b == ", num_comp))
    ) |>
    group_by(time, sim_f, num_f) |>
    summarize(
      ni_mean = mean(.data[[paste0("no_ints_", stat)]]),
      ni_se   = sd(.data[[paste0("no_ints_", stat)]]) / sqrt(n()),
      it_mean = mean(.data[[paste0("ints_", stat)]]),
      it_se   = sd(.data[[paste0("ints_", stat)]]) / sqrt(n()),
      .groups = "drop"
    ) |>
    mutate(
      ni_lower = ni_mean - 2 * ni_se,
      ni_upper = ni_mean + 2 * ni_se,
      it_lower = it_mean - 2 * it_se,
      it_upper = it_mean + 2 * it_se
    )
}

# Per-condition time series of the two models' means with shaded +/-2 SE bands
# (no-intercepts solid, with-intercepts dashed).
plot_error_bands <- function(df, y_label) {
  ggplot(data = df) +
    geom_ribbon(aes(x = time, ymin = ni_lower, ymax = ni_upper), alpha = 0.2, fill = "#858585") +
    geom_ribbon(aes(x = time, ymin = it_lower, ymax = it_upper), alpha = 0.2, fill = "#858585") +
    geom_line(aes(x = time, y = ni_mean)) +
    geom_line(aes(x = time, y = it_mean), linetype = "dashed") +
    facet_grid(vars(num_f), vars(sim_f), scales = "free_y") +
    xlab("Post-Treatment Time") +
    ylab(y_label) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(strip.background = element_rect(fill = "white", color = "black"))
}

sim_study_std_err <- summarize_error("absz")
std_err_plot <- plot_error_bands(
  sim_study_std_err, "Average Standardized Error of Posterior\n Expected Treatment Effect"
)
ggsave(std_err_plot, device = "pdf", width = 5, height = 4, file = "../figs/ints_std_err.pdf", create.dir = TRUE)

sim_study_abs_err <- summarize_error("mean")
abs_err_plot <- plot_error_bands(
  sim_study_abs_err, "Average Absolute Error of Posterior\n Expected Treatment Effect"
)
ggsave(abs_err_plot, device = "pdf", width = 5, height = 4, file = "../figs/ints_abs_err.pdf", create.dir = TRUE)

## Overfitting plot

sim_study_overfit <- abs(sim_study_ints) |>
  pivot_longer(
    cols = !c(num_comp, sim),
    names_to = c("model", ".value"),
    names_transform = list(time = as.integer),
    names_pattern = "^(no_ints|ints)_(.*)$"
  ) |>
  filter(num_comp == 3) |>
  mutate(
    errspur = (cor_sq_3 + cor_sq_4 + cor_sq_5) / 3
  ) |>
  select(!starts_with("cor_sq")) |>
  # Second pivot only over the time-indexed columns (name ends in _<digit>);
  # the scalar per-model stats (pred_perc, pred_mad, loc_cor_pval) are left
  # untouched, so new scalar stats can be added without breaking this pivot.
  pivot_longer(
    cols = matches("_\\d+$"),
    names_to = c(".value", "time"),
    names_transform = list(time = as.integer),
    names_pattern = "(.*)_(\\d+)$"
  ) |>
  mutate(
    model = fct_recode(as.factor(model),
      `No Int` = "no_ints",
      `With Int` = "ints"
    ),
    time = paste0("Time ", time)
  ) |>
  group_by(time, model, sim) |>
  summarize(
    mean_absz = mean(absz),
    mean_err = mean(errspur),
    .groups = "drop"
  )

overfit_plot <- ggplot(data = sim_study_overfit) +
  geom_line(aes(x = mean_err, y = mean_absz), linewidth = 0.8) +
  geom_label(aes(label = model, x = mean_err, y = mean_absz), size = 3) +
  facet_grid(vars(sim), vars(time), scales = "free") +
  xlab("Modeled Long-Run Correlation (Treated vs Spuriously Correlated Units)") +
  ylab("Average Standardized Error of Posterior\n Expected Treatment Effect") +
  scale_x_continuous(expand = expansion(mult = 0.6), n.breaks = 3) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(strip.background = element_rect(fill = "white", color = "black"))

ggsave(overfit_plot, device = "pdf", width = 6, height = 3.5, file = "../figs/ints_overfit.pdf", create.dir = TRUE)
