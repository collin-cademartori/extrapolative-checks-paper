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
    mean_perc_nint = mean(nint_pred_perc),
    mean_perc_ints = mean(ints_pred_perc),
    .groups = "drop"
  )

cat("\n99% posterior-predictive interval coverage by condition (no-int vs with-int):\n")
print(perc_summary)

## Statistic S2 predictive p-value (see paper Section 5), averaged separately for
## each data generating condition (num_comp x sim).
loc_cor_summary <- sim_study_ints |>
  group_by(num_comp, sim) |>
  summarize(
    mean_loc_cor_nint = mean(nint_loc_cor_pval),
    mean_loc_cor_ints = mean(ints_loc_cor_pval),
    .groups = "drop"
  )

cat("\nS2 location-correlation predictive p-value by condition (no-int vs with-int):\n")
print(loc_cor_summary)

sim_study_std_err <- abs(sim_study_ints) |>
  # select(contains("absz_")) |>
  pivot_longer(
    cols = contains("absz_"),
    names_to = c(".value", "time"),
    names_transform = list(time = as.integer),
    names_pattern = "(.*)_(\\d+)$"
  ) |>
  mutate(
    sim_f = as.factor(paste0("lambda == ", sim)),
    num_f = as.factor(paste0("b == ", num_comp))
  ) |>
  group_by(time, sim_f, num_f) |>
  summarize(
    ni_mean = mean(nint_absz),
    ni_se = sd(nint_absz) / sqrt(n()),
    it_mean = mean(ints_absz),
    it_se = sd(ints_absz) / sqrt(n()),
    .groups = "drop"
  ) |>
  mutate(
    ni_lower = ni_mean - 2 * ni_se,
    ni_upper = ni_mean + 2 * ni_se,
    it_lower = it_mean - 2 * it_se,
    it_upper = it_mean + 2 * it_se
  )


std_err_plot <- ggplot(data = sim_study_std_err) +
  geom_ribbon(
    aes(x = time, ymin = ni_lower, ymax = ni_upper),
    alpha = 0.2, fill = "#858585"
  ) + 
  geom_ribbon(
    aes(x = time, ymin = it_lower, ymax = it_upper),
    alpha = 0.2, fill = "#858585"
  ) +
  geom_line(aes(x = time, y = ni_mean)) +
  geom_line(aes(x = time, y = it_mean), linetype = "dashed") +
  facet_grid(vars(num_f), vars(sim_f), scales = "free_y", labeller = label_parsed) + 
  xlab("Post-Treatment Time") +
  ylab("Mean Standardized Error (Posterior Mean)") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(strip.background = element_rect(fill = "white", color = "black"))

ggsave(std_err_plot, device = "pdf", width = 5, height = 4, file = "../figs/ints_std_err.pdf", create.dir = TRUE)

## Overfitting plot

sim_study_overfit <- abs(sim_study_ints) |>
  pivot_longer(
    cols = !c(num_comp, sim),
    names_to = c("model", ".value"),
    names_transform = list(time = as.integer),
    names_pattern = "^([^_]+)_(.*)$"
  ) |>
  filter(num_comp == 3) |>
  mutate(
    errspur = (acor_3 + acor_4 + acor_5) / 3
  ) |>
  select(!starts_with("acor")) |>
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
      `No Int` = "nint",
      `With Int` = "ints"),
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
  facet_grid(vars(sim), vars(time), scales = "free", labeller = label_bquote(rows = lambda == .(sim))) +
  xlab("Population Absolute Correlation with Spurious") +
  ylab("Mean Standardized Error (Posterior Mean)") +
  scale_x_continuous(expand = expansion(mult = 0.6)) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(strip.background = element_rect(fill = "white", color = "black"))

ggsave(overfit_plot, device = "pdf", width = 6, height = 3.5, file = "../figs/ints_overfit.pdf", create.dir = TRUE)
