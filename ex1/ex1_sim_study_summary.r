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

sim_study_abs_err <- abs(sim_study_stat) |>
  select(contains("mean_")) |>
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "time"),
    names_transform = list(time = as.integer),
    names_pattern = "(.*)_(\\d+)$"
  ) |>
  group_by(time) |>
  summarize(
    ns_mean = mean(nonstat_mean),
    ns_se = sd(nonstat_mean) / sqrt(n()),
    st_mean = mean(stat2_mean),
    st_se = sd(stat2_mean) / sqrt(n())
  ) |>
  mutate(
    ns_lower = ns_mean - 2 * ns_se,
    ns_upper = ns_mean + 2 * ns_se,
    st_lower = st_mean - 2 * st_se,
    st_upper = st_mean + 2 * st_se
  )

abs_mad_plot <- ggplot(data = sim_study_abs_err) +
  geom_ribbon(
    aes(x = time, ymin = ns_lower, ymax = ns_upper),
    alpha = 0.1, color = "grey"
  ) + 
  geom_ribbon(
    aes(x = time, ymin = st_lower, ymax = st_upper),
    alpha = 0.1, color = "grey"
  ) +
  geom_line(aes(x = time, y = ns_mean)) +
  geom_line(aes(x = time, y = st_mean), linetype = "dashed") +
  xlab("Post-Treatment Time") +
  ylab("Mean Absolute Error (Posterior Mean)") +
  theme_bw()

ggsave(abs_mad_plot, device = "pdf", width = 5, height = 4, file = "stat_abs_err.pdf")

sim_study_std_err <- abs(sim_study_stat) |>
  select(contains("absz_")) |>
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "time"),
    names_transform = list(time = as.integer),
    names_pattern = "(.*)_(\\d+)$"
  ) |>
  group_by(time) |>
  summarize(
    ns_mean = mean(nonstat_absz),
    ns_se = sd(nonstat_absz) / sqrt(n()),
    st_mean = mean(stat2_absz),
    st_se = sd(stat2_absz) / sqrt(n())
  ) |>
  mutate(
    ns_lower = ns_mean - 2 * ns_se,
    ns_upper = ns_mean + 2 * ns_se,
    st_lower = st_mean - 2 * st_se,
    st_upper = st_mean + 2 * st_se
  )


std_err_plot <- ggplot(data = sim_study_std_err) +
  geom_ribbon(
    aes(x = time, ymin = ns_lower, ymax = ns_upper),
    alpha = 0.1, color = "grey"
  ) + 
  geom_ribbon(
    aes(x = time, ymin = st_lower, ymax = st_upper),
    alpha = 0.1, color = "grey"
  ) +
  geom_line(aes(x = time, y = ns_mean)) +
  geom_line(aes(x = time, y = st_mean), linetype = "dashed") +
  xlab("Post-Treatment Time") +
  ylab("Mean Standardized Error (Posterior Mean)") +
  theme_bw()

ggsave(std_err_plot, device = "pdf", width = 5, height = 4, file = "stat_std_err.pdf")

sim_study_overfit <- abs(sim_study_stat) |>
  pivot_longer(
    cols = everything(),
    names_to = c("model", ".value"),
    names_transform = list(time = as.integer),
    names_pattern = "^([^_]+)_(.*)$"
  ) |>
  rename(
    perc = pred_perc,
    mad = pred_mad
  ) |>
  pivot_longer(
    cols = contains("_"),
    names_to = c(".value", "time"),
    names_transform = list(time = as.integer),
    names_pattern = "(.*)_(\\d+)$"
  ) |>
  group_by(time, model) |>
  summarize(
    mean_absz = mean(absz),
    mean_mad = mean(mad)
  ) |>
  mutate(
    model = fct_recode(as.factor(model),
      "Nonstat" = "nonstat",
      "Weaker" = "stat2",
      "Stronger" = "stat1"),
    time = paste0("Time ", time)
  )

overfit_plot <- ggplot(data = sim_study_overfit) +
  # geom_point(aes(x = mean_mad, y = mean_absz, group = time)) +
  geom_line(aes(x = mean_mad, y = mean_absz, group = time), linewidth = 0.8) +
  geom_label(aes(label = model, x = mean_mad, y = mean_absz), size = 3) +
  facet_wrap(vars(time), ncol = 1, strip.position = "right") +
  xlab("Average Error (MAD)") +
  ylab("Mean Standardized Error (Posterior Mean)") +
  scale_y_continuous(limits = c(0.43, 0.8)) +
  scale_x_continuous(limits = c(1.25, 1.85)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(strip.background = element_rect(fill = "white", color = "black")) +
  ggtitle("(B)")

ggsave(overfit_plot, device = "pdf", width = 4, height = 5, file = "stat_overfit.pdf")
