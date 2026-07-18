library(dplyr)
library(forcats)
library(ggplot2)
library(tidyr)

source("../plotting.r")
load("sim_study_ints.RData")

ints_summary <- sim_study_ints |>
  group_by(num_comp, sim) |>
  summarize(
    absz_diff1 = mean(ints_absz_1 - nint_absz_1),
    absz_diff_se1 = sd(ints_absz_1 - nint_absz_1) / sqrt(n()),
    nint_absz_mean1 = mean(nint_absz_1),
    ints_absz_mean1 = mean(ints_absz_1)
  )

print(ints_summary)

perc_summary <- sim_study_ints |>
  group_by(num_comp, sim) |>
  summarize(
    mean_perc_nint = mean(nint_pred_perc),
    mean_perc_ints = mean(ints_pred_perc)
  )

print(perc_summary)

cors_summary <- sim_study_ints |>
  group_by(num_comp, sim) |>
  summarize(
    mean_cor_nint = mean(nint_acor_4),
    mean_cor_ints = mean(ints_acor_4),
    mean_cerr_nint = mean(nint_acor_err_4),
    mean_cerr_ints = mean(ints_acor_err_4)
  )

print(cors_summary)

## Statistic S2 predictive p-value (see paper Section 5), averaged separately for
## each of the four data generating conditions (num_comp x sim).
loc_cor_summary <- sim_study_ints |>
  group_by(num_comp, sim) |>
  summarize(
    mean_loc_cor_nint = mean(nint_loc_cor_pval),
    mean_loc_cor_ints = mean(ints_loc_cor_pval)
  )

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
    it_se = sd(ints_absz) / sqrt(n())
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
    alpha = 0.1, fill = "#858585"
  ) + 
  geom_ribbon(
    aes(x = time, ymin = it_lower, ymax = it_upper),
    alpha = 0.1, fill = "#858585"
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
  select(-ends_with("loc_cor_pval")) |>
  pivot_longer(
    cols = !c(num_comp, sim),
    names_to = c("model", ".value"),
    names_transform = list(time = as.integer),
    names_pattern = "^([^_]+)_(.*)$"
  ) |>
  rename(
    perc = pred_perc,
    mad = pred_mad,
    ncomp = num_comp
  ) |>
  filter(ncomp == 3) |>
  mutate(
    # simf = as.factor(paste0("lambda == ", sim)),
    errspur = (acor_3 + acor_4 + acor_5) / 3
  ) |>
  select(!starts_with("acor")) |>
  pivot_longer(
    cols = contains("_"),
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
    mean_err = mean(errspur)
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
