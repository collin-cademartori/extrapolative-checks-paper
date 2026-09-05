library(dplyr)
library(forcats)
library(ggplot2)
library(tidyr)

source("../plotting.r")
# Optional CLI argument selects the results file, so a fast-mode run can be summarized without
# touching the full study's output: Rscript ex2_sim_study_summary.r sim_study_ints_fast.RData
.res_file <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(.res_file)) .res_file <- "sim_study_ints.RData"
cat(sprintf("\nLoading %s\n", .res_file))
load(.res_file)

# Drop the driver's bookkeeping columns (rep, failed, error): they are not per-arm statistics and
# would otherwise be swept into the pivots that split names on "<model>_<statistic>". `level` and
# `num_comp` stay, being the study's design factors. Failed tasks are reported before being
# excluded, so a study that lost tasks is not summarized as though it were whole.
if ("failed" %in% names(sim_study_ints)) {
  n_failed <- sum(sim_study_ints$failed, na.rm = TRUE)
  cat(sprintf("\nTasks loaded: %d", nrow(sim_study_ints)))
  if (n_failed > 0) {
    cat(sprintf("  --  %d FAILED, excluded from every summary below:\n", n_failed))
    for (m in unique(sim_study_ints$error[which(sim_study_ints$failed)])) cat("    - ", m, "\n")
    sim_study_ints <- sim_study_ints[!sim_study_ints$failed, , drop = FALSE]
  } else {
    cat("  (no failures)\n")
  }
}
sim_study_ints <- sim_study_ints |> select(-any_of(c("rep", "failed", "error")))

# Posterior-predictive interval coverage, averaged per condition. The standardized error and
# correlation summaries are conveyed by the plots below.
perc_summary <- sim_study_ints |>
  group_by(num_comp, level) |>
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
## each data generating condition (num_comp x level).
loc_cor_summary <- sim_study_ints |>
  group_by(num_comp, level) |>
  summarize(
    q5_loc_cor_no_ints = quantile(no_ints_loc_cor_pval, 0.05),
    q95_loc_cor_no_ints = quantile(no_ints_loc_cor_pval, 0.95),
    q5_loc_cor_ints = mean(ints_loc_cor_pval, 0.05),
    q95_loc_cor_ints = mean(ints_loc_cor_pval, 0.95),
    .groups = "drop"
  )

cat("\nS2 location-correlation predictive p-value by condition (no-int vs with-int):\n")
print(loc_cor_summary, width = Inf)

# The study sweeps a grid, but the figures below show a single condition. The numeric summaries
# above still cover every cell.
PLOT_LEVEL <- 5
PLOT_NUM_COMP <- 3

# Signed relative bias, mean_k / sd_k. The recorded absz_k is already an absolute value, so it
# cannot be used here: the true effect is 0, so the signed quantity is what carries the bias.
for (.a in c("no_ints", "ints")) for (.k in 1:5) {
  sim_study_ints[[paste0(.a, "_relbias_", .k)]] <-
    sim_study_ints[[paste0(.a, "_mean_", .k)]] / sim_study_ints[[paste0(.a, "_sd_", .k)]]
}

# Mean and +/-2 SE bands over post-treatment time for a per-time statistic, for the no-intercepts
# and with-intercepts models. stat is "mean" (absolute bias) or "relbias" (relative bias). Values
# are SIGNED: the true effect is 0, so the mean of the posterior mean is the bias.
summarize_error <- function(stat) {
  sim_study_ints |>
    filter(level == PLOT_LEVEL, num_comp == PLOT_NUM_COMP) |>
    pivot_longer(
      cols = contains(paste0(stat, "_")),
      names_to = c(".value", "time"),
      names_transform = list(time = as.integer),
      names_pattern = "(.*)_(\\d+)$"
    ) |>
    group_by(time) |>
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
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "#999999") +
    xlab("Post-Treatment Time") +
    ylab(y_label) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(strip.background = element_rect(fill = "white", color = "black"))
}

sim_study_std_err <- summarize_error("relbias")
std_err_plot <- plot_error_bands(
  sim_study_std_err, "Average Relative Bias of Posterior\n Expected Treatment Effect"
)
ggsave(std_err_plot, device = "pdf", width = 5, height = 4, file = "../figs/ints_std_err.pdf", create.dir = TRUE)

sim_study_abs_err <- summarize_error("mean")
abs_err_plot <- plot_error_bands(
  sim_study_abs_err, "Average Absolute Bias of Posterior\n Expected Treatment Effect"
)
ggsave(abs_err_plot, device = "pdf", width = 5, height = 4, file = "../figs/ints_abs_err.pdf", create.dir = TRUE)

## Overfitting plot

sim_study_overfit <- sim_study_ints |>
  filter(level == PLOT_LEVEL, num_comp == PLOT_NUM_COMP) |>
  pivot_longer(
    cols = !c(num_comp, level),
    names_to = c("model", ".value"),
    names_transform = list(time = as.integer),
    names_pattern = "^(no_ints|ints)_(.*)$"
  ) |>
  # Separation: how far the model tells true comparators from spurious ones. Units 1:2 are the
  # true comparators and 3:(2 + num_comp) the spurious ones, in the DGP's generating order.
  mutate(
    sep = (cor_sq_1 + cor_sq_2) / 2 -
      rowMeans(pick(num_range("cor_sq_", 3:(2 + PLOT_NUM_COMP))))
  ) |>
  # Drop the per-unit correlation stats: they end in a digit but index units, not time.
  select(!starts_with("cor_sq") & !starts_with("acor_err")) |>
  # Second pivot only over the time-indexed columns; scalar per-model stats are left untouched.
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
  group_by(time, model) |>
  summarize(
    mean_relbias = mean(relbias),
    mean_sep = mean(sep),
    .groups = "drop"
  )

overfit_plot <- ggplot(data = sim_study_overfit) +
  geom_line(aes(x = mean_sep, y = mean_relbias), linewidth = 0.8) +
  geom_label(aes(label = model, x = mean_sep, y = mean_relbias), size = 3) +
  facet_wrap(vars(time), nrow = 1, scales = "free") +
  xlab("Modeled Separation of True from Spurious Comparators") +
  ylab("Average Relative Bias of Posterior\n Expected Treatment Effect") +
  scale_x_continuous(expand = expansion(mult = 0.6), n.breaks = 4) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(strip.background = element_rect(fill = "white", color = "black"))

ggsave(overfit_plot, device = "pdf", width = 6, height = 3.5, file = "../figs/ints_overfit.pdf", create.dir = TRUE)
