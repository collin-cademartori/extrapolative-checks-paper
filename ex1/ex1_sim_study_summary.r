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
# Optional CLI argument selects the results file, so a fast-mode run can be summarized without
# touching the full study's output: Rscript ex1_sim_study_summary.r sim_study_ns_fast.RData
.res_file <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(.res_file)) .res_file <- "sim_study_ns.RData"
cat(sprintf("\nLoading %s\n", .res_file))
load(.res_file)

# The checkpointing driver adds four bookkeeping columns (rep, unit, failed, error) that are not
# per-arm statistics. They must come off before anything below: `error` is character, so the
# abs(sim_study_stat) calls used to build the error and overfit frames fail outright on it, and all
# four would otherwise be swept into the pivots that split names on "<model>_<statistic>".
# Report failed reps first -- a study that lost reps to an error should not be summarized silently
# as though it were complete.
if ("failed" %in% names(sim_study_stat)) {
  n_failed <- sum(sim_study_stat$failed, na.rm = TRUE)
  cat(sprintf("\nReps loaded: %d", nrow(sim_study_stat)))
  if (n_failed > 0) {
    cat(sprintf("  --  %d FAILED, excluded from every summary below:\n", n_failed))
    for (m in unique(sim_study_stat$error[which(sim_study_stat$failed)])) cat("    - ", m, "\n")
    sim_study_stat <- sim_study_stat[!sim_study_stat$failed, , drop = FALSE]
  } else {
    cat("  (no failures)\n")
  }
}
sim_study_stat <- sim_study_stat |> select(-any_of(c("rep", "unit", "failed", "error")))

# Numeric-only results (the standardized-error curves and overfit trade-off are
# shown in the plots below): 95% posterior-predictive interval coverage -- the
# nominal level is met or exceeded by every model, so the check cannot rule out
# the misspecified model
# (paper Section 5) -- and the S1 time-correlation predictive p-value, averaged
# over the study for each model.

perc_summary <- sim_study_stat |>
  summarize(
    q5_perc_nonstat = quantile(nonstat_pred_perc, 0.05),
    q5_perc_stat_weak = quantile(stat_weak_pred_perc, 0.05),
    q5_perc_stat_strong = quantile(stat_strong_pred_perc, 0.05),
    q95_perc_nonstat = quantile(nonstat_pred_perc, 0.95),
    q95_perc_stat_weak = quantile(stat_weak_pred_perc, 0.95),
    q95_perc_stat_strong = quantile(stat_strong_pred_perc, 0.95),
    mean_perc_nonstat = mean(nonstat_pred_perc),
    mean_perc_stat_weak = mean(stat_weak_pred_perc),
    mean_perc_stat_strong = mean(stat_strong_pred_perc),
  )

cat("\nPer-model 95% interval coverage:\n")
print(as_tibble(perc_summary), width = Inf)

pval_summary <- sim_study_stat |>
  summarize(
    q5_pval_nonstat = quantile(nonstat_time_cor_pval, 0.05),
    q5_pval_stat_weak = quantile(stat_weak_time_cor_pval, 0.05),
    q5_pval_stat_strong = quantile(stat_strong_time_cor_pval, 0.05),
    q95_pval_nonstat = quantile(nonstat_time_cor_pval, 0.95),
    q95_pval_stat_weak = quantile(stat_weak_time_cor_pval, 0.95),
    q95_pval_stat_strong = quantile(stat_strong_time_cor_pval, 0.95),
    mean_pval_nonstat = mean(nonstat_time_cor_pval),
    mean_pval_stat_weak = mean(stat_weak_time_cor_pval),
    mean_pval_stat_strong = mean(stat_strong_time_cor_pval),
  )

cat("\nS1 time-correlation predictive p-value:\n")
print(as_tibble(pval_summary), width = Inf)

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
    # Overfitting of the treated unit's pre-treatment window (the extrapolation basis the
    # counterfactual is built from): the fraction of that unit's noise absorbed into its fitted
    # signal. Replaces pred_mad on the overfit plot -- pred_mad is mean|fitted - observed| over ALL
    # units, so it averages the treated unit together with seven that do not feed the delta estimate,
    # and it is scale-dependent. See the note in ex1_sim_study.r for the comparison behind this.
    mean_noise_abs = mean(noise_abs_tr),
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
  geom_line(aes(x = mean_noise_abs, y = mean_absz, group = time), linewidth = 0.8) +
  geom_label(aes(label = model, x = mean_noise_abs, y = mean_absz), size = 3) +
  facet_wrap(vars(time), ncol = 1, strip.position = "right") +
  xlab("Fraction of Treated Unit's Pre-Treatment\n Noise Absorbed by the Fitted Signal") +
  ylab("Average Standardized Error of Posterior Expected Treatment Effect") +
  scale_x_continuous(expand = expansion(mult = 0.3)) +
  scale_y_continuous(expand = expansion(mult = 0.2)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(strip.background = element_rect(fill = "white", color = "black")) +
  ggtitle("(B)")

ggsave(overfit_plot, device = "pdf", width = 2.5, height = 5, file = "../figs/stat_overfit.pdf", create.dir = TRUE)
