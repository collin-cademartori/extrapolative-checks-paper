plot_data_units <- function(data, data_comp = NULL, unit, samples = 16, hide_y = TRUE,
                            n_x_breaks = NULL) {
  ys <- data$ys

  num_samples <- nrow(ys)
  plot_units <- sample.int(num_samples, size = samples)
  ys_plot <- ys[plot_units, , unit]

  if (!is.null(data_comp)) {
    ys_plot <- rbind(data_comp[, unit], ys_plot)
  }

  ys_long <- as.data.frame.table(ys_plot)
  names(ys_long) <- c("sample", "time", "obs")
  ys_long$time <- as.integer(ys_long$time)
  ys_long$sample <- paste0("Sample ", ys_long$sample)

  plot <- ggplot(data = ys_long) +
    geom_line(aes(x = time, y = obs)) +
    facet_wrap(vars(sample), ncol = 7, scales = "free_y") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(strip.background = element_rect(fill = "white", color = "black")) +
    xlab("Time") +
    ylab("Outcome")

  # Thin the x-axis ticks (e.g. for long horizons where the default breaks collide).
  if (!is.null(n_x_breaks)) {
    plot <- plot + scale_x_continuous(n.breaks = n_x_breaks)
  }

  if (hide_y) {
    plot <- plot +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }

  return(plot)
}

# Observed series (grey) and the three posterior-mean fits (nonstationary black;
# stationary weak/strong blue solid/dashed) for a single unit. data, post_ns, post_2,
# post_1 are all time x unit matrices; `unit` selects the column to show.
# pred_rep: optionally, a SINGLE posterior predictive replicate for `unit`, drawn on top in red.
# Pass NULL (the default) to switch it off.
#
# Why it is worth seeing: the other three lines are posterior means of Y_latent, the noiseless
# latent signal. Statistic S1 is instead computed on Y_pred = that signal plus observation noise of
# scale tau * sigma. The two can look very different, so a model whose fitted mean tracks a trend
# closely can still fail S1 because the trend in its REPLICATES is diluted by noise. Only the
# replicate shows what S1 actually sees. One is drawn rather than several deliberately -- the
# pattern is meant to be read across datasets, and an envelope of draws would clutter the panel.
plot_post_fits_stat <- function(data, post_ns, post_2, post_1, unit, pred_rep = NULL) {
  series <- function(m) data.frame(time = seq_len(nrow(m)), obs = m[, unit])

  plot <- ggplot(mapping = aes(x = time, y = obs)) +
    geom_line(data = series(data), color = "grey") +
    geom_line(data = series(post_ns), color = "black") +
    geom_line(data = series(post_2), color = "blue") +
    geom_line(data = series(post_1), color = "blue", linetype = "dashed") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    xlab("Time") +
    ylab("Outcome")

  # Added last so it draws on top of the means.
  if (!is.null(pred_rep)) {
    plot <- plot + geom_line(
      data = data.frame(time = seq_along(pred_rep), obs = as.numeric(pred_rep)),
      color = "red", alpha = 0.8
    )
  }

  return(plot)
}

# Plot each unit in a separate frame with posterior predictive samples
plot_data_matrix_post <- function(ys, post_ys) {
  ys_long <- as.data.frame.table(ys)
  names(ys_long) <- c("time", "unit", "obs")
  ys_long$time <- as.integer(ys_long$time)
  ys_long$sample <- 0

  post_ys <- as.data.frame.table(post_ys)
  names(post_ys) <- c("sample", "time", "unit", "obs")
  post_ys$time <- as.integer(post_ys$time)
  post_ys$sample <- as.integer(post_ys$sample)

  post_ys <- rbind(post_ys, ys_long)
  post_ys$unit <- paste0("Unit ", post_ys$unit)
  ys_long$unit <- paste0("Unit ", ys_long$unit)

  post_ys <- post_ys |>
    group_by(unit, time) |>
    summarize(
      y_min = quantile(obs, 0.005),
      y_max = quantile(obs, 0.995)
    )

  plot <- ggplot() +
    geom_ribbon(
      data = post_ys,
      aes(ymin = y_min, ymax = y_max, x = time),
      fill = "grey", alpha = 0.5
    ) +
    geom_line(
      data = ys_long,
      aes(x = time, y = obs),
      linewidth = 0.1, color = "black"
    ) +
    facet_wrap(vars(unit), ncol = 4, scales = "free_y") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(strip.background = element_rect(fill = "white", color = "black")) +
    xlab("Time") +
    ylab("Outcome")

  return(plot)
}

# Comparator-status plot for the intercepts simulation, for a single model fit. Each
# simulated dataset is one line per unit. 
plot_intercepts_fits <- function(test_ys, cor_sq, groups, num_treated) {
  T_times <- nrow(test_ys)
  T_pre <- T_times - num_treated
  # Color switches (and the marker sits) at the last pre-treatment time, so the
  # gradient spans only fully-pre-treatment segments; the segment crossing into
  # treatment takes the post (true-group) color.
  boundary <- T_pre

  # Blue gradient over a fixed [0, 1] correlation domain, so shade intensity is
  # comparable across datasets and models. The top of the ramp matches the "true"
  # group color, so a correctly-identified comparator holds one blue across the
  # boundary while a spurious one flips blue -> vermilion (the reveal).
  blue_ramp <- colorRamp(c("#d6e2ff", "#0072b2"))
  grad_color <- function(x) {
    m <- blue_ramp(pmin(pmax(x, 0), 1))
    rgb(m[, 1], m[, 2], m[, 3], maxColorValue = 255)
  }
  # Okabe-Ito-based colorblind-safe groups: treated black, true blue (matches the
  # gradient top), spurious vermilion (warm reveal against the blue), uncorrelated grey.
  group_color <- c(
    treated = "#000000", true = "#0072b2",
    spurious = "#d55e00", uncorrelated = "#999999"
  )

  df <- as.data.frame.table(test_ys)
  names(df) <- c("time", "unit", "obs")
  df$time <- as.integer(df$time)
  df$unit <- as.integer(df$unit)
  df$period <- ifelse(df$time <= boundary, "pre", "post")

  # Default color is the unit's true-group color (used post-treatment and for the
  # treated unit throughout); untreated pre-treatment rows overwrite it with the
  # model-inferred correlation gradient.
  df$color <- unname(group_color[groups[df$unit]])
  pre_un <- df$period == "pre" & groups[df$unit] != "treated"
  df$color[pre_un] <- grad_color(cor_sq[df$unit[pre_un] - 1])

  # Repeat the boundary time in the post segment (true-group color) so each line
  # connects across the switch rather than breaking at it.
  bridge <- df[df$time == boundary, ]
  bridge$period <- "post"
  bridge$color <- unname(group_color[groups[bridge$unit]])
  seg <- rbind(df, bridge)

  plot <- ggplot(seg) +
    geom_line(aes(
      x = time, y = obs,
      group = interaction(unit, period), color = color
    )) +
    geom_vline(xintercept = boundary, color = "grey40", linetype = "solid") +
    scale_color_identity() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    xlab("Time") +
    ylab("Outcome")

  return(plot)
}

# Plot all series in one frame: grey for the rest, blue for the units most
# correlated with the treated unit, and green for the treated unit itself.
plot_data_highlight <- function(data, cor_perc = 0.95, use_exp = TRUE, num_samples = 9) {
  trans <- ifelse(use_exp, exp, function(x) x)

  yis <- sample.int(dim(data$ys)[1], size = num_samples)
  ys <- data$ys[yis, , ]

  ys_long <- as.data.frame.table(ys)
  names(ys_long) <- c("sample", "time", "unit", "obs")

  ys_long$time <- as.integer(ys_long$time)

  y1 <- ys_long |>
    filter(unit == "A") |>
    select(time, sample, obs1 = obs)

  ys_long <- ys_long |>
    left_join(y1, by = c("time", "sample")) |>
    group_by(unit, sample) |>
    mutate(
      cor_y1 = cor(obs, obs1),
    ) |>
    ungroup()

  cor_cuts <- ys_long |>
    group_by(sample) |>
    summarize(cor_cut = quantile(cor_y1, cor_perc)) |>
    mutate(sample_name = paste0("Sample ", sample, " (", round(cor_cut, 2), ")"))

  ys_long <- ys_long |>
    left_join(cor_cuts, by = "sample") |>
    mutate(
      is_comp = case_when(
        unit == "A" ~ "treated", # the treated unit itself
        cor_y1 >= cor_cut ~ "comp", # top cor_perc% most correlated with treated
        TRUE ~ "notcomp"
      )
    ) |>
    mutate(obs = trans(obs)) |>
    select(-obs1) |>
    select(-cor_cut)

  levels(ys_long$sample) <- as.character(cor_cuts$sample_name)

  plot <- ggplot() +
    geom_line(data = filter(ys_long, is_comp == "notcomp"), aes(
      x = time, y = obs,
      group = unit
    ), color = "grey", alpha = 0.3) +
    geom_line(data = filter(ys_long, is_comp == "comp"), aes(
      x = time, y = obs,
      group = unit
    ), color = "#3a3aff", alpha = 1) +
    geom_line(data = filter(ys_long, is_comp == "treated"), aes(
      x = time, y = obs,
      group = unit
    ), color = "#009e73", linewidth = 0.8, alpha = 1) +
    facet_wrap(vars(sample), nrow = 2, scales = "free_y") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(strip.background = element_rect(fill = "white", color = "black")) +
    xlab("Time") +
    ylab("Outcome")

  return(plot)
}
