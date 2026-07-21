plot_data_units <- function(data, data_comp = NULL, unit, samples = 16, hide_y = TRUE) {

  ys <- data$ys

  num_samples <- nrow(ys)
  plot_units <- sample.int(num_samples, size = samples)
  ys_plot <- ys[plot_units, , unit]

  if(!is.null(data_comp)) {
    ys_plot <- rbind(data_comp[, unit], ys_plot)
  }

  ys_long <- as.data.frame.table(ys_plot)
  names(ys_long) <-  c("sample", "time", "obs")
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

  if (hide_y) {
    plot <- plot +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }

  return(plot)
}

plot_post_fits_all <- function(data, post_ns, post_2, post_1) {
  ys_long <- as.data.frame.table(data)
  names(ys_long) <- c("time", "unit", "obs")
  ys_long$time <- as.integer(ys_long$time)

  means_ns <- as.data.frame.table(post_ns)
  names(means_ns) <- c("time", "unit", "obs")
  means_ns$time <- as.integer(means_ns$time)

  means_2 <- as.data.frame.table(post_2)
  names(means_2) <- c("time", "unit", "obs")
  means_2$time <- as.integer(means_2$time)

  means_1 <- as.data.frame.table(post_1)
  names(means_1) <- c("time", "unit", "obs")
  means_1$time <- as.integer(means_1$time)

  plot <- ggplot() +
    geom_line(
      data = ys_long,
      aes(x = time, y = obs), color = "grey"
    ) + 
    geom_line(
      data = means_ns,
      aes(x = time, y = obs), color = "black"
    ) +
    geom_line(
      data = means_2,
      aes(x = time, y = obs), color = "blue",
    ) +
    geom_line(
      data = means_1,
      aes(x = time, y = obs), color = "blue", linetype = "dashed"
    ) +
    facet_wrap(vars(unit), ncol = 2, scales = "free_y") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(strip.background = element_rect(fill = "white", color = "black")) +
    xlab("Time") +
    ylab("Outcome")

  return(plot)
}

# Plot each unit in a separate frame with posterior predictive samples
plot_data_matrix_post <- function(ys, post_ys) {
  ys_long <- as.data.frame.table(ys)
  names(ys_long) <-  c("time", "unit", "obs")
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
      linewidth = 0.1, color="black"
    ) +
    facet_wrap(vars(unit), ncol = 4, scales = "free_y") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(strip.background = element_rect(fill = "white", color = "black")) +
    xlab("Time") +
    ylab("Outcome")

  return(plot)
}

# Plot all time series in one frame with 
# highlighting for those most similar to the treated unit
plot_data_highlight <- function(data, cor_perc = 0.95, use_exp = TRUE, num_samples = 9) {
  trans <- ifelse(use_exp, exp, function(x) x)

  yis <- sample.int(dim(data$ys)[1], size = num_samples)
  ys <- data$ys[yis, , ]

  # cor_level <- quantile(cor(ys)[1,], cor_perc)
  cor_level <- 0.8

  ys_long <- as.data.frame.table(ys)
  names(ys_long) <-  c("sample", "time", "unit", "obs")

  ys_long$time <- as.integer(ys_long$time)
  #levels(ys_long$unit) <- paste0("Unit ", levels(ys_long$unit))
  # levels(ys_long$sample) <- paste0("Sample ", levels(ys_long$sample))

  y1 <- ys_long |> filter(unit == "A") |>
        select(time, sample, obs1 = obs)

  ys_long <- ys_long |>
    left_join(y1, by = c("time", "sample")) |>
    group_by(unit, sample) |>
    mutate(
      cor_y1 = cor(obs, obs1),
    ) |>
    ungroup()

  cor_cuts <- ys_long |> group_by(sample) |> 
    summarize(cor_cut = quantile(cor_y1, cor_perc)) |>
    mutate(sample_name = paste0("Sample ", sample, " (", round(cor_cut, 2), ")"))

  ys_long <- ys_long |>
    left_join(cor_cuts, by = "sample") |>
    mutate(
      is_comp = ifelse(cor_y1 >= cor_cut, "comp", "notcomp")
    ) |>
    # ungroup() |>
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
    facet_wrap(vars(sample), nrow = 2, scales = "free_y") +
    # scale_y_continuous(limits = c(
    #   ifelse(use_exp, 0, -6),
    #   ifelse(use_exp, 160, 6)
    # )) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(strip.background = element_rect(fill = "white", color = "black")) +
    xlab("Time") +
    ylab("Outcome")

  return(plot)
}