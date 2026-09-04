## Data-generating process for the intercepts example, shared by ex2_sim_study.r and
## ex2_derive_scales.r so the study and the derivation of its constants cannot describe
## different data. Requires ex2_config.r to be sourced first.

ruv <- function(d) {
  v <- rnorm(d)
  uv <- v / sqrt(sum(v * v))
  return(uv)
}

# Untreated units split into three groups differing in their pre-treatment correlation with the
# treated unit and in their long-run average, so location and correlation are entangled in a way
# the unit-intercepts model wrongly assumes independent.
#
# Returns Y (T x N, the fit orientation), the ground-truth group of each column, and the
# observation-noise sd used.
sim_model_intercepts <- function(
    N_unc = 2, N_comp_true = DGP_N_COMP_TRUE, N_comp_spur = 2,
    T_times = DGP_T_TIMES, T_treated = DGP_T_TREATED,
    K_unc = DGP_K_UNC, sim = DGP_SIM, level_offset) {
  N_units <- 1 + N_comp_true + N_comp_spur + N_unc
  K_gen <- 2 + K_unc

  f_treat <- level_offset + arima.sim(model = list(ar = 0.9), n = T_times)

  # f_alt matches the treated factor pre-treatment, then diverges downward over the treatment
  # window -- the driver of the "spurious" comparators.
  f_alt <- (f_treat - level_offset) +
    c(rep(0, T_times - T_treated), rep(-DGP_F_TREAT_SD, T_treated))

  # Reject until the "uncorrelated" factors are genuinely uncorrelated with the treated factor.
  f_unc <- matrix(nrow = K_unc, ncol = T_times)
  cor_unc <- Inf
  while (cor_unc > 0.01) {
    for (k in seq_len(K_unc)) {
      f_unc[k, ] <- rnorm(1, 0, 2) + arima.sim(model = list(ar = 0.9), n = T_times)
    }
    cor_unc <- max(abs(cor(t(f_unc), t(t(f_treat)))))
  }

  facs <- rbind(f_treat, f_alt, f_unc)
  loads <- matrix(nrow = N_units, ncol = K_gen)
  loads[1, ] <- c(1, rep(0, K_gen - 1))
  # True comparators load on the treated factor: genuine correlation throughout.
  for (n in seq_len(N_comp_true)) {
    loads[1 + n, ] <- c(sqrt(sim), 0, sqrt(1 - sim) * ruv(K_gen - 2))
  }
  # Spurious comparators load on f_alt: pre-treatment correlation only.
  for (n in seq_len(N_comp_spur)) {
    loads[1 + N_comp_true + n, ] <- c(0, sqrt(sim), sqrt(1 - sim) * ruv(K_gen - 2))
  }
  # Uncorrelated units carry no shared-factor signal.
  for (n in seq_len(N_unc)) {
    loads[1 + N_comp_true + N_comp_spur + n, ] <- c(0, 0, ruv(K_gen - 2))
  }

  # Treated and true comparators sit high (via f_treat), spurious low (via f_alt), so correlation
  # with the treated unit does not determine location -- contrary to what the intercepts model
  # assumes.
  lat <- loads %*% facs
  noise_sd <- DGP_NOISE_FRAC * mean(apply(lat, 1, sd))
  Y <- t(lat + rnorm(nrow(lat) * ncol(lat), sd = noise_sd))

  groups <- c(
    "treated",
    rep("true", N_comp_true),
    rep("spurious", N_comp_spur),
    rep("uncorrelated", N_unc)
  )

  return(list(Y = Y, groups = groups, noise_sd = noise_sd))
}
