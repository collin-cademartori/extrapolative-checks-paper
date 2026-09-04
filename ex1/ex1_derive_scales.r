## Derivation of the two sigma multiples used by ex1:
##
##     SIGMA_MULT_NONSTAT = 1/6.2     SIGMA_MULT_STAT = 2
##
## Both mulitply the observed RMS of each dataset to give the value of sigma used by each model.
## This file exists to ensure these constants satisfy the desired constraints
## (see ex1_config.r for details).
##
## The ETA_FRAC_STAT is not derived, since this quantity is justified directly by a
## prior predictive check of the correlation between outcome and time in the stationary model.
##
## Derivation of constants is performed by simulating from the prior predictive distribution of
## ife_named.stan.

set.seed(20260831)
N_DRAWS <- 20000L

## Constants read from config, which is also used in predictive checks and simulation study,
## ensuring derived constraints apply the models run in those files.
source("ex1_config.r")

M_UNITS <- N_UNITS
# Inverse gamma shape parameter for diagonal loadings prior shared across DGP and both models
DGP_ALPHA_DIAG <- ALPHA_DIAG
FIT_ALPHA_DIAG <- ALPHA_DIAG

# The nonstationary arm's error, as a multiple of its own sigma: the DGP's ratio, 2/1. Fixing it
# here is what makes that arm's multiple a one-line solve rather than a root-find. RHO_NONSTAT
# equals DGP_RHO in the config, deliberately -- with a common rho prior, "reproduce the RMS of the
# data" and "receive the DGP's own sigma" are the SAME condition, so the nonstationary multiple has
# a single unambiguous target instead of two competing ones.
ETA_OVER_SIGMA_NONSTAT <- DGP_ETA / DGP_SIGMA

## --- the prior simulation -----------------------------------------------------------------------

# Lambda: Lower triangular with positive diagonal. Simulation here matches `Loadings prior` block
# of ife_named.stan. 
draw_lambda <- function(M, K, alpha_diag) {
  L <- matrix(0, M, K)
  for (n in seq_len(M)) {
    if (alpha_diag > 2 && n <= K) {
      if (n > 1) L[n, 1:(n - 1)] <- rnorm(n - 1, 0, 1 / sqrt(n))
      # Stan's inv_gamma(a, b) has density prop. to x^(-a-1) exp(-b/x), i.e. 1/X ~ Gamma(a, rate=b).
      L[n, n] <- 1 / rgamma(1, shape = alpha_diag,
                            rate = sqrt((alpha_diag - 1) * (alpha_diag - 2) / n))
    } else {
      m <- min(K, n)
      L[n, 1:m] <- rnorm(m, 0, 1 / sqrt(m))
    }
  }
  L
}

# Phi: K independent AR(1) columns, each scaled to unit marginal long-run variance. Matches
# definition of factors using AR process function ar_process() in ife_named.stan.
draw_phi <- function(T_times, K, a_rho, b_rho) {
  Phi <- matrix(0, T_times, K)
  for (k in seq_len(K)) {
    rho <- rbeta(1, a_rho, b_rho)
    e <- rnorm(T_times)
    x <- numeric(T_times); x[1] <- e[1]
    s <- sqrt(1 - rho^2)
    for (t in 2:T_times) x[t] <- rho * x[t - 1] + s * e[t]
    Phi[, k] <- x
  }
  Phi
}

# Draw dataset from prior predictive. Matches implementation logic of ife_named.stan
draw_dataset <- function(sigma, eta, a_rho, b_rho, nonstationary, alpha_diag,
                         M = M_UNITS, T_times = T_TIMES, K = K_LATENT) {
  signal <- draw_phi(T_times, K, a_rho, b_rho) %*% t(draw_lambda(M, K, alpha_diag))
  signal <- signal * rep(sigma, each = T_times)
  if (nonstationary) signal <- apply(signal, 2, cumsum)
  signal + matrix(rnorm(T_times * M, 0, eta), T_times, M)
}

# Estimate E[RMS(y_n)] and E[sd(y_n)] for each unit from prior predictive samples.
prior_moments <- function(sigma, eta, a_rho, b_rho, nonstationary, alpha_diag,
                          n_draws = N_DRAWS) {
  rms <- sds <- numeric(n_draws)
  for (i in seq_len(n_draws)) {
    Y <- draw_dataset(sigma, eta, a_rho, b_rho, nonstationary, alpha_diag)
    rms[i] <- mean(apply(Y, 2, function(y) sqrt(mean(y^2))))
    sds[i] <- mean(apply(Y, 2, sd))
  }
  c(rms = mean(rms), sd = mean(sds))
}

## --- prior predictive moments --------------------------------------------------------------------
## Both constant multiples are chosen so that the corresponding model's prior predictive 
## reproduces an RMS or SD of the data. These are estimated from samples here.

dgp <- prior_moments(DGP_SIGMA, DGP_ETA, DGP_RHO[1], DGP_RHO[2], nonstationary = TRUE,
                     alpha_diag = DGP_ALPHA_DIAG)
dgp_sd_over_rms <- unname(dgp["sd"] / dgp["rms"])

cat("\n=== The datasets ex1 fits (ex1_sim_study.r's DGP) ===\n")
cat(sprintf("  E[RMS(y)] = %.3f     E[sd(y)] = %.3f     sd/RMS = %.4f\n",
            dgp["rms"], dgp["sd"], dgp_sd_over_rms))
cat("  Only the RATIO matters below: every multiple is expressed per unit of RMS(y), so the\n")
cat("  derivation is invariant to the overall scale of the data.\n")

## --- SIGMA_MULT_NONSTAT -------------------------------------------------------------------------
## Simulation study defines sigma = c * RMS(y) for each dataset y. Here we take RMS = 1 and solve
## for c. In this model we also have eta = 2 * sigma, so both scales are proportional to c.

ns1 <- prior_moments(1, ETA_OVER_SIGMA_NONSTAT * 1, RHO_NONSTAT[1], RHO_NONSTAT[2],
                     nonstationary = TRUE, alpha_diag = FIT_ALPHA_DIAG)
mult_nonstat <- 1 / unname(ns1["rms"])

cat("\n=== SIGMA_MULT_NONSTAT: match the RMS ===\n")
cat(sprintf("  at sigma = 1 (eta = %.1f x sigma) the prior predictive gives E[RMS(y)] = %.3f\n",
            ETA_OVER_SIGMA_NONSTAT, ns1["rms"]))
cat(sprintf("  RMS is proportional to sigma here, so the self-consistent multiple is 1/%.2f = %.4f\n",
            ns1["rms"], mult_nonstat))
cat(sprintf("  committed value = %.4f\n", SIGMA_MULT_NONSTAT))
cat("\n  Because this arm shares the DGP's rho prior, that one solve satisfies both criteria at\n")
cat("  once: the arm reproduces the RMS of the data it is fitted to, AND it receives the DGP's own\n")
cat(sprintf("  sigma and eta -- %.3f and %.3f against the true 1.000 and 2.000. There is no\n",
            mult_nonstat * unname(ns1["rms"]), ETA_OVER_SIGMA_NONSTAT * mult_nonstat * unname(ns1["rms"])))
cat("  trade-off to adjudicate here, unlike the stationary multiple below.\n")
cat("\n  Why the multiple is so far below 1: sigma scales the DIFFERENCED signal in this branch\n")
cat("  (Y_means = sigma * Lambda_Phi, then cumulative_sum), while the anchor is measured on the\n")
cat("  level. Integrating a T = 20 window inflates the scale by roughly an order of magnitude.\n")

## --- SIGMA_MULT_STAT ----------------------------------------------------------------------------
## Because stationary model is misspecified, no multiple reproduces the RMS and the SD
## at once.
## Here, eta is a fixed fraction of the RMS, so the solution is found by root finding rather 
## than division.

stat_moments <- function(c) prior_moments(c, ETA_FRAC_STAT, RHO_STAT[1], RHO_STAT[2],
                                          nonstationary = FALSE, alpha_diag = FIT_ALPHA_DIAG,
                                          n_draws = N_DRAWS %/% 4L)
mult_stat_rms <- uniroot(function(c) stat_moments(c)["rms"] - 1, c(0.5, 5), tol = 1e-3)$root
mult_stat_sd <- uniroot(function(c) stat_moments(c)["sd"] - dgp_sd_over_rms, c(0.5, 5), tol = 1e-3)$root

at_rms <- stat_moments(mult_stat_rms)
at_sd <- stat_moments(mult_stat_sd)

cat("\n=== SIGMA_MULT_STAT: match the SD, not the RMS ===\n")
cat(sprintf("  multiple that reproduces the data's RMS : %.2f\n", mult_stat_rms))
cat(sprintf("  multiple that reproduces the data's SD  : %.2f     <- committed value 2\n", mult_stat_sd))
cat("\n  They cannot both be met -- the model is misspecified -- so compare what each gives up:\n")
cat(sprintf("    at c = %.2f :  E[RMS] = %.3f (%3.0f%% of the data's)   E[sd] = %.3f (%3.0f%%)\n",
            mult_stat_rms, at_rms["rms"], 100 * at_rms["rms"] / 1,
            at_rms["sd"], 100 * at_rms["sd"] / dgp_sd_over_rms))
cat(sprintf("    at c = %.2f :  E[RMS] = %.3f (%3.0f%%)               E[sd] = %.3f (%3.0f%%)\n",
            mult_stat_sd, at_sd["rms"], 100 * at_sd["rms"] / 1,
            at_sd["sd"], 100 * at_sd["sd"] / dgp_sd_over_rms))
cat("\n  Matching the RMS would cap the prior predictive SD near 60% of the data's -- the prior\n")
cat("  would PROHIBIT the dispersion the data actually show. Matching the SD instead overshoots\n")
cat("  the RMS by about half. Excluding a realised moment is a strong assumption made silently;\n")
cat("  exceeding one is a mild and visible assumption, so the SD fixed point is the conservative\n")
cat("  choice. The SD is also the quantity that governs how far the fitted factors can track the\n")
cat("  error term, which is the overfitting this arm exists to exhibit.\n")

## --- Guard --------------------------------------------------------------------------------------

## --- SD_PER_SIGMA -------------------------------------------------------------------------------
## The sample SD of a T_TIMES window divded by sigma.

eta_over_sigma <- ETA_FRAC_STAT / SIGMA_MULT_STAT
sd_per_sigma <- unname(prior_moments(1, eta_over_sigma, RHO_STAT[1], RHO_STAT[2],
                                     nonstationary = FALSE, alpha_diag = FIT_ALPHA_DIAG)["sd"])

cat("\n=== SD_PER_SIGMA: how much dispersion a short window actually realises ===\n")
cat(sprintf("  config, anchor-free: eta = %.2f x sigma\n", eta_over_sigma))
cat(sprintf("  E[realised sd(y)] / sigma = %.3f   over T = %d at rho ~ Beta(%g, %g)\n",
            sd_per_sigma, T_TIMES, RHO_STAT[1], RHO_STAT[2]))
cat(sprintf("  -> to realise an expected SD of S, set sigma = %.2f x S\n", 1 / sd_per_sigma))
cat(sprintf("  committed value %.3f\n", SD_PER_SIGMA))

TOL <- 0.10
# Committed values read from ex1_config.r.
committed <- c(SIGMA_MULT_NONSTAT = SIGMA_MULT_NONSTAT, SIGMA_MULT_STAT = SIGMA_MULT_STAT,
               SD_PER_SIGMA = SD_PER_SIGMA)
derived <- c(SIGMA_MULT_NONSTAT = mult_nonstat, SIGMA_MULT_STAT = mult_stat_sd,
             SD_PER_SIGMA = sd_per_sigma)
gap <- derived / committed - 1

cat("\n=== Committed constants ===\n")
for (nm in names(committed)) {
  cat(sprintf("  %-19s committed %.4f   derived %.4f   %+.1f%%   %s\n", nm,
              committed[[nm]], derived[[nm]], 100 * gap[[nm]],
              if (abs(gap[[nm]]) < TOL) "ok" else "OUT OF TOLERANCE"))
}
if (any(abs(gap) >= TOL)) {
  stop(sprintf("derived multiples differ from ex1_sim_study.r's by more than %.0f%%: %s. ",
               100 * TOL, paste(names(gap)[abs(gap) >= TOL], collapse = ", ")),
       "Either the model config changed (T_times, K_latent, alpha_diag, a rho prior) or the ",
       "constants need updating -- do not widen this tolerance to make it pass.")
}
cat(sprintf("\nBoth within %.0f%% of the committed values.\n", 100 * TOL))
