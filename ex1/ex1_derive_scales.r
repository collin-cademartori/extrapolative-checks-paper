## Derivation of the two sigma multiples used by ex1:
##
##     SIGMA_MULT_NONSTAT = 1/7     SIGMA_MULT_STAT = 2
##
## Both convert the anchor RMS(y_n) -- measured from each dataset -- into the `overall_scales`
## (sigma) that ife_named.stan receives. This file exists so those constants are not bare numbers in
## ex1_sim_study.r: run it and it prints where each comes from. The study treats the results as
## FROZEN CONSTANTS and never calls this file, so that the prior depends on the data in exactly one
## stated place (through RMS(y_n)) rather than twice.
##
## Not derived here: ETA_FRAC_WEAK = 0.1 and ETA_FRAC_STRONG = 0.05. Those are not scale
## conversions -- they are prior-predictive-check choices about how much a linear trend may be
## attenuated, and ex1_rsq_priors.r is where they are justified.
##
## Method: forward-simulate the prior of ife_named.stan in plain R. No Stan and no MCMC -- the
## prior predictive is a forward model, so it can be sampled directly. Everything here is a faithful
## transcription of that file's `parameters`, `transformed parameters` and `generated quantities`
## blocks; see transcription notes at each function.

set.seed(20260831)
N_DRAWS <- 20000L

## Every constant below comes from ex1_config.r, the single source the study and the figures read
## too. This file previously kept its own copies and they drifted: ETA_FRAC_STAT sat at 0.1 here
## after the study had moved to 0.2, so the guard at the foot of this file was checking the
## committed multiples against a configuration nobody was fitting.
source("ex1_config.r")

M_UNITS <- N_UNITS
# The DGP and every fitted arm share alpha_diag, so the loadings prior is one setting rather than
# two. (It was not always: run_sim_study_stat used to pass none and get the default half-normal
# diagonal while the arms got the zero-avoiding inverse-gamma. Both have E||Lambda[n,]||^2 = 1, so
# the gap was only a few percent, but it meant the "correctly specified" reference arm differed
# from the DGP in a way nothing recorded.)
DGP_ALPHA_DIAG <- ALPHA_DIAG
FIT_ALPHA_DIAG <- ALPHA_DIAG

# The nonstationary arm's error, as a multiple of its own sigma: the DGP's ratio, 2/1. Fixing it
# here is what makes that arm's multiple a one-line solve rather than a root-find. RHO_NONSTAT
# equals DGP_RHO in the config, deliberately -- with a common rho prior, "reproduce the RMS of the
# data" and "receive the DGP's own sigma" are the SAME condition, so the nonstationary multiple has
# a single unambiguous target instead of two competing ones.
ETA_OVER_SIGMA_NONSTAT <- DGP_ETA / DGP_SIGMA

## --- the model's prior, transcribed ------------------------------------------------------------

# Lambda: cholesky_factor_cov[M, K] -- lower trapezoidal, positive diagonal. Prior from the
# `Loadings prior` block of ife_named.stan. With alpha_diag > 2 the first K diagonal entries take a
# zero-avoiding inverse-gamma whose second moment is 1/n, matching the half-normal it replaces, so
# every row keeps E||Lambda[n,]||^2 = 1 either way.
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

# Phi: K independent AR(1) columns, each scaled to UNIT MARGINAL VARIANCE (ar_process() in
# ife_named.stan uses scale_inno = sqrt(1 - rho^2) and seeds the chain at innovations[1]).
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

# One prior-predictive dataset. Y_means = sigma[n] * (Phi Lambda')[, n]; observation error is iid
# N(0, eta) on the LEVEL in both branches. In the nonstationary branch ife_named.stan works on
# differences and reconstructs with cumulative_sum -- its differenced error covariance
# (2 * eta^2 * errors_cov) is exactly Cov of diff(e) for iid level e, so cumsum returns iid level
# errors and simulating them directly here is equivalent, not an approximation.
draw_dataset <- function(sigma, eta, a_rho, b_rho, nonstationary, alpha_diag,
                         M = M_UNITS, T_times = T_TIMES, K = K_LATENT) {
  signal <- draw_phi(T_times, K, a_rho, b_rho) %*% t(draw_lambda(M, K, alpha_diag))
  signal <- signal * rep(sigma, each = T_times)
  if (nonstationary) signal <- apply(signal, 2, cumsum)
  signal + matrix(rnorm(T_times * M, 0, eta), T_times, M)
}

# E[RMS(y_n)] and E[sd(y_n)] over units and draws -- the two shape statistics the multiples are
# solved against.
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

## --- what the data look like --------------------------------------------------------------------
## Both multiples are solved so the fitted arm's prior predictive reproduces a shape statistic of
## the datasets the study actually fits. So measure those first.

dgp <- prior_moments(DGP_SIGMA, DGP_ETA, DGP_RHO[1], DGP_RHO[2], nonstationary = TRUE,
                     alpha_diag = DGP_ALPHA_DIAG)
dgp_sd_over_rms <- unname(dgp["sd"] / dgp["rms"])

cat("\n=== The datasets ex1 fits (ex1_sim_study.r's DGP) ===\n")
cat(sprintf("  E[RMS(y)] = %.3f     E[sd(y)] = %.3f     sd/RMS = %.4f\n",
            dgp["rms"], dgp["sd"], dgp_sd_over_rms))
cat("  Only the RATIO matters below: every multiple is expressed per unit of RMS(y), so the\n")
cat("  derivation is invariant to the overall scale of the data.\n")

## --- SIGMA_MULT_NONSTAT -------------------------------------------------------------------------
## Work in units where the anchor RMS(y) = 1 and solve for c in sigma = c * RMS(y).
## For this arm eta = 2 * sigma, so BOTH scales are proportional to c and therefore so is RMS(y):
## simulate once at c = 1 and read the multiple straight off. No root-find, no target beyond
## self-consistency -- the arm reproduces the RMS of the data it is fitted to.

ns1 <- prior_moments(1, ETA_OVER_SIGMA_NONSTAT * 1, RHO_NONSTAT[1], RHO_NONSTAT[2],
                     nonstationary = TRUE, alpha_diag = FIT_ALPHA_DIAG)
mult_nonstat <- 1 / unname(ns1["rms"])

cat("\n=== SIGMA_MULT_NONSTAT: match the RMS ===\n")
cat(sprintf("  at sigma = 1 (eta = %.1f x sigma) the prior predictive gives E[RMS(y)] = %.3f\n",
            ETA_OVER_SIGMA_NONSTAT, ns1["rms"]))
cat(sprintf("  RMS is proportional to sigma here, so the self-consistent multiple is 1/%.2f = %.4f\n",
            ns1["rms"], mult_nonstat))
cat(sprintf("  committed value 1/7 = %.4f\n", 1 / 7))
cat("\n  Because this arm shares the DGP's rho prior, that one solve satisfies both criteria at\n")
cat("  once: the arm reproduces the RMS of the data it is fitted to, AND it receives the DGP's own\n")
cat(sprintf("  sigma and eta -- %.3f and %.3f against the true 1.000 and 2.000. There is no\n",
            mult_nonstat * unname(ns1["rms"]), ETA_OVER_SIGMA_NONSTAT * mult_nonstat * unname(ns1["rms"])))
cat("  trade-off to adjudicate here, unlike the stationary multiple below.\n")
cat("\n  Why the multiple is so far below 1: sigma scales the DIFFERENCED signal in this branch\n")
cat("  (Y_means = sigma * Lambda_Phi, then cumulative_sum), while the anchor is measured on the\n")
cat("  level. Integrating a T = 20 window inflates the scale by roughly an order of magnitude.\n")

## --- SIGMA_MULT_STAT ----------------------------------------------------------------------------
## The stationary model is MISSPECIFIED for these data, so no multiple reproduces the RMS and the SD
## at once and the derivation has to choose. Solve for both and compare.
##
## Here eta is a fixed fraction of the anchor, not of sigma, so RMS(y) is not exactly proportional
## to c and the solve is a root-find rather than a division.

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

## Why the two fixed points are so far apart -- a closed form, for interpretation only; nothing
## above depends on it. A near-unit-root AR(1) realises far less dispersion over a short window
## than its long-run marginal sd, because the sample mean absorbs the low-frequency wandering.
ar1_sd_shrinkage <- function(rho, T_times) {
  k <- 1:(T_times - 1)
  sqrt(((T_times - 1) - (2 / T_times) * sum((T_times - k) * rho^k)) / (T_times - 1))
}
cat(sprintf("\n  (E[sample sd]/sigma_longrun for one AR(1) at T = %d: %.3f at rho = 0.97,\n",
            T_TIMES, ar1_sd_shrinkage(0.97, T_TIMES)))
cat(sprintf("   %.3f averaged over rho ~ Beta(%g, %g). The RMS is unaffected by this shrinkage,\n",
            mean(sapply(rbeta(20000, RHO_STAT[1], RHO_STAT[2]), ar1_sd_shrinkage, T_times = T_TIMES)),
            RHO_STAT[1], RHO_STAT[2]))
cat("   which is why the two fixed points separate.)\n")

## --- guard --------------------------------------------------------------------------------------
## The constants depend on T_times, K_latent, alpha_diag and both rho priors. Fail loudly if a
## change to any of those has moved them, rather than letting ex1_sim_study.r drift silently.

## --- SD_PER_SIGMA -------------------------------------------------------------------------------
## The realised sd of a T_TIMES window per unit of sigma, at the stationary configuration stated
## WITHOUT an anchor: eta = (ETA_FRAC_STAT / SIGMA_MULT_STAT) * sigma. ex1_sd_priors.r turns this
## one number into its figure -- sigma set to the expected SD implies series only this fraction as
## dispersed as intended -- so it needs the same protection as the multiples above.

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
# Committed values read from ex1_config.r, so this guard checks the constants the study and the
# figures actually use rather than a second copy that could drift from them.
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
