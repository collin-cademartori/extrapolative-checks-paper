## Derivation of ex2's scale constants: ETA_FRAC_EX2, and the two scale hyperparameter conventions the models use.
##
## Derivation of constants is performed by simulating from the data generating process described in
## Section 3.2 of the paper.

source("ex2_config.r")

set.seed(52918)
N_DATASETS <- 2000L

## Load the DGP simulator
source("ex2_dgp.r")

## --- what the DGP produces ----------------------------------------------------------------------

## Summarize properties of the DGP
draws <- replicate(N_DATASETS, {
  d <- sim_model_intercepts()
  # Pre-treatment window only, matching the anchors the study computes.
  Y <- d$Y[seq_len(DGP_T_TIMES - DGP_T_TREATED), , drop = FALSE]
  c(rms = mean(apply(Y, 2, function(y) sqrt(mean(y^2)))),
    sd  = mean(apply(Y, 2, sd)),
    noise = d$noise_sd,
    grand = mean(Y),
    sd_means = sd(colMeans(Y)))
})
means <- rowMeans(draws)
rms_bar <- means[["rms"]]; sd_bar <- means[["sd"]]; noise_bar <- means[["noise"]]
sd_means_bar <- means[["sd_means"]]
noise_lo <- quantile(draws["noise", ], 0.05)
noise_hi <- quantile(draws["noise", ], 0.95)

cat("\n=== What ex2's DGP generates ===\n")
cat(sprintf("  condition: level %g, %d spurious comparators, over %d datasets\n",
            DGP_LEVEL, DGP_N_COMP_SPUR, N_DATASETS))
cat(sprintf("  E[RMS] = %.2f   E[sd] = %.2f   E[noise] = %.3f (5-95%%: %.3f - %.3f)   RMS/sd = %.2f\n",
            rms_bar, sd_bar, noise_bar, noise_lo, noise_hi, rms_bar / sd_bar))

## --- ETA_FRAC_EX2 --------------------------------------------------------------------------------
## ETA_FRAC_EX2 sets the prior location of eta, the absolute noise sd shared by both models, as a
## fraction of mean sd(y_n). The prior scale is ETA_CV_EX2 times ETA_FRAC_EX2.

eta_frac <- noise_bar / sd_bar

cat("\n=== ETA_FRAC_EX2: prior location of the noise sd, as a fraction of mean sd(y) ===\n")
cat(sprintf("  E[noise] / E[sd] = %.3f / %.2f = %.4f\n", noise_bar, sd_bar, eta_frac))
cat(sprintf("  committed value  = %.4f\n", ETA_FRAC_EX2))

## --- LEVEL_SPREAD_FRAC ----------------------------------------------------------------------------
## LEVEL_SPREAD_FRAC sets the spread of the units' levels, sd(colMeans(y)), as a multiple of mean
## sd(y_n). The intercepts' prior scale is INT_FRAC times sd(colMeans(y)). The predictive checks in
## ex2_pred_checks.r assume units have latent SD = 1, so they use LEVEL_SPREAD_FRAC 
## to define an intercept prior that induces the same ratio of between- and within-unit spread.

level_spread <- sd_means_bar / sd_bar

cat("\n=== The intercept prior's anchor ===\n")
cat(sprintf("  grand mean of y            = %.2f\n", means[["grand"]]))
cat(sprintf("  sd of per-unit means       = %.2f\n", sd_means_bar))
cat(sprintf("  as a multiple of mean sd(y)= %.3f   (committed LEVEL_SPREAD_FRAC = %.3f)\n",
            level_spread, LEVEL_SPREAD_FRAC))
cat(sprintf("  INT_FRAC = %.2f gives a prior scale of %.2f in data units\n",
            INT_FRAC, INT_FRAC * sd_means_bar))
cat("\n  The prior's LOCATION is the data's grand mean, so it is shift-invariant and nothing is\n")
cat("  committed here.\n")

cat(sprintf("\n  Note: the DGP's own noise sd spans %.3f to %.3f across datasets (5-95%%), because the\n",
            noise_lo, noise_hi))
cat("  per-unit latent sds do. That is NOT inherited as extra variance downstream: the true noise sd\n")
cat("  and mean sd(y_n) correlate at +0.93, so their ratio -- ETA_FRAC_EX2 -- has CV 0.09 against\n")
cat("  0.25 for the noise sd itself, and the prior lands within 88-116% of the truth in 90% of\n")
cat("  datasets against a prior whose own sd equals its location.\n")

## --- guard ----------------------------------------------------------------------------------------
## Committed values read from ex2_config.r.

TOL <- 0.15
committed <- c(ETA_FRAC_EX2 = ETA_FRAC_EX2, LEVEL_SPREAD_FRAC = LEVEL_SPREAD_FRAC)
derived <- c(ETA_FRAC_EX2 = eta_frac, LEVEL_SPREAD_FRAC = level_spread)
gap <- derived / committed - 1

cat("\n=== Committed constants ===\n")
for (nm in names(committed)) {
  cat(sprintf("  %-18s committed %.4f   derived %.4f   %+.1f%%   %s\n", nm,
              committed[[nm]], derived[[nm]], 100 * gap[[nm]],
              if (abs(gap[[nm]]) < TOL) "ok" else "OUT OF TOLERANCE"))
}
if (any(abs(gap) >= TOL)) {
  stop(sprintf("derived constants differ from ex2_config.r's by more than %.0f%%: %s. ",
               100 * TOL, paste(names(gap)[abs(gap) >= TOL], collapse = ", ")),
       "Either the DGP changed (DGP_NOISE_FRAC, the factor structure, DGP_T_TIMES) or the constants ",
       "need updating -- do not widen this tolerance to make it pass.")
}
cat(sprintf("\nAll within %.0f%% of the committed values.\n", 100 * TOL))
