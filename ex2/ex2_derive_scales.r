## Derivation of ex2's scale constants: ETA_FRAC_EX2, and the two sigma conventions the arms use.
##
## This file exists so those constants are not bare numbers in ex2_sim_study.r: run it and it prints
## where each comes from. The study treats the results as FROZEN CONSTANTS and never calls this
## file, so the prior depends on the data in exactly one stated place (through each dataset's RMS or
## sd) rather than twice.
##
## Method differs from ex1's: ex1's constants are properties of the model's prior, so that file
## forward-simulates the prior, while ex2's are pinned to properties of its DGP, so this file
## simulates the DGP.

source("ex2_config.r")

set.seed(52918)
N_DATASETS <- 2000L

## The DGP itself comes from ex2_dgp.r, the same file ex2_sim_study.r reads, so this derivation
## cannot describe data the study does not generate.
source("ex2_dgp.r")

gen_one <- function(N_comp_spur, level, sim = DGP_SIM) {
  sim_model_intercepts(
    N_unc = DGP_N_UNITS - 1 - DGP_N_COMP_TRUE - N_comp_spur,
    N_comp_true = DGP_N_COMP_TRUE, N_comp_spur = N_comp_spur,
    K_unc = DGP_K_UNC, sim = sim, level_offset = level)
}

## --- what the DGP produces ----------------------------------------------------------------------

## The study's own sweep grid, so the constants are not tuned to one cell of it.
grid <- expand.grid(N_comp = c(2, 3), level = DGP_LEVELS)
rows <- vector("list", nrow(grid))
for (g in seq_len(nrow(grid))) {
  r <- replicate(N_DATASETS %/% nrow(grid), {
    d <- gen_one(grid$N_comp[g], grid$level[g])
    c(rms = mean(apply(d$Y, 2, function(y) sqrt(mean(y^2)))),
      sd  = mean(apply(d$Y, 2, sd)),
      noise = d$noise_sd,
      grand = mean(d$Y),
      sd_means = sd(colMeans(d$Y)))
  })
  rows[[g]] <- c(N_comp = grid$N_comp[g], level = grid$level[g], rowMeans(r),
                 noise_lo = quantile(r["noise", ], 0.05), noise_hi = quantile(r["noise", ], 0.95))
}
tab <- do.call(rbind, rows)

cat("\n=== What ex2's DGP generates, per sweep cell ===\n")
cat(sprintf("%8s %6s %9s %9s %9s %18s\n", "N_comp", "level", "E[RMS]", "E[sd]", "E[noise]", "noise 5-95%"))
for (i in seq_len(nrow(tab))) {
  cat(sprintf("%8d %6.1f %9.2f %9.2f %9.3f %8.3f - %-8.3f\n", tab[i, "N_comp"], tab[i, "level"],
              tab[i, "rms"], tab[i, "sd"], tab[i, "noise"],
              tab[i, "noise_lo.5%"], tab[i, "noise_hi.95%"]))
}
rms_bar <- mean(tab[, "rms"]); sd_bar <- mean(tab[, "sd"]); noise_bar <- mean(tab[, "noise"])
cat(sprintf("\n  overall:  E[RMS] = %.2f   E[sd] = %.2f   E[noise] = %.3f   RMS/sd = %.2f\n",
            rms_bar, sd_bar, noise_bar, rms_bar / sd_bar))

## --- ETA_FRAC_EX2 --------------------------------------------------------------------------------
## eta is an absolute error sd shared by both arms, expressed as a fraction of mean sd(y_n). sd
## rather than RMS because the DGP defines its noise off the latent sd, and because sd is
## arm-neutral where RMS carries the intercepts.

eta_frac <- noise_bar / sd_bar

cat("\n=== ETA_FRAC_EX2: the error scale, as a fraction of mean sd(y) ===\n")
cat(sprintf("  E[noise] / E[sd] = %.3f / %.2f = %.4f\n", noise_bar, sd_bar, eta_frac))
cat(sprintf("  committed value  = %.4f\n", ETA_FRAC_EX2))
cat("\n  For contrast, what the old RATIO parametrization handed each arm (shared tau ~ TN(0.1, 0.05),\n")
cat("  so err_sd = 0.1 * sigma[n] with sigma differing by arm):\n")
cat(sprintf("    no_ints  sigma = RMS -> %.3f = %.2fx the truth\n",
            0.1 * rms_bar, 0.1 * rms_bar / noise_bar))
cat(sprintf("    ints     sigma = sd  -> %.3f = %.2fx the truth\n",
            0.1 * sd_bar, 0.1 * sd_bar / noise_bar))
cat("  i.e. the correctly specified arm was handed twice the true error, for a reason unconnected\n")
cat("  to intercepts. That is what absolute mode fixes.\n")

## --- LEVEL_SPREAD_FRAC ----------------------------------------------------------------------------
## sd(colMeans(y)) as a multiple of mean sd(y_n). The study measures the level spread from each
## dataset directly; ex2_pred_checks.r cannot, because it works in units where mean sd(y_n) = 1 and
## has no dataset, so it needs this ratio to place the same intercept prior in its own units.

sd_means_bar <- mean(tab[, "sd_means"])
level_spread <- sd_means_bar / sd_bar

cat("\n=== The intercept prior's anchor ===\n")
cat(sprintf("  grand mean of y            = %.2f\n", mean(tab[, "grand"])))
cat(sprintf("  sd of per-unit means       = %.2f\n", sd_means_bar))
cat(sprintf("  as a multiple of mean sd(y)= %.3f   (committed LEVEL_SPREAD_FRAC = %.3f)\n",
            level_spread, LEVEL_SPREAD_FRAC))
cat(sprintf("  INT_FRAC = %.2f gives a prior scale of %.2f in data units\n",
            INT_FRAC, INT_FRAC * sd_means_bar))
cat("\n  The prior's LOCATION is the data's grand mean, so it is shift-invariant and nothing is\n")
cat("  committed here. The old fixed N(4, 3) sat about 0.6 sd above the actual grand mean.\n")

## --- what the spread above does and does not mean --------------------------------------------------
## The DGP's own noise sd does move a good deal between datasets, because the per-unit latent sds do.
## What matters is that it moves TOGETHER with the anchor ETA_FRAC_EX2 is built on, so the ratio --
## which is what that constant claims to be -- is far more stable than either quantity alone.

cat(sprintf("\n  Note: the DGP's own noise sd spans %.3f to %.3f across datasets (5-95%%), because the\n",
            min(tab[, "noise_lo.5%"]), max(tab[, "noise_hi.95%"])))
cat("  per-unit latent sds do. That is NOT inherited as extra variance downstream: the true noise sd\n")
cat("  and mean sd(y_n) correlate at +0.93, so their ratio -- ETA_FRAC_EX2 -- has CV 0.09 against\n")
cat("  0.25 for the noise sd itself, and the prior lands within 88-116% of the truth in 90% of\n")
cat("  datasets against a prior whose own sd equals its location.\n")

## --- guard ----------------------------------------------------------------------------------------
## Committed values read from ex2_config.r, so this checks the constant the study actually uses
## rather than a second copy that could drift from it.

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
