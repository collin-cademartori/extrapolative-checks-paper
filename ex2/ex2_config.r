## Shared model configuration for example 2.
##
## Sourced by other files in /ex2 so that common model assumptions are consistent throughout.
##
## Some model hyperparameters are expressed here as dimensionless ratios. In inference, these
## are multiplied by empirical point estimates of the denominators.

## ---- data-generating process ---------------------------------------------------------------
## Untreated units split into groups differing in their pre-treatment correlation with the treated
## unit AND in their long-run average, so location and correlation are entangled in exactly the way
## the unit-intercepts model assumes they are not.
DGP_T_TIMES <- 30
DGP_T_TREATED <- 5
DGP_N_UNITS <- 8           # 1 treated + N_COMP_TRUE + N_comp_spur + N_unc, held fixed
DGP_N_COMP_TRUE <- 2
DGP_K_UNC <- 1

## Post-treatment divergence of the spurious factor f_alt, as a multiple of the realised sd of the
## treated factor in that dataset. Anchored rather than absolute so the divergence is the same size
## relative to the series it perturbs in every dataset, as every other scale here is.
DGP_F_TREAT_FRAC <- 1.0

## Swept axis 1: the gap between the treated and true group's long-run average and the spurious
## group's. This is the quantity that makes location informative about correlation, which is the
## entanglement the unit-intercepts model assumes away, so it is swept rather than fixed.
##
## The effect runs through no_ints, which has no gamma: a unit's level must come out of
## Lambda * Phi_means, the same loadings that generate correlation, so separating levels forces the
## spurious loadings away from the treated ones. ints puts the level in gamma and is unaffected.
DGP_LEVELS <- c(5, 10)

## Swept axis 2 is N_comp, the number of spurious comparators, set in ex2_sim_study.r's call.

## Strength of the spurious units' pre-treatment correlation. Fixed rather than swept: it has
## almost no effect on the contrast between the arms, unlike the level gap.
DGP_SIM <- 0.9

## Observation-noise sd, as a fraction of the average unit's latent sd. Anchoring it this way holds
## the signal-to-noise ratio roughly constant across datasets, which a fixed noise level would not:
## the true noise sd and mean sd(y_n) correlate strongly, so their ratio is far more stable than
## either quantity alone.
DGP_NOISE_FRAC <- 0.121

## ---- model shape -----------------------------------------------------------------------------
K_LATENT <- 3
NUM_TREATED <- DGP_T_TREATED

## Shape parameter for zero-avoiding inverse-gamma on the loading diagonal.
## Reduces probability of minor modes at zero and expresses belief that all factors should be used
## to explain data.
ALPHA_DIAG <- 20

## ---- autocorrelation prior --------------------------------------------------------------------
RHO_EX2 <- c(90, 10)       # mean 0.9

## ---- unit-intercept prior ---------------------------------------------------------------------
## The ints arm's gamma[n] ~ N(location, scale). The independence of a unit's level from its
## correlation with the treated unit is the assumption this example is about; what is set here is
## only where the prior sits, not that assumption.
##
## Anchored rather than fixed, so that the prior is shift-invariant and the same rule serves both
## the study and the prior predictive figure:
##     location = the data's grand mean          (0 in the figure's centred units)
##     scale    = INT_FRAC * sd(colMeans(y))
##
## INT_FRAC = 1: the prior scale on unit intercepts is the observed spread of unit means, with no
## calibration factor.
##
## CAVEAT: sd(colMeans(y)) is estimated from DGP_N_UNITS = 8 numbers, so its own sampling CV is
## ~27%. The prior's width will wobble between datasets for reasons unrelated to the level spread.
INT_FRAC <- 1.0

## sd(colMeans(y)) as a multiple of mean sd(y_n), averaged over the sweep. Needed only by
## ex2_pred_checks.r, which works in units where mean sd(y_n) = 1 and so cannot measure the level
## spread itself. Coupled to DGP_LEVELS, since the numerator scales with the level gap; derived and
## guarded by ex2_derive_scales.r.
LEVEL_SPREAD_FRAC <- 2.50

## ---- error scale ------------------------------------------------------------------------------
## eta is a single absolute observation-error sd on the data's own scale, shared by both arms. Their
## sigma differs by design -- no_ints produces each series' level from its factor means and anchors
## on RMS, ints has gamma for the level and anchors on sd -- and eta being absolute is what stops
## that justified difference leaking into the error scale, where nothing would justify it.
##
## Anchored on mean sd(y_n): the DGP defines its noise off the latent sd, and sd is arm-neutral
## where RMS carries the intercepts.
##
## ETA_CV_EX2 is the prior sd of eta over its prior location, and matches ex1's ETA_CV_EX1, so that
## neither study asserts the error scale more confidently than the other.
ETA_FRAC_EX2 <- 0.12
ETA_CV_EX2 <- 1.0

## Prior scale for the treatment effect, multiplied by the average data SD.
## Expresses idea that treatment effects are typically smaller than the scale of the data.
## Shared by both arms, ensuring identical treatment effect priors so that differences in posterior
## inference are not confounded with differences in the prior. Anchored on mean sd(y_n) rather than
## RMS because an effect is a change in the series, not a distance from zero.
DELTA_FRAC_EX2 <- 0.5
