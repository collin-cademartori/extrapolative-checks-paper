## Shared model configuration for example 2.
##
## Sourced by other files in /ex2 so that common model assumptions are consistent throughout.
##
## Some model hyperparameters are expressed here as dimensionless ratios. In inference, these
## are multiplied by empirical point estimates of the denominators.

## ---- data-generating process ---------------------------------------------------------------
## Untreated units split into groups differing in their latent correlation with the treated
## unit and in their long-run average.
DGP_T_TIMES <- 30
DGP_T_TREATED <- 5
DGP_N_UNITS <- 8
# Number of units with strong latent correlation with treated
DGP_N_COMP_TRUE <- 2
# Number of factors to be used to generate the units which have small latent 
# correlation with treated
DGP_K_UNC <- 1

## Post-treatment divergence of the factor defining the spuriously correlated units, as a multiple 
## of the sample sd of the factor defining the treated (and truly correlated) units.
DGP_F_TREAT_FRAC <- 1.0

## Level gap between factor defining spurious and factor defining treated.
## DGP_LEVELS records the range explored while establishing that the effect is not specific to one
## configuration. STUDY_LEVEL and STUDY_N_COMP are the single condition the study runs, and are
## what ex2_derive_scales.r calibrates against -- the between-unit spread scales with the level
## gap while the within-unit sd does not, so a constant averaged over a wider grid misstates it.
DGP_LEVELS <- c(5, 10)
STUDY_LEVEL <- 10
STUDY_N_COMP <- 3

## Strength of the spurious units' pre-treatment correlation.
DGP_SIM <- 0.9

## SD of iid noise, as a fraction of the average SD of units' latent components,
## holding the signal-to-noise ratio roughly constant across datasets.
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
## For the unit-level intercepts model, INT_FRAC defines the prior SD of the intercepts as a
## of the SD of the units' sample means, i.e. prior scale    = INT_FRAC * sd(colMeans(y)).
INT_FRAC <- 1.0

## SD of unit's sample means divided by average of units' sample SDs. Needed by
## ex2_pred_checks.r, which simulates units with long-run SDs equalt to 1. 
## Coupled to DGP_LEVELS, since the numerator scales with the level gap; value checked
## in derivation script ex2_derive_scales.r.
LEVEL_SPREAD_FRAC <- 3.35

## ---- error scale ------------------------------------------------------------------------------
## The iid noise SD in both models is assigned a truncated normal prior based on the mean of 
## units' sample SDs. The prior location parameter is ETA_FRAC_EX2 times this SD, and the prior
## scale parameter is ETA_CV_EX2 * ETA_FRAC_EX2 times this SD. 
ETA_FRAC_EX2 <- 0.128
ETA_CV_EX2 <- 1.0

## Prior scale for the treatment effect, multiplied by the treated unit's PRE-treatment SD.
## Expresses idea that treatment effects are typically smaller than the variation already present
## in the series they act on. Pre-treatment only, so that a large effect cannot widen its own prior.
## Both arms compute it from the data identically, so differences in posterior inference are not
## confounded with differences in the prior.
DELTA_FRAC_EX2 <- 0.5
