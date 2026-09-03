## Shared model configuration for example 1.
##
## Sourced by other files in /ex1 so that common model assumptions are consistent throughout.
##
## Some model hyperparameters are expressed here as dimensionless ratios. In inference, these
## are multiplied by empirical point estimates of the denominators.

## ---- data-generating process ----------------------------------------------------------------------
## Ground truth for simulation studies: nonstationary factors with sigma = 1 and an iid error sd = 2.
## In simulation study, nonstationary model inherits DGP_RHO and ALPHA_DIAG, matching the DGP.
DGP_SIGMA <- 1
DGP_ETA <- 2
DGP_RHO <- c(7, 3)

## ---- model shape ----------------------------------------------------------------------------------
N_UNITS <- 8
T_TIMES <- 20
K_LATENT <- 4
NUM_TREATED <- 5

## Shape parameter for zero-avoiding inverse-gamma on the loading diagonal.
## Reduces probability of minor modes at zero and expresses belief that all factors should be used
## to explain data.
ALPHA_DIAG <- 20

## ---- autocorrelation priors -----------------------------------------------------------------------
RHO_NONSTAT <- DGP_RHO        # nonstationary model uses DGP's factor autocorrelation prior
RHO_STAT <- c(98, 2)          # stationary model gets large autocorrelations to mimic nonstationarity

## ---- sigma multiples ------------------------------------------------------------------------------
## In simulation studies, sigma is a multiple of the data root mean square.
## Parameter sigma has different meanings in the two models, so multiples are used to refer
## sigma to a common baseline and ensure scaling is comparable across models.
## This ensures that comparisons between models reflect differences between the factor process
## assumptions, not incidental differences in the scaling of their latent components.
## Script ex1_derive_scales.r derives these constants from scaling invariants and errors
## if the constants recorded here no longer satisfy the desired invariants.
##
##   SIGMA_MULT_NONSTAT  In the nonstationary branch sigma scales the differenced
##                       signal; integrating a T = 20 window increases the outcomes scale
##                       by about an order of magnitude.
##                       Coupled to DGP_RHO -- a DGP with lower rho integrates to a
##                       smaller RMS, so the multiple rises to compensate.
##
##   SIGMA_MULT_STAT     When rho is high in the stationary model, the sample SD of the
##                       latent factors over T = 20 points is substantially less than the
##                       long-run SD, so sigma is inflated to compensate.
##                       This tends to overestimate the RMS of the data, but not to the
##                       point of ruling it out.
##                       Coupled to ETA_FRAC_STAT, which scales the error term and also
##                       affects the outcome SD.
SIGMA_MULT_NONSTAT <- 1 / 6.2
SIGMA_MULT_STAT <- 2

## Error scales, as fractions of average RMS. Eta is expressed on the outcome scale in both
## models, not on the differenced scale in the nonstationary case like sigma.
##
##   ETA_FRAC_NONSTAT  the DGP's own iid error sd, 2 x this arm's sigma.
##
##   ETA_FRAC_STAT     Justified by prior predictive check of the correlation between time
##                     and outcome. Large error levels attenuate this correlation, so must
##                     be small enough for stationary model to mimic nonstationarity in the
##                     short run and pass predictive check.
##                     Coupled to SIGMA_MULT_STAT since the combination of these two scales
##                     determines the overall outcome scale.
ETA_FRAC_NONSTAT <- 2 * SIGMA_MULT_NONSTAT
ETA_FRAC_STAT <- 0.1

## Spread of the truncated-normal prior on eta, as a coefficient of variation, shared by
## stationary and nonstationary models.
ETA_CV_EX1 <- 1.0

## Prior scale for the treatment effect, multiplied by the average data SD.
## Expresses idea that treatment effects are typically smaller than the scale of the data.
## Shared by stationary and nonstationary models, ensuring identical treatment effect priors
## so that differences in posterior inference are not confounded with differences in the prior.
DELTA_FRAC <- 0.5

## Realised outcome SD over sigma, for the stationary configuration over T = 20 window.
## Outcome SD is substantially lower than long-run SD in stationary model when
## autocorrelation is high due to large dependence on starting point.
## This is a measured property of the model which is verified in the derivation script
## ex1_derive_scales.r.
## Coupled to T_TIMES, K_LATENT, ALPHA_DIAG, RHO_STAT and ETA_FRAC_STAT.
SD_PER_SIGMA <- 0.322
