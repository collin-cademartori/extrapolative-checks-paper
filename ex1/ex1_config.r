## Shared model configuration for example 1.
##
## Sourced by ex1_sim_study.r (and its workers), ex1_pred_checks.r, ex1_rsq_priors.r and
## ex1_derive_scales.r, so the study, the figures and the derivation cannot disagree about what the
## models are. They had, three separate ways: the nonstationary arm's rho prior drifted to Beta(7,3)
## while the DGP kept Beta(8,2); the DGP ran at alpha_diag = 0 while every fit used 20; and
## ex1_derive_scales.r checked its constants against ETA_FRAC = 0.1 after the study had moved to
## 0.2. Each was invisible until measured.
##
## WHAT LIVES HERE: dimensionless multiples and the model's shape -- quantities that mean the same
## thing in every context.
##
## WHAT DOES NOT: anything carrying the units of a particular dataset. The study multiplies these by
## per-dataset rms_y / sd_y; the figures apply them in units where the anchor is 1; the derivation
## applies them to its own forward simulator. Sampler settings (iterations, chains, the escalation
## ladder) are the study's alone and stay there.

## ---- data-generating process ------------------------------------------------------------------
## The truth the study measures its arms against: an integrated latent structure with sigma = 1 and
## an iid level-error sd of 2. The nonstationary arm shares DGP_RHO and ALPHA_DIAG deliberately --
## being the correctly specified reference arm means matching the DGP in kind, not merely in
## functional form.
DGP_SIGMA <- 1
DGP_ETA <- 2
DGP_RHO <- c(7,3)

## ---- model shape ------------------------------------------------------------------------------
N_UNITS <- 8
T_TIMES <- 20
K_LATENT <- 4
NUM_TREATED <- 5

## Zero-avoiding inverse-gamma on the loading diagonal, repelling the collapsed-loading minor modes.
## An identification and sampling device, not a modelling claim -- a prior predictive check cannot
## see it, so it is not justified by one.
ALPHA_DIAG <- 20

## ---- autocorrelation priors ---------------------------------------------------------------------
RHO_NONSTAT <- DGP_RHO        # the reference arm is handed the DGP's own persistence prior
RHO_STAT <- c(98, 2)          # near-unit-root: mean 0.98

## ---- scale multiples ----------------------------------------------------------------------------
## Every scale is a fraction of an anchor measured on the data's own scale, so the arms are
## comparable to each other and to the truth. See ex1_derive_scales.r for where each comes from and
## for the guard that fails loudly if the model config moves out from under them.
##
##   SIGMA_MULT_NONSTAT  a units conversion. In the nonstationary branch sigma scales the
##                       DIFFERENCED signal, so a level anchor must be converted; integrating a
##                       T = 20 window costs about an order of magnitude. It is the self-consistency
##                       fixed point, and because this arm shares the DGP's rho prior the same value
##                       also hands it the DGP's own sigma and eta.
##
##   SIGMA_MULT_STAT     NOT a conversion. Under this misspecification no multiple matches the
##                       data's RMS and its SD at once, so it chooses: it matches the SD. A
##                       near-unit-root AR(1) realises only ~40% of its long-run SD over T = 20, so
##                       setting sigma to the expected SD implies series far flatter than intended.
##                       Matching the SD overshoots the RMS; matching the RMS would PROHIBIT the
##                       dispersion the analyst expects, and excluding a moment is a stronger
##                       assumption than exceeding one.
SIGMA_MULT_NONSTAT <- 1 / 7
SIGMA_MULT_STAT <- 2

## Error scales, as fractions of the anchor. eta is on the LEVEL scale in both branches --
## ife_named.stan applies the sqrt(2) differencing inflation itself -- so these take no differencing
## correction, and a caller-side one would double-count.
##
##   ETA_FRAC_NONSTAT  the DGP's own level-error sd, 2 x this arm's sigma.
##   ETA_FRAC_STAT     PROVISIONAL, still being settled. It is not a scale conversion but a
##                     prior-predictive-check choice about how much a linear trend may be
##                     attenuated: a large error flattens any trend, so only a small one lets a
##                     stationary model put enough prior mass on high |cor(y, t)| to be plausible.
##                     ex1_rsq_priors.r is where it is justified.
##
##                     It is COUPLED to SIGMA_MULT_STAT: sd(y)^2 ~ sd(signal)^2 + eta^2, so a larger
##                     error contributes realised dispersion and less sigma is needed to match the
##                     expected SD. Measured, the SD-matching multiple runs 1.98 / 1.95 / 1.85 /
##                     1.67 as ETA_FRAC_STAT goes 0.05 / 0.1 / 0.2 / 0.3 -- nearly flat up to 0.2,
##                     so the pair is well identified there, but a move to 0.3 would put
##                     SIGMA_MULT_STAT 17% off and trip the guard in ex1_derive_scales.r. The two
##                     constants are justified jointly, never one at a time.
ETA_FRAC_NONSTAT <- 2 * SIGMA_MULT_NONSTAT
ETA_FRAC_STAT <- 0.1

## Prior scale for the treatment effect, SHARED by every arm so they differ in their model of the
## DATA, not in their prior over the ESTIMAND. Each arm used to inherit its own sigma[1], which is
## not a common scale -- sigma is a differenced innovation scale in the nonstationary branch and a
## long-run marginal SD in the stationary one. That gave the nonstationary arm 1.00 against the
## stationary arms' 14.05 (19.9x once the cumsum construction is accounted for), a gap on the very
## quantity the study compares them on; and since the true effect here is exactly zero, the tighter
## prior shrank the reference arm toward the truth for free.
##
## Anchored on mean sd(y_n), not RMS: an effect is a CHANGE in the series, so the yardstick is how
## much a typical unit moves, not how far it sits from zero.
DELTA_FRAC <- 1.0

## Realised sd(y) per unit of sigma, for the stationary configuration over a T_TIMES window. A
## MEASURED property of the model, not a choice: a near-unit-root AR(1) realises far less dispersion
## over a short window than its long-run marginal SD, because the sample mean absorbs the
## low-frequency wandering. It is what makes "set sigma to the SD you expect" wrong by a factor of
## three, and it is the constant ex1_sd_priors.r turns into a figure.
##
## Derived and guarded by ex1_derive_scales.r, so it cannot silently go stale if T_TIMES, K_LATENT,
## ALPHA_DIAG, RHO_STAT or ETA_FRAC_STAT move -- it depends on all of them.
SD_PER_SIGMA <- 0.322
