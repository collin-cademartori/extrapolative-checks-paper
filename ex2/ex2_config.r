## Shared model configuration for example 2.
##
## Sourced by ex2_sim_study.r (and its workers), ex2_pred_checks.r and ex2_derive_scales.r, so the
## study, the figure and the derivation cannot disagree about what the models are. They did: the
## prior predictive figure drew intercepts from N(0, 5) while the study fitted N(4, 3), ran at
## alpha_diag = 0 against the study's 20, and fixed the error scale where the study estimated it.
## None of that was visible without putting the two files side by side. See ex1/ex1_config.r for
## the same pattern and the same history.
##
## WHAT LIVES HERE: dimensionless multiples and the model's shape -- quantities meaning the same
## thing in every context.
##
## WHAT DOES NOT: anything carrying a particular dataset's units. The study multiplies these by
## per-dataset RMS / sd; the figure applies them in its own units. Sampler settings (iterations,
## chains, the escalation ladder) and the study's sweep grid stay with the study.

## ---- data-generating process ---------------------------------------------------------------
## The adversarial DGP of paper Section 5: untreated units split into groups differing in their
## pre-treatment correlation with the treated unit AND in their long-run average, so location and
## correlation are entangled in exactly the way the unit-intercepts model assumes they are not.
DGP_T_TIMES <- 30
DGP_T_TREATED <- 5
DGP_N_UNITS <- 8           # 1 treated + N_COMP_TRUE + N_comp_spur + N_unc, held fixed
DGP_N_COMP_TRUE <- 2
DGP_K_UNC <- 1

## Post-treatment divergence of the "spurious" factor f_alt, in absolute units. NOTE this is a fixed
## amount rather than one tied to each dataset's realised sd -- see EX2_PLAN.md section 5, still an
## open question.
DGP_F_TREAT_SD <- 1.9

## Level of the treated factor: f_treat = DGP_LEVEL_OFFSET + AR(0.9), and f_alt = f_treat -
## DGP_LEVEL_OFFSET, so this IS the level gap between the treated/true group and the spurious group.
##
## It was a bare literal 6 appearing twice inside sim_model_intercepts, which understated its role:
## it is effectively a tuning parameter for this example's headline contrast. Measured by sweeping it
## (2 datasets per level, fast settings), cor_sq with the SPURIOUS units responds very differently in
## the two arms:
##
##     level gap        no_ints        ints
##        ~1             0.61          0.64 - 0.73
##        ~3             0.48 - 0.61    0.64 - 0.74
##        ~6             0.58 - 0.61    0.65
##       ~10             0.39          0.64
##     slope           -0.0221/unit   -0.0049/unit   (cor -0.77 vs -0.40)
##
## no_ints falls 4.5x faster, because it has no gamma: a unit's level must come out of
## Lambda * Phi_means, the same loadings that generate correlation, so separating levels forces
## Lambda_spurious away from Lambda_treated. ints puts the level in gamma and is largely immune.
## The same mechanism shows on the other side -- no_ints' cor_sq with the TRUE comparators RISES with
## the gap (0.845 -> 0.979) while ints stays flat near 0.82.
##
## So at a small offset the two arms nearly agree and the example has little to show; at a large one
## they separate sharply. 6 sits mid-range. That is a defensible choice but it is a CHOICE, and the
## write-up should say so rather than let a reader assume the separation is purely a model property.
DGP_LEVEL_OFFSET <- 6

## Observation-noise sd, as a fraction of the LARGEST unit's latent sd. `max` is the least stable
## anchor available: measured over 300 datasets the realised noise sd ranges 0.103 to 0.446 (mean
## 0.201), so the truth itself varies four-fold rep to rep and everything downstream inherits that
## as extra variance. Also open in EX2_PLAN.md section 5.
DGP_NOISE_FRAC <- 0.1

## ---- model shape -----------------------------------------------------------------------------
K_LATENT <- 3
NUM_TREATED <- DGP_T_TREATED

## Zero-avoiding inverse-gamma on the loading diagonal. An identification and sampling device, not a
## modelling claim -- a prior predictive check cannot see it, so it is not justified by one.
ALPHA_DIAG <- 20

## ---- autocorrelation prior --------------------------------------------------------------------
RHO_EX2 <- c(90, 10)       # mean 0.9

## ---- unit-intercept prior ---------------------------------------------------------------------
## The `ints` arm's gamma[n] ~ N(location, scale). This is the assumption the example is about: the
## intercepts model treats a unit's LEVEL as independent of its CORRELATION with the treated unit,
## which the DGP deliberately violates. What is set here is only where that prior sits, not the
## independence assumption itself.
##
## ANCHORED, not fixed. It used to be N(4, 3): the 4 aimed at where the largest latent factors live
## (~6), but the DGP's actual grand mean is 2.01, because the low-level spurious and uncorrelated
## units pull it down -- so the prior sat about 0.6 sd above the levels it was meant to cover. A
## fixed location is also not shift-invariant: add 100 to every outcome and the model is badly
## misspecified for no reason.
##
## The rule is one rule, used in BOTH the study and the prior predictive figure:
##     location = the data's grand mean          (0 in the figure's centred units)
##     scale    = INT_FRAC * sd(colMeans(y))
## Using one rule matters: a figure drawn under a different prior than the study fits would not be
## checking the fitted model, which is the trap ex2_pred_checks.r was already in.
##
## INT_FRAC = 0.9 reproduces the old scale of 3 (the DGP's sd of per-unit means is 3.41), so this
## changes where the prior sits, not how wide it is.
##
## CAVEAT: sd(colMeans(y)) is estimated from DGP_N_UNITS = 8 numbers, so its own sampling CV is
## ~27%. The prior's width will wobble between datasets for reasons unrelated to the level spread.
INT_FRAC <- 0.9

## sd(colMeans(y)) as a multiple of mean sd(y_n), measured on the DGP. Needed only by
## ex2_pred_checks.r, which works in units where mean sd(y_n) = 1 and so cannot measure the level
## spread itself. Derived and guarded by ex2_derive_scales.r.
LEVEL_SPREAD_FRAC <- 2.05

## ---- error scale ------------------------------------------------------------------------------
## eta is a single ABSOLUTE observation-error sd on the data's own scale, shared by both arms.
##
## Both arms' sigma differs BY DESIGN -- no_ints must produce each series' level from its factor
## means, so it anchors on RMS; ints has gamma for the level, so it anchors on sd. Under the old
## ratio parametrization (err_sd = tau[n] * sigma[n], shared tau prior) that justified difference
## leaked into the error scale, where nothing justifies it:
##
##     E[mean RMS(y)] = 3.93,  E[mean sd(y)] = 1.66   ->  a 2.37x gap
##     true DGP noise sd = 0.201
##       no_ints got 0.1 * 3.93 = 0.393  ->  1.95x the truth
##       ints    got 0.1 * 1.66 = 0.166  ->  0.83x the truth
##
## i.e. the CORRECTLY SPECIFIED arm was handed twice the true error, for a reason unconnected to
## intercepts. Absolute mode closes that: sigma keeps its per-arm difference, eta cannot inherit it.
##
## Anchored on mean sd(y_n): the DGP defines its noise off the latent SD, and sd is arm-neutral
## where RMS carries the intercepts. 0.201 / 1.66 = 0.121, hence 0.12. The CV matches what the old
## tau ~ TN(0.1, 0.05) had.
ETA_FRAC_EX2 <- 0.12
ETA_CV_EX2 <- 0.5

## Prior scale for the treatment effect, SHARED by both arms so they differ in their model of the
## DATA, not in their prior over the ESTIMAND. Without it each arm inherits its own sigma[1], and
## since those anchor on RMS and sd respectively the arms get effect priors 2.37x apart -- the same
## leak as the error scale, on the quantity the study actually compares them on. See ife_named.stan's
## delta_scale, and ex1_config.r's DELTA_FRAC for the ex1 version of this argument.
##
## Anchored on mean sd(y_n), not RMS: an effect is a CHANGE in the series, so the yardstick is how
## much a typical unit moves, not how far it sits from zero.
DELTA_FRAC_EX2 <- 1.0
