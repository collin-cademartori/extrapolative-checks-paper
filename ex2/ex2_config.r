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

## SWEPT AXIS 1: the level gap. f_treat = level + AR(0.9) and f_alt = f_treat - level, so this IS
## the gap between the treated/true group's long-run average and the spurious group's. It is the
## quantity that makes location informative about correlation, which is exactly the entanglement
## the unit-intercepts model assumes away -- so it belongs on the sweep rather than being a fixed
## literal, which is what it was (a bare 6, appearing twice inside sim_model_intercepts).
##
## MEASURED by ex2_level_delta_screen.r: 30 datasets, each fitted under both levels and both effect
## priors under common random numbers. Posterior mean cor_sq with the SPURIOUS units:
##
##     level gap      no_ints        ints          gap (ints - no_ints)
##        1            0.635         0.684              0.049
##       10            0.469         0.689              0.220
##
##     level effect on the gap:  +0.171   se 0.021   95% CI [+0.130, +0.212]
##
## no_ints carries the whole effect; ints does not move. no_ints has no gamma, so a unit's level
## must come out of Lambda * Phi_means -- the same loadings that generate correlation -- and
## separating levels forces Lambda_spurious away from Lambda_treated. ints puts the level in gamma
## and is immune. The same mechanism shows on the other side: over the same range no_ints' cor_sq
## with the TRUE comparators RISES +0.101 (se 0.008) while ints moves +0.001 (se 0.001).
##
## 5 and 10, replacing the old fixed 6. 5 keeps a moderate cell close to what the study used before
## so earlier runs stay roughly comparable; 10 is where the effect was actually measured.
DGP_LEVELS <- c(5, 10)

## SWEPT AXIS 2 is N_comp, the number of spurious comparators (set in ex2_sim_study.r's call).
##
## The `sim` axis it REPLACED -- the strength of the spurious units' pre-treatment correlation,
## swept over {0.7, 0.9} -- was dropped because it does almost nothing. Measured on the archived
## 200-rep 2x2 (200 reps per cell), main effect of each axis on the ints - no_ints contrast:
##
##                            sim 0.7 -> 0.9        num_comp 1 -> 3
##     |delta| error gap     -0.0002  p 0.97       +0.0099  p 0.084
##     |z| gap               -0.0002  p 0.99       +0.0286  p 0.074
##     acor_err gap          -0.0047  p 0.010      -0.0137  p <0.0001
##
## sim is indistinguishable from zero on the estimand and three times weaker where both register,
## while the level gap above dwarfs them both. So sim is now FIXED at the stronger of its two old
## values. (Caveat: that 2x2 predates absolute error mode, the shared delta scale and the anchored
## intercept prior, so it ranks the axes rather than measuring them under the current config.)
DGP_SIM <- 0.9

## Observation-noise sd, as a fraction of the AVERAGE unit's latent sd.
##
## It was 0.1 * max_n sd(latent_n), flagged as unstable because the realised noise sd swings about
## 2x across datasets (5-95%: 0.133 to 0.296, mean 0.202, CV 0.25). That flag was half right and its
## CONSEQUENCE was overstated. Measured over 1500 datasets on the sweep grid:
##
##   * `max` is barely less stable than `mean` -- CV 0.25 against 0.22. The swing is not an artifact
##     of using an order statistic; the per-unit latent sds themselves vary that much, because the
##     AR realisations and the random loadings do.
##   * and it very largely CANCELS. The true noise sd and mean sd(y_n) -- the anchor ETA_FRAC_EX2 is
##     built on -- correlate at +0.93, so the RATIO between them, which is what ETA_FRAC_EX2 claims
##     to be, has CV 0.09 rather than 0.25. The prior located at 0.12 * mean sd(y_n) lands within
##     88-116% of the true noise sd in 90% of datasets, against a prior whose own sd equals its
##     location (CV 1.0). It is calibrated per dataset an order of magnitude better than its width.
##
## So the anchoring is doing real work: it holds the signal-to-noise ratio roughly constant across
## datasets. Fixing the noise at a constant (ex1's approach) would REMOVE that and make the SNR
## swing 2x instead, which is worse for a study comparing arms across datasets.
##
## `mean` rather than `max` anyway, for consistency with every other anchor here and because "12% of
## the average unit's latent variation" is easier to state than "of the largest unit's". The
## fraction rises 0.100 -> 0.121 to hold the realised noise sd at its previous 0.202, so this is a
## change of anchor, not of noise level.
DGP_NOISE_FRAC <- 0.121

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
##
## AVERAGED OVER THE SWEEP, so it now covers both level cells rather than one fixed offset. It moved
## 2.05 -> 2.50 when the level gap went from a fixed 6 to a swept {5, 10}: the numerator scales with
## the level gap and the mean gap rose to 7.5. ETA_FRAC_EX2 did not move, because mean sd(y_n) is
## invariant to adding a constant to a unit's level. The guard in ex2_derive_scales.r caught this
## automatically when the sweep changed -- which is what it is for.
LEVEL_SPREAD_FRAC <- 2.50

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
##
## The VALUE is not doing much work. ex2_level_delta_screen.r brackets it, fitting the same 30
## datasets under 0.5 and 2.0: widening the prior raises the ints - no_ints gap on cor_sq(spurious)
## by +0.023 (se 0.004), real but roughly seven times smaller than the level effect above. It costs
## accuracy on the estimand it prices -- mean |delta| error rises 0.221 -> 0.293 (no_ints) and
## 0.325 -> 0.537 (ints) over that range -- so there is no case for widening past 1.0.
##
## 0.5, matching ex1's DELTA_FRAC and for the same reason: an effect of a full outcome SD would be
## enormous in most applied settings, so half an SD is the honest elicitation. Justified as a belief
## about effect sizes, NOT by what it does to the study's contrast -- the true effect is zero here,
## so any tightening of a zero-centred prior flatters both arms automatically. Measured, widening it
## slightly INCREASES the ints - no_ints gap (see above), so this change costs a little contrast and
## is made anyway because the elicitation is what should decide it.
DELTA_FRAC_EX2 <- 0.5
