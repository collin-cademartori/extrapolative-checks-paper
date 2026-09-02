## Setup shared by the master process and every worker of ex1_k_eta_screen.r.
##
## A FILE rather than a string passed to clusterCall. The first version carried this code in a
## quoted string and the nested escaping broke on the first inner double quote; sourcing a file has
## no escaping level at all. It is also the study's own pattern -- ex1_sim_study.r does
## clusterCall(cl, source, "ex1_config.r") for the same reason, so that a worker cannot end up
## running different constants than the master.

source("../sample_model.r"); source("../plotting.r"); source("../pathfinder_init.r")
source("ex1_config.r")

## run_sim_stat, anchor_order and worker_progress are lifted out of ex1_sim_study.r by line range
## rather than copied, so this screen cannot describe a fit the study does not perform. Bounds are
## resolved explicitly and checked: an unchecked grep() pair silently produced a BACKWARDS range
## once already, which extracts nothing and fails far from the cause.
src <- readLines("ex1_sim_study.r")
grab <- function(from, to, tweak = identity) {
  a <- grep(from, src); stopifnot(length(a) == 1)
  b <- grep(to, src); b <- b[b > a][1]; stopifnot(!is.na(b), b - 1 > a)
  eval(parse(text = tweak(paste(src[a:(b - 1)], collapse = "\n"))), envir = globalenv())
}
## Each substitution asserts it matched exactly once, so a drift in the study fails loudly here
## rather than silently leaving the screen measuring something else.
sub1 <- function(t, o, n) {
  stopifnot(sum(gregexpr(o, t, fixed = TRUE)[[1]] > 0) == 1)
  sub(o, n, t, fixed = TRUE)
}

grab("^worker_progress <- function", "^# Escalation ladder|^EX1_LADDER")
grab("^anchor_order <- function", "^run_sim_stat <- function")
## Two surgical edits to the stationary arm only: return its rho draws, and record the per-factor
## posterior median. Needed because the Beta(1,1) configs ask whether the factors SPLIT -- some
## staying smooth to carry the trend, others going rough to absorb error -- which the aggregate
## metrics cannot show.
grab("^run_sim_stat <- function", "^run_sim_study_stat <- function", function(t) {
  t <- sub1(t,
    '      n_chains = 3, pathfinder_init = TRUE\n    ),\n    seeds = fit_seeds[2, ]',
    '      n_chains = 3, pathfinder_init = TRUE, return_draws = "rho"\n    ),\n    seeds = fit_seeds[2, ]')
  ## Point-mass option for the stationary arm's error scale, mirroring the study's own
  ## NONSTAT_FIX_ETA branch: tau_val > 0 fixes err_sd in ife_named.stan, m_tau/s_tau > 0 estimates
  ## it. Needed because the ESTIMATED version is a weak instrument -- halving the prior location
  ## moves the realised eta only 12% (pass-through 0.24), so a config cannot deliver the error scale
  ## it nominally asks for. Fixing makes the manipulation exact.
  t <- sub1(t,
    '      overall_scales = overall_scales_stat, err_scale = 0, absolute_error = TRUE,\n      err_scale_mean = ETA_FRAC_STAT * eta_anchor,\n      err_scale_sd = ETA_FRAC_STAT * eta_anchor,',
    '      overall_scales = overall_scales_stat, absolute_error = TRUE,\n      err_scale      = if (STAT_FIX_ETA) ETA_FRAC_STAT * eta_anchor else 0,\n      err_scale_mean = if (STAT_FIX_ETA) 0 else ETA_FRAC_STAT * eta_anchor,\n      err_scale_sd   = if (STAT_FIX_ETA) 0 else ETA_FRAC_STAT * eta_anchor,')
  sub1(t, "      return(res)", paste0(
    '      if (!is.null(pfit$draws)) {\n',
    '        rq <- apply(posterior::as_draws_matrix(pfit$draws), 2, median)\n',
    '        res[paste0("rho_", seq_along(rq))] <- as.list(sort(rq, decreasing = TRUE))\n',
    '      }\n      return(res)'))
})

## Sampler settings. Longer than the study's fast mode (500/500): K = 6 is more weakly identified
## and rhat_M is what certifies delta, so buying mixing here is the difference between a paired
## contrast and a noise measurement. No escalation ladder -- rhat is recorded per fit and reported.
EX1_ITER <- 1500L; EX1_WARM <- 1000L
EX1_LADDER <- escalation_ladder(integer(0), integer(0))
ESCALATE_MAX <- EX1_LADDER$max_rounds
PLOT_REPS <- 0L                      # no per-rep figures from a screen
STAT_FIX_ETA <- FALSE                # default; each config overrides it in fit_one's enclosure
