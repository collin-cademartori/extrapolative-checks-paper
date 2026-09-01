## Companion to ex1_rsq_priors.r: the implied prior predictive distribution of the REALISED
## standard deviation of a series, under two ways of setting sigma.
##
## Like that figure, this one contains no data, by design. A prior predictive check in this paper's
## sense asks whether a configuration implies what the analyst EXPECTS, before any data exists; the
## reference here is an elicited expected SD, written S, not a measured one.
##
## What it shows: sigma is the LONG-RUN marginal SD of a unit (verified -- at T = 4000 the realised
## sd approaches sigma times E||Lambda[n,]|| = 0.95). Over a short window a near-unit-root AR(1)
## realises far less than that, because the sample mean absorbs the low-frequency wandering. So
## setting sigma to the SD one expects to see implies series far flatter than intended. The check
## catches that before any data arrives, which is the whole claim.
##
## Stated without reference to any anchor: the study's stationary configuration is eta = 0.10 x
## sigma (ETA_FRAC_STAT / SIGMA_MULT_STAT). That is deliberate here -- the study anchors its scales
## on RMS(y_n) measured per dataset, but this figure has no dataset, and the calibration argument
## does not need one. It is a statement about the model, not about any data.
##
## Read alongside ex1_rsq_priors.r rather than separately. The two statistics -- dispersion here,
## trendiness there -- do NOT pin one hyperparameter each: sd(y)^2 ~ sd(signal)^2 + eta^2, so both
## depend on sigma and eta together. They identify the PAIR. Measured, the SD-matching multiple
## moves only 1.98 -> 1.85 as ETA_FRAC_STAT goes 0.05 -> 0.2, so the pair is well identified over
## that range rather than balanced on a knife edge.

source("../sample_model.r")
source("../plotting.r")
source("ex1_config.r")

seed <- 81207
set.seed(seed)
pp_seed_naive <- sample.int(.Machine$integer.max, 1)
pp_seed_calib <- sample.int(.Machine$integer.max, 1)

# Work in units where the analyst's expected SD is S = 1.
S <- 1
# The study's stationary configuration, expressed without an anchor.
ETA_OVER_SIGMA <- ETA_FRAC_STAT / SIGMA_MULT_STAT

# SD_PER_SIGMA comes from ex1_config.r and is derived and guarded by ex1_derive_scales.r, so this
# figure cannot silently go stale if T_TIMES, K_LATENT, ALPHA_DIAG, RHO_STAT or ETA_FRAC_STAT move.
SIGMA_CALIB <- S / SD_PER_SIGMA   # ~3.1 x the expected SD

prior_sds <- function(sigma, pp_seed) {
  d <- sample_model(
    N_units = N_UNITS, T_times = T_TIMES, K_latent = K_LATENT,
    overall_scales = rep(sigma, N_UNITS), alpha_diag = ALPHA_DIAG,
    err_scale = ETA_OVER_SIGMA * sigma, absolute_error = TRUE,
    data = NULL,
    autocor_a = RHO_STAT[1], autocor_b = RHO_STAT[2],
    nonstationary = FALSE, num_treated = 0,
    type = "prior_pred", iter = 6000, n_chains = 1, seed = pp_seed, quiet = TRUE
  )
  apply(d$ys[, , 1], 1, sd)
}

sds <- list(
  `sigma = expected SD` = prior_sds(S, pp_seed_naive),
  `sigma = 3.1 x expected SD` = prior_sds(SIGMA_CALIB, pp_seed_calib)
)

for (nm in names(sds)) {
  cat(sprintf("%-28s E[realised sd] / S = %.2f\n", nm, mean(sds[[nm]]) / S))
}

# The expected SD is drawn as a reference line, not as data: it is the elicited quantity the
# configuration has to reproduce, and the panel is read by asking whether the mass sits on it.
plot_prior_sd <- function(sds, expected) {
  df <- stack(sds)
  names(df) <- c("sd", "config")

  ggplot(df) +
    facet_wrap(vars(config), ncol = 1, scales = "free_y", strip.position = "left") +
    geom_histogram(aes(x = sd), fill = "#aaaaaa", color = "#aaaaaa", bins = 40) +
    geom_vline(xintercept = expected, linetype = "dashed") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(strip.background = element_rect(fill = "white", color = "black")) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    xlab("Realized Standard Deviation\n (relative to the expected SD)") +
    ylab("") +
    ggtitle("(B)")
}

sd_hists <- plot_prior_sd(sds, S)
ggsave(sd_hists, file = "../figs/stat_sd_hists.pdf", device = "pdf",
  width = 2.5, height = 5, create.dir = TRUE)
