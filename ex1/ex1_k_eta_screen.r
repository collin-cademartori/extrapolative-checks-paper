## Paired screen on ex1's stationary arm: which prior choices move the stat/nonstat contrast, and
## does the OVERFITTING statistic actually predict the DELTA error?
##
## PRIMARY COMPARISONS, fixed before this run:
##   (1) SIGMA LADDER, rho and K held at the committed values, sigma = 1.0 / 1.5 / 2.0 / 3.0:
##       is noise_abs_tr monotone in sigma, and is |delta| error monotone too, or U-shaped with its
##       minimum at the dispersion-matching sigma = 2.0? Monotone-together supports "prior
##       predictive sensibility -> overfitting -> worse estimation"; U-shaped says overfitting and
##       UNDERfitting both damage the counterfactual and noise_abs_tr sees only one of them.
##   (2) CROSS-DATASET CORRELATION within a single config, cor(noise_abs_tr, |delta| error) over
##       datasets, for the stationary arm. This is the mechanism at the level it should operate --
##       no configuration is compared to any other, so no single config can carry it. Measured at
##       n = 24 it was +0.25 to +0.36 in all five configs but never individually significant;
##       n = 62 resolves r = 0.25 at p < 0.05.
##   (3) C2 - C0 on the |z| gap, which at n = 24 was +0.062 (one-sided p 0.065). n = 38 resolves it.
## Everything else is secondary and labelled as such.
##
## WHY BOTH. The cross-config slope of the delta gap on the overfitting gap is +0.35 (se 0.24) over
## all five configs but flips to -0.67 with C2 removed, because C4 -- the only config with a TIGHTER
## signal prior -- has the lowest overfitting and the highest delta error. So the cross-config
## evidence cannot settle the question: it confounds overfitting with underfitting. (1) fixes that
## by moving sigma with rho held, so every rung is the same failure mode in different amounts, and
## (2) sidesteps configurations entirely.
##
## CONFIGURATIONS. C0 is the committed one. Only C0 and S3-at-2.0 reproduce the data's dispersion;
## the other ladder rungs are deliberately mis-scaled (prior predictive sd about 49%, 73% and 147%
## of the data's at sigma 1.0, 1.5, 3.0), so the ladder is a MECHANISM PROBE and not a menu of
## defensible configurations. That sigma = 2.0 is the unique dispersion-matching point is itself
## part of the argument.
##
## DATASETS come in two BATCHES from the same DGP, because the dataset draw depends on the rep count
## (test_data is drawn with iter = 2 * n) -- so raising n in one batch would have changed every
## existing dataset and thrown away 240 completed fits. Batch 1 reproduces the original 24 exactly;
## batch 2 adds 38 more from its own seed. Comparisons are within-dataset throughout, so pooling
## across batches is just pooling iid draws.
##
## USAGE (run from the ex1/ directory):
##     Rscript ex1_k_eta_screen.r [n_cores]
##     TIME_ONE=1 TIME_CFG=S3 Rscript ex1_k_eta_screen.r      # one dataset, one config, then stop
## Results accumulate per (batch, dataset, config) in ckpt_k_eta/ and combine into K_ETA_ALL.rds,
## so an interrupted run resumes rather than restarting.
suppressPackageStartupMessages({library(foreach); library(doParallel); library(purrr)})

.args <- commandArgs(trailingOnly = TRUE)
.cores <- suppressWarnings(as.integer(.args[1])); if (is.na(.cores)) .cores <- 4L
.nds <- suppressWarnings(as.integer(.args[2]))
TIME_ONE <- nzchar(Sys.getenv("TIME_ONE"))

## Same tmpfs precaution as ex1_sim_study.r; TMPDIR must be set before makeCluster.
scratch_root <- Sys.getenv("STUDY_SCRATCH", file.path(getwd(), ".scratch"))
dir.create(scratch_root, showWarnings = FALSE, recursive = TRUE)
if (!nzchar(Sys.getenv("CMDSTAN_OUTPUT_DIR")))
  Sys.setenv(CMDSTAN_OUTPUT_DIR = file.path(scratch_root, "cmdstan"))
if (!nzchar(Sys.getenv("STUDY_KEEP_TMPDIR"))) {
  tmp_dir <- file.path(scratch_root, "rtmp"); dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir)
}
cat(sprintf("  scratch: CMDSTAN_OUTPUT_DIR=%s  TMPDIR=%s\n",
  Sys.getenv("CMDSTAN_OUTPUT_DIR"), Sys.getenv("TMPDIR")))

## Master and workers run the SAME setup file; see ex1_k_eta_setup.r for why it is a file.
source("ex1_k_eta_setup.r")

## Two batches of datasets; see the header for why the rep count cannot simply be raised.
## Batch 1 reproduces the original 24 bit-for-bit (same seed, same count, so same draw).
BATCHES <- list(list(seed = 40318L, n = 24L), list(seed = 40319L, n = 38L))

## C0 is the committed configuration.
##   C1  K only            C2  K + tighter eta      C5  tighter eta only, at the committed K
##   C3  free rho          C4  free rho + its own sigma fixed point
##   S1/S2/S3  the sigma ladder, everything else at the committed values
##
## C5 is new: C2 moves K and eta together, so the eta effect was only reachable by subtracting C1.
## C5 measures it directly at the committed K.
CONFIGS <- list(
  C0 = list(K = 4L, eta = 0.10, rho = RHO_STAT, smult = 2.00),
  C1 = list(K = 6L, eta = 0.10, rho = RHO_STAT, smult = 2.00),
  C2 = list(K = 6L, eta = 0.05, rho = RHO_STAT, smult = 2.00),
  C3 = list(K = 4L, eta = 0.10, rho = c(1, 1),  smult = 2.00),
  C4 = list(K = 4L, eta = 0.10, rho = c(1, 1),  smult = 0.74),
  C5 = list(K = 4L, eta = 0.05, rho = RHO_STAT, smult = 2.00),
  S1 = list(K = 4L, eta = 0.10, rho = RHO_STAT, smult = 1.00),
  S2 = list(K = 4L, eta = 0.10, rho = RHO_STAT, smult = 1.50),
  S3 = list(K = 4L, eta = 0.10, rho = RHO_STAT, smult = 3.00),
  ## The HIGH end of the error ladder, which the sigma ladder and C2/C5 do not reach. The prior
  ## predictive check is ONE-SIDED in eta -- it imposes a ceiling and says nothing below it:
  ## P(|cor(y,t)| > 0.7) runs 0.412 / 0.400 / 0.306 / 0.195 at eta_frac 0.05 / 0.10 / 0.286 / 0.50
  ## against the DGP's 0.433. So 0.286 and 0.50 are values the check RULES OUT.
  ##   E1  eta 0.286 is the DGP's own error scale (0.286 x mean RMS ~ 1.77 against a truth of 2.0)
  ##       and the value the NONSTATIONARY arm uses. If overfitting falls from 0.361 toward
  ##       nonstat's 0.211 and the |z| gap closes here, then the error scale the check forbids is
  ##       exactly the safe one -- which is the example's thesis, demonstrated rather than argued.
  ##   E2  eta 0.50, well past the ceiling, to show the far end of the dose-response.
  E1 = list(K = 4L, eta = 0.286, rho = RHO_STAT, smult = 2.00),
  E2 = list(K = 4L, eta = 0.50,  rho = RHO_STAT, smult = 2.00),
  ## POINT-MASS eta at two locations. Estimating eta is a weak instrument (pass-through 0.24), so
  ## F1/F2 pin it exactly. F1 vs C0 and F2 vs E1 are the fixed-versus-estimated contrasts.
  ##
  ## What fixing COSTS, and it is worth stating: the prior-data conflict diagnostic goes dark. That
  ## diagnostic is what currently flags eta = 0.05 as indefensible (posterior at the 3e-5 tail in 22
  ## of 24 datasets), i.e. it is the signal that the posterior can catch a misspecification the prior
  ## predictive check missed. eta_prior_tail is recorded as NA for these two configs rather than
  ## reporting the constant that a point mass trivially produces.
  F1 = list(K = 4L, eta = 0.10,  rho = RHO_STAT, smult = 2.00, fix = TRUE),
  F2 = list(K = 4L, eta = 0.286, rho = RHO_STAT, smult = 2.00, fix = TRUE)
)
## Default every config to an ESTIMATED error scale unless it says otherwise.
CONFIGS <- lapply(CONFIGS, function(c_) { if (is.null(c_$fix)) c_$fix <- FALSE; c_ })
## C1 and C3 stay at batch 1 only: the K effect on delta is 0.003 (se 0.015) and the rho effect is
## 0.007 (se 0.022), both settled as null at n = 24, so more datasets would buy nothing. Everything
## else runs on both batches.
CONFIG_BATCHES <- list(C0 = 1:2, C1 = 1L, C2 = 1:2, C3 = 1L, C4 = 1:2,
                       C5 = 1:2, S1 = 1:2, S2 = 1:2, S3 = 1:2,
                       E1 = 1:2, E2 = 1:2, F1 = 1:2, F2 = 1:2)

stopifnot(CONFIGS$C0$K == K_LATENT, abs(CONFIGS$C0$eta - ETA_FRAC_STAT) < 1e-12,
          identical(CONFIGS$C0$rho, RHO_STAT), CONFIGS$C0$smult == SIGMA_MULT_STAT,
          setequal(names(CONFIGS), names(CONFIG_BATCHES)))

## ONLY=C0,E1,E2,F1,F2 restricts the run to those configs. Everything is checkpointed per
## (batch, dataset, config), so a machine that runs a subset writes cells another machine can read:
## copy ckpt_k_eta/ contents together and the summary sees one combined design.
.only <- Sys.getenv("ONLY")
if (nzchar(.only)) {
  keep <- trimws(strsplit(.only, ",")[[1]])
  stopifnot(all(keep %in% names(CONFIGS)))
  CONFIGS <- CONFIGS[keep]; CONFIG_BATCHES <- CONFIG_BATCHES[keep]
  cat(sprintf("  ONLY set: running %s\n", paste(keep, collapse = ", ")))
}
stopifnot(setequal(names(CONFIGS), names(CONFIG_BATCHES)))

## ---- datasets ------------------------------------------------------------------------------------
make_batch <- function(seed, n) {
  set.seed(seed)
  units <- sample.int(2 * n, size = n, replace = FALSE)
  pp_seed <- sample.int(.Machine$integer.max, 1)
  td <- sample_model(
    overall_scales = rep(DGP_SIGMA, N_UNITS), err_scale = DGP_ETA, absolute_error = TRUE,
    alpha_diag = ALPHA_DIAG, autocor_a = DGP_RHO[1], autocor_b = DGP_RHO[2],
    nonstationary = TRUE, num_treated = 0,
    type = "prior_pred", K_latent = K_LATENT, iter = 2 * n, seed = pp_seed)
  list(td = td, units = units)
}
batches <- lapply(BATCHES, function(b) make_batch(b$seed, b$n))
tasks <- do.call(rbind, lapply(seq_along(batches),
  function(bi) data.frame(bi = bi, unit = batches[[bi]]$units)))
N_TOTAL <- nrow(tasks)

cat(sprintf("\n=== %d datasets (batch 1: %d, batch 2: %d) from the study DGP ===\n",
  N_TOTAL, BATCHES[[1]]$n, BATCHES[[2]]$n))
ds_summ <- do.call(rbind, lapply(seq_len(N_TOTAL), function(k) {
  y <- batches[[tasks$bi[k]]]$td$ys[tasks$unit[k], , ]
  data.frame(batch = tasks$bi[k], unit = tasks$unit[k],
             rms = round(mean(apply(y, 2, function(v) sqrt(mean(v^2)))), 2),
             sd = round(mean(apply(y, 2, sd)), 2))
}))
cat(sprintf("  distinct datasets: %d of %d   mean RMS %.2f   mean sd %.2f\n",
  nrow(unique(ds_summ[, c("rms", "sd")])), N_TOTAL, mean(ds_summ$rms), mean(ds_summ$sd)))
cat(sprintf("  per batch: RMS %.2f / %.2f   sd %.2f / %.2f  (same DGP, so these should agree)\n",
  mean(ds_summ$rms[ds_summ$batch == 1]), mean(ds_summ$rms[ds_summ$batch == 2]),
  mean(ds_summ$sd[ds_summ$batch == 1]),  mean(ds_summ$sd[ds_summ$batch == 2])))
cat(sprintf("  cells to fill: %d\n", sum(sapply(names(CONFIGS),
  function(cn) sum(tasks$bi %in% CONFIG_BATCHES[[cn]])))))

fit_one <- function(bi, s, cfg_name) {
  cfg <- CONFIGS[[cfg_name]]
  env <- new.env(parent = globalenv())
  assign("ETA_FRAC_STAT", cfg$eta, envir = env)
  assign("RHO_STAT", cfg$rho, envir = env)
  assign("SIGMA_MULT_STAT", cfg$smult, envir = env)
  assign("STAT_FIX_ETA", isTRUE(cfg$fix), envir = env)
  f <- run_sim_stat; environment(f) <- env
  ## Batch 1 keeps its ORIGINAL seed rule, so the configs added in this run share common random
  ## numbers with the ones already checkpointed. Changing it would silently unpair them.
  set.seed(if (bi == 1L) 90000 + s else 200000 + s)
  r <- f(batches[[bi]]$td, s, K_latent = cfg$K, progress_log = NULL)
  keep <- c("noise_abs_tr", "pred_perc", "pred_width", "eta_med", "eta_prior_tail",
            "time_cor_pval", "pred_mad", "rhat_M", "rhat_max", "n_div", "ess_delta1", "mean_1", "absz_1")
  ## What the fit ACTUALLY resolved, not what was requested -- see the note above.
  out <- data.frame(batch = bi, ds = s, config = cfg_name, K = cfg$K, eta_frac = cfg$eta,
                    eta_frac_used = evalq(ETA_FRAC_STAT, environment(f)),
                    rho_a_used = evalq(RHO_STAT, environment(f))[1],
                    rho_b_used = evalq(RHO_STAT, environment(f))[2],
                    smult_used = evalq(SIGMA_MULT_STAT, environment(f)))
  for (arm in c("nonstat", "stat")) for (k in keep) {
    v <- r[[paste0(arm, "_", k)]]
    out[[paste0(arm, "_", k)]] <- if (is.null(v)) NA_real_ else as.numeric(v)
  }
  for (nm in grep("^stat_rho_", names(r), value = TRUE)) out[[nm]] <- as.numeric(r[[nm]])
  ## A point mass has no tail to sit in; the recorded number would be a constant, not a diagnostic.
  if (isTRUE(cfg$fix)) out$stat_eta_prior_tail <- NA_real_
  out$stat_fix_eta <- isTRUE(cfg$fix)
  out
}

if (TIME_ONE) {
  t0 <- Sys.time(); x <- fit_one(1L, batches[[1]]$units[1], Sys.getenv("TIME_CFG", "C2"))
  cat(sprintf("\none dataset, one config (2 fits) took %.1f min\n",
    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  print(t(x)); quit(save = "no")
}

CK <- file.path(getwd(), "ckpt_k_eta"); dir.create(CK, showWarnings = FALSE)
cat(sprintf("\n=== %d cells outstanding, %d workers ===\n",
  sum(sapply(seq_len(N_TOTAL), function(k) sum(sapply(names(CONFIGS), function(cn)
    tasks$bi[k] %in% CONFIG_BATCHES[[cn]] &&
    !file.exists(file.path(CK, sprintf("b%d_ds%03d_%s.rds", tasks$bi[k], tasks$unit[k], cn))))))),
  .cores))
cl <- makeCluster(.cores, outfile = ""); registerDoParallel(cl)
invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({
  library(cmdstanr); library(posterior); library(forcats); library(dplyr); library(ggplot2); library(purrr)})))
invisible(clusterCall(cl, source, "ex1_k_eta_setup.r"))
t0 <- Sys.time()
res <- foreach(k = seq_len(N_TOTAL), .combine = dplyr::bind_rows,
  .export = c("batches", "tasks", "CONFIGS", "CONFIG_BATCHES", "fit_one", "CK"),
  .packages = c("cmdstanr", "posterior", "purrr", "dplyr", "ggplot2")) %dopar% {
  bi <- tasks$bi[k]; s <- tasks$unit[k]
  tf <- Sys.time()
  ## Checkpointed per (batch, dataset, config), so a config that has to be refitted -- as C2 did,
  ## after the ETA_FRAC_STAT plumbing bug -- costs only its own cells.
  out <- dplyr::bind_rows(lapply(names(CONFIGS), function(cn) {
    if (!(bi %in% CONFIG_BATCHES[[cn]])) return(NULL)
    ck <- file.path(CK, sprintf("b%d_ds%03d_%s.rds", bi, s, cn))
    if (file.exists(ck)) return(readRDS(ck))
    r <- tryCatch(fit_one(bi, s, cn), error = function(e) {
      cat(sprintf("  b%d ds %d %s FAILED: %s\n", bi, s, cn, conditionMessage(e))); NULL })
    if (!is.null(r)) saveRDS(r, ck)
    r
  }))
  cat(sprintf("  b%d ds %d done, %d configs, %.1f min\n", bi, s,
    length(unique(out$config)), as.numeric(difftime(Sys.time(), tf, units = "mins"))))
  out
}
stopCluster(cl)
saveRDS(res, "K_ETA_ALL.rds")
cat(sprintf("\n%d rows in %.1f min -> K_ETA_ALL.rds\n", nrow(res),
  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
