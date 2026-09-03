## Paired factorial screen: does the LEVEL OFFSET or the EFFECT PRIOR change the ints - no_ints gap
## on cor_sq(spurious)?
##
## PRIMARY COMPARISONS, fixed before running:
##   (1) paired difference in (ints - no_ints) between L = 1 and L = 10
##   (2) paired difference in (ints - no_ints) between DF = 0.5 and DF = 2.0
## Everything else is secondary and labelled as such.
##
## DESIGN. Between-dataset sd of cor_sq(spurious) is ~0.07, which is the size of the effects in
## question -- which is why every earlier unpaired screen was uninformative. So every dataset is
## measured under ALL FOUR conditions and the comparisons are within-dataset:
##   L  -- common random numbers. f_alt = f_treat - L cancels the offset and per-unit sd is
##         level-invariant, so changing L shifts unit n by exactly L * Lambda[n,1] and nothing else.
##   DF -- fit-side, so literally the same data under two priors.
##
## SETTINGS. 500/500 rather than the study's 1500/1000, deliberately: measured on the fast study,
## rhat_cor_sq had median 1.005 / max 1.014 and ess ~1000, so the Monte Carlo error on a cor_sq
## posterior mean is ~0.003 -- twenty times smaller than the between-dataset sd this design has to
## beat. Iterations are the wrong place to spend the budget here; datasets are the right place.
## USAGE (run from the ex2/ directory):
##     Rscript ex2_level_delta_screen.r [n_cores] [n_datasets]
##     SUMMARY_ONLY=1 Rscript ex2_level_delta_screen.r      # print dataset summaries and stop
## Results accumulate in ckpt_level_delta/ (one .rds per dataset) and are combined into ALL.rds, so
## an interrupted run resumes rather than restarting. Delete that directory to force a fresh run.
suppressPackageStartupMessages({library(foreach); library(doParallel)})

.args <- commandArgs(trailingOnly = TRUE)
.cores <- suppressWarnings(as.integer(.args[1])); if (is.na(.cores)) .cores <- 4L
.nds <- suppressWarnings(as.integer(.args[2]))

## Keep CmdStan's CSVs and R's temp staging OFF /tmp, which is a TMPFS on the study machine -- the
## failure mode that killed two long runs and, in this session, three scratch harnesses. Same
## mechanism as ex2_sim_study.r; see the long note there. TMPDIR must be set before makeCluster so
## the PSOCK workers inherit it.
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
source("../sample_model.r"); source("../pathfinder_init.r"); source("ex2_config.r")
src <- readLines("ex2_sim_study.r")
eval(parse(text = paste(src[grep("^anchor_order <- function", src):(grep("^unpermute_untreated <- function", src)-1)], collapse="\n")))
eval(parse(text = paste(src[grep("^unpermute_untreated <- function", src):(grep("^unpermute_untreated <- function", src)+5)], collapse="\n")))
ruv <- function(d) { v <- rnorm(d); v/sqrt(sum(v*v)) }

N_DATASETS <- if (!is.na(.nds)) .nds else 30L
LEVELS <- c(1, 10)
DFS <- c(0.5, 2.0)
SIM <- 0.9; N_SPUR <- 2
SUMMARY_ONLY <- nzchar(Sys.getenv("SUMMARY_ONLY"))

## ---- generate all datasets up front, each from its own explicit seed ---------------------------
make_components <- function(seed) {
  set.seed(seed)
  Tt <- DGP_T_TIMES; N_unc <- DGP_N_UNITS - 1 - DGP_N_COMP_TRUE - N_SPUR; K <- 2 + DGP_K_UNC
  ar_treat <- as.numeric(arima.sim(model=list(ar=0.9), n=Tt))
  f_alt <- ar_treat + c(rep(0, Tt-DGP_T_TREATED), rep(-DGP_F_TREAT_SD, DGP_T_TREATED))
  f_unc <- matrix(nrow=DGP_K_UNC, ncol=Tt); cu <- Inf
  while (cu > 0.01) { for (k in seq_len(DGP_K_UNC)) f_unc[k,] <- rnorm(1,0,2)+arima.sim(model=list(ar=0.9),n=Tt)
    cu <- max(abs(cor(t(f_unc), t(t(ar_treat))))) }
  lo <- matrix(nrow=DGP_N_UNITS, ncol=K); lo[1,] <- c(1, rep(0,K-1))
  for (n in seq_len(DGP_N_COMP_TRUE)) lo[1+n,] <- c(sqrt(SIM),0,sqrt(1-SIM)*ruv(K-2))
  for (n in seq_len(N_SPUR)) lo[1+DGP_N_COMP_TRUE+n,] <- c(0,sqrt(SIM),sqrt(1-SIM)*ruv(K-2))
  for (n in seq_len(N_unc)) lo[1+DGP_N_COMP_TRUE+N_SPUR+n,] <- c(0,0,ruv(K-2))
  lat0 <- lo %*% rbind(ar_treat, f_alt, f_unc)
  noise <- matrix(rnorm(DGP_N_UNITS*Tt, sd=DGP_NOISE_FRAC*mean(apply(lat0,1,sd))), DGP_N_UNITS, Tt)
  list(lat0=lat0, lo=lo, noise=noise, Tt=Tt)
}
build_Y <- function(cp, L) t(cp$lat0 + L * outer(cp$lo[,1], rep(1, cp$Tt)) + cp$noise)

comps <- lapply(seq_len(N_DATASETS), function(i) make_components(70000 + i))

## ---- summaries, computed SEPARATELY AT EACH LEVEL so the shift is verified, not assumed --------
summ <- do.call(rbind, lapply(seq_len(N_DATASETS), function(i)
  do.call(rbind, lapply(LEVELS, function(L) {
    Y <- build_Y(comps[[i]], L); m <- colMeans(Y)
    data.frame(ds=i, L=L, treated_mean=round(m[1],2), grand=round(mean(Y),2),
      spur_gap=round(mean(abs(m[4:5]-m[1])),2), sd_colmeans=round(sd(m),2),
      mean_sd_y=round(mean(apply(Y,2,sd)),3), mean_rms=round(mean(apply(Y,2,function(y) sqrt(mean(y^2)))),2))
  }))))
cat("=== dataset summaries, at each level ===\n")
print(summ[order(summ$L, summ$ds), ], row.names = FALSE)

cat("\n=== validity checks ===\n")
for (L in LEVELS) {
  keys <- sapply(seq_len(N_DATASETS), function(i) paste(round(colMeans(build_Y(comps[[i]], L)), 8), collapse=","))
  cat(sprintf("  L=%-3d  distinct datasets: %d of %d\n", L, length(unique(keys)), N_DATASETS))
  stopifnot(length(unique(keys)) == N_DATASETS)
}
d1 <- summ[summ$L==LEVELS[1],]; d2 <- summ[summ$L==LEVELS[2],]
cat(sprintf("  level shift applied: treated_mean rises in %d of %d datasets (mean rise %.2f)\n",
    sum(d2$treated_mean > d1$treated_mean), N_DATASETS, mean(d2$treated_mean - d1$treated_mean)))
cat(sprintf("  spur_gap:  L=%d mean %.2f  ->  L=%d mean %.2f\n",
    LEVELS[1], mean(d1$spur_gap), LEVELS[2], mean(d2$spur_gap)))
cat(sprintf("  mean_sd_y (should be UNCHANGED by L): %.3f vs %.3f  -> max |diff| %.4f\n",
    mean(d1$mean_sd_y), mean(d2$mean_sd_y), max(abs(d1$mean_sd_y - d2$mean_sd_y))))
if (SUMMARY_ONLY) { cat("\nSUMMARY_ONLY set -- stopping before the fits.\n"); quit(save="no") }

## ---- fits: 30 datasets x 2 L x 2 DF x 2 arms = 240, checkpointed per dataset -------------------
CK <- file.path(getwd(), "ckpt_level_delta"); dir.create(CK, showWarnings=FALSE)
cat(sprintf("\n=== %d datasets x %d levels x %d effect priors x 2 arms = %d fits, %d workers ===\n",
  N_DATASETS, length(LEVELS), length(DFS), N_DATASETS * length(LEVELS) * length(DFS) * 2, .cores))
cl <- makeCluster(.cores, outfile=""); registerDoParallel(cl)
invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, {suppressPackageStartupMessages({library(cmdstanr); library(posterior)})}))
t0 <- Sys.time()
res <- foreach(i = seq_len(N_DATASETS), .combine=rbind,
  .export=c("comps","build_Y","anchor_order","unpermute_untreated","sample_model","ife_mod",
            "pathfinder_inits","draw_to_init","PF_PARAM_BASES","LEVELS","DFS","CK",
            "DGP_N_UNITS","DGP_T_TIMES","K_LATENT","ETA_FRAC_EX2","ETA_CV_EX2","RHO_EX2",
            "NUM_TREATED","ALPHA_DIAG","INT_FRAC"),
  .packages=c("cmdstanr","posterior")) %dopar% {
  ck <- file.path(CK, sprintf("ds_%03d.rds", i)); if (file.exists(ck)) return(readRDS(ck))
  out <- NULL
  for (L in LEVELS) {
    Y <- build_Y(comps[[i]], L); perm <- anchor_order(Y, K_LATENT); fy <- Y[, perm]
    sd_y <- apply(fy,2,sd); ea <- mean(sd_y); rms <- apply(fy,2,function(x) sqrt(mean(x^2)))
    for (DF in DFS) for (arm in c("no_ints","ints")) {
      a <- list(N_units=DGP_N_UNITS, T_times=DGP_T_TIMES, K_latent=K_LATENT,
        overall_scales = if (arm=="no_ints") rms else sd_y,
        err_scale=0, absolute_error=TRUE, err_scale_mean=ETA_FRAC_EX2*ea,
        err_scale_sd=ETA_CV_EX2*ETA_FRAC_EX2*ea, data=fy,
        autocor_a=RHO_EX2[1], autocor_b=RHO_EX2[2], nonstationary=FALSE,
        num_treated=NUM_TREATED, delta_scale=DF*ea, fit_scales=0, alpha_diag=ALPHA_DIAG,
        pathfinder_init=TRUE, type="posterior", quiet=TRUE, ad=0.8,
        iter=500, iter_warm=500, n_chains=4, parallel_chains=2, seed=800000 + 1000*i + 10*L + DF*10)
      if (arm=="no_ints") a$include_factor_means <- TRUE else {
        a$include_ints <- TRUE; a$int_scale <- INT_FRAC*sd(colMeans(fy)); a$int_loc <- mean(fy) }
      f <- do.call(sample_model, a)
      cs <- unpermute_untreated(f$cor_sq, perm)
      out <- rbind(out, data.frame(ds=i, L=L, DF=DF, arm=arm, dscale=round(DF*ea,3),
        true=mean(cs[1:2]), spurious=mean(cs[3:4]), unc=mean(cs[5:7]),
        mae_delta=mean(abs(f$effect_means)), rhat_cor_sq=f$sampler_diag$rhat_cor_sq,
        rhat_M=f$sampler_diag$rhat_M, n_div=f$sampler_diag$n_div))
    }
  }
  saveRDS(out, ck); out
}
stopCluster(cl)
saveRDS(res, file.path(CK, "ALL.rds"))
cat(sprintf("\n=== %d fits in %.1f min ===\n", nrow(res), as.numeric(difftime(Sys.time(), t0, units="mins"))))
