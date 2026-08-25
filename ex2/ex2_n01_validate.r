## Moderate validation (NOT shipped): does the normal(0,1) sigma_raw prior recover the ex2 science while
## keeping est-sigma, and are the high rhats minor-mode contamination (fixable) or genuine multimodality
## (systematic fix)? normal(0,1) model for ALL fits (no per-task switching -> no env bug). Both models,
## sim{0.7,0.9} x N_comp{2,3} x 2 reps. 6 chains, 800 warmup/600 -- longer, to separate transient from
## persistent. Captures science (cs by group, absz, multiplier), min_diag (zero-avoidance check), and the
## mode signature: per-chain lp, offmode, lp_gap, and rhat over ON-MODE chains vs ALL. Crash-safe rows.
library(foreach); library(doParallel); library(posterior)
source("../sample_model.r"); source("../pathfinder_init.r")
# CRITICAL: set ife_mod in the MASTER so sample_model (whose closure is this global env) actually
# uses the intended model in %dopar%. A clusterEvalQ override alone is IGNORED -- sample_model
# resolves ife_mod through its master closure, not the worker global. See diag: worker-global was
# correct but sample_model_sees was ife_named. So compile here and export ife_mod to the workers.
ife_mod <- cmdstan_model("../ife_named_unitnorm_ig_s1.stan")   # normal(0,1) sigma_raw
stopifnot(basename(ife_mod$stan_file()) == "ife_named_unitnorm_ig_s1.stan")

ruv <- function(d) { v <- rnorm(d); v / sqrt(sum(v * v)) }
anchor_order <- function(y, K) {
  N <- ncol(y); yc <- scale(y, center = TRUE, scale = FALSE); sel <- 1L
  remaining <- setdiff(seq_len(N), sel)
  while (length(sel) < K && length(remaining) > 0) {
    Q <- qr.Q(qr(yc[, sel, drop = FALSE]))
    resid <- yc[, remaining, drop = FALSE] - Q %*% crossprod(Q, yc[, remaining, drop = FALSE])
    pick <- remaining[which.max(colSums(resid^2))]; sel <- c(sel, pick); remaining <- setdiff(remaining, pick)
  }
  c(sel, remaining)
}
unpermute_untreated <- function(v, perm) { out <- numeric(length(v)); out[perm[-1] - 1] <- v; out }
sim_model_intercepts <- function(N_unc = 2, N_comp_true = 2, N_comp_spur = 2, T_times = 30, T_treated = 5, K_unc = 1, sim = 0.9) {
  N_units <- 1 + N_comp_true + N_comp_spur + N_unc; K_gen <- 2 + K_unc
  f_treat <- 6 + arima.sim(model = list(ar = 0.9), n = T_times); f_treat_sd <- 1.9
  f_alt <- (f_treat - 6) + c(rep(0, T_times - T_treated), rep(-f_treat_sd, T_treated))
  f_unc <- matrix(nrow = K_unc, ncol = T_times); cor_unc <- Inf
  while (cor_unc > 0.01) { for (k in 1:K_unc) f_unc[k, ] <- rnorm(1, 0, 2) + arima.sim(model = list(ar = 0.9), n = T_times); cor_unc <- max(abs(cor(t(f_unc), t(t(f_treat))))) }
  facs <- rbind(f_treat, f_alt, f_unc); loads <- matrix(nrow = N_units, ncol = K_gen); loads[1, ] <- c(1, rep(0, K_gen - 1))
  for (n in seq_len(N_comp_true)) loads[1 + n, ] <- c(sqrt(sim), 0, sqrt(1 - sim) * ruv(K_gen - 2))
  for (n in seq_len(N_comp_spur)) loads[1 + N_comp_true + n, ] <- c(0, sqrt(sim), sqrt(1 - sim) * ruv(K_gen - 2))
  for (n in seq_len(N_unc)) loads[1 + N_comp_true + N_comp_spur + n, ] <- c(0, 0, ruv(K_gen - 2))
  lat <- loads %*% facs
  list(Y = t(lat + rnorm(nrow(lat) * ncol(lat), sd = 0.1 * max(apply(lat, 1, sd)))), groups = c("treated", rep("true", N_comp_true), rep("spurious", N_comp_spur), rep("uncorrelated", N_unc)))
}

K_LAT <- 3L; NPT <- 5L; run_seed <- 88213; NCH <- 6L
GRID <- expand.grid(rep = 1:2, N_comp = c(2, 3), sim = c(0.7, 0.9), model = c("no_ints", "ints"), stringsAsFactors = FALSE)
SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
ROWDIR <- file.path(SP, "ex2_n01_rows"); dir.create(ROWDIR, showWarnings = FALSE)

fit_one <- function(row) {
  rp <- GRID$rep[row]; N_comp <- GRID$N_comp[row]; sim <- GRID$sim[row]; model <- GRID$model[row]
  wr <- function(df) { saveRDS(df, file.path(ROWDIR, sprintf("row_%02d.rds", row))); df }
  tryCatch({
    set.seed(run_seed + 1000 * rp + 10 * N_comp + round(100 * sim))
    gen <- sim_model_intercepts(N_unc = 5 - N_comp, N_comp_true = 2, N_comp_spur = N_comp, K_unc = 1, sim = sim, T_times = 30)
    test_ys <- gen$Y; groups <- gen$groups; N_ <- ncol(test_ys)
    perm <- anchor_order(test_ys, K_LAT); fit_ys <- test_ys[, perm]
    ints <- identical(model, "ints")
    os <- if (ints) apply(fit_ys, 2, sd) else apply(fit_ys, 2, function(x) sqrt(mean(x^2)))
    extra <- if (ints) list(include_ints = TRUE, int_scale = 3, int_loc = 4) else list(include_factor_means = TRUE)
    secs <- system.time(f <- do.call(sample_model, c(list(
      N_units = N_, T_times = nrow(fit_ys), K_latent = K_LAT, overall_scales = os,
      err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.05, data = fit_ys, autocor_a = 90, autocor_b = 10,
      nonstationary = FALSE, num_treated = NPT, fit_scales = TRUE, alpha_diag = 10, pathfinder_init = TRUE,
      type = "posterior", quiet = TRUE, ad = 0.8, iter = 600, iter_warm = 800, max_treedepth = 12, n_chains = NCH, seed = 42,
      return_draws = c("lp__", "tau", "delta", "sigma", "sigma_raw", "Lambda")), extra)))["elapsed"]
    d <- f$draws
    absz <- mean(abs(f$effect_means / f$effect_sds))
    cs <- unpermute_untreated(f$cor_sq, perm); ug <- groups[2:N_]
    cs_by <- sapply(c("true", "spurious", "uncorrelated"), function(g) if (any(ug == g)) mean(cs[ug == g]) else NA)
    mult <- mean(sapply(1:N_, function(n) mean(as.numeric(extract_variable_matrix(d, sprintf("sigma_raw[%d]", n))))))
    min_diag <- min(sapply(1:K_LAT, function(k) mean(as.numeric(extract_variable_matrix(d, sprintf("Lambda[%d,%d]", k, k))))))
    # mode signature: per-chain lp, on-mode subset rhat
    lpc <- colMeans(extract_variable_matrix(d, "lp__")); gap <- max(lpc) - lpc; onmode <- which(gap <= 5)
    idv <- c(grep("^tau", dimnames(d)$variable, value = TRUE), grep("^delta\\[", dimnames(d)$variable, value = TRUE),
      grep("^sigma\\[", dimnames(d)$variable, value = TRUE), grep("^Lambda\\[", dimnames(d)$variable, value = TRUE))
    rhat_all <- suppressWarnings(max(summarise_draws(subset_draws(d, variable = idv), "rhat")$rhat, na.rm = TRUE))
    rhat_on <- if (length(onmode) >= 2) suppressWarnings(max(summarise_draws(subset_draws(d, variable = idv, chain = onmode), "rhat")$rhat, na.rm = TRUE)) else NA
    sd_ <- f$sampler_diag
    used_model <- basename(get("ife_mod", envir = environment(sample_model))$stan_file())   # PROVE the model
    wr(data.frame(rep = rp, N_comp = N_comp, sim = sim, model = model, used = used_model, absz = round(absz, 3),
      cs_true = round(cs_by[1], 3), cs_spur = round(cs_by[2], 3), cs_unc = round(cs_by[3], 3), mult = round(mult, 2),
      min_diag = round(min_diag, 3), rhat_all = round(rhat_all, 2), rhat_on = round(rhat_on, 2), n_on = length(onmode),
      offmode = sd_$n_offmode, lp_gap = round(sd_$lp_gap_max, 1), div = sd_$n_div, secs = round(secs),
      lpc = paste(sprintf("%.0f", sort(lpc, decreasing = TRUE)), collapse = ","), row.names = NULL))
  }, error = function(e) wr(data.frame(rep = rp, N_comp = N_comp, sim = sim, model = model, used = NA, absz = NA, cs_true = NA,
    cs_spur = NA, cs_unc = NA, mult = NA, min_diag = NA, rhat_all = NA, rhat_on = NA, n_on = NA, offmode = NA,
    lp_gap = NA, div = NA, secs = NA, lpc = conditionMessage(e), row.names = NULL)))
}

cl <- makeCluster(4, outfile = "")
registerDoParallel(cl); invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, { library(cmdstanr); library(posterior) }))
# Export the MASTER's sample_model (+ its closure deps) and ife_mod, so %dopar% uses the intended model.
clusterExport(cl, c("fit_one", "sample_model", "ife_mod", "pathfinder_inits", "draw_to_init", "PF_PARAM_BASES",
  "ruv", "anchor_order", "unpermute_untreated", "sim_model_intercepts", "GRID", "K_LAT", "NPT", "run_seed", "NCH", "ROWDIR"))
res <- tryCatch(foreach(row = seq_len(nrow(GRID)), .combine = "rbind", .packages = c("cmdstanr", "posterior")) %dopar% fit_one(row), error = function(e) NULL)
tryCatch(stopCluster(cl), error = function(e) NULL)
if (is.null(res) || nrow(res) < nrow(GRID)) { rfs <- list.files(ROWDIR, pattern = "^row_.*\\.rds$", full.names = TRUE); res <- do.call(rbind, lapply(rfs, readRDS)) }
save(res, file = file.path(SP, "ex2_n01_validate.RData"))
res <- res[order(res$model, res$sim, res$N_comp, res$rep), ]
cat("\nex2 normal(0,1) sigma_raw | est sigma | 6 chains 800/600 | both models | science + mode diagnosis\n\n")
print(res[, setdiff(names(res), "lpc")], row.names = FALSE)
cat("\n=== by model ===\n")
for (m in c("no_ints", "ints")) { r <- res[res$model == m & is.finite(res$cs_spur), ]
  if (nrow(r)) cat(sprintf("%-8s: cs_spur=%.2f cs_true=%.2f | mult=%.2f min_diag=%.3f | rhat_all(max)=%.2f rhat_on(max)=%.2f offmode=%d div=%d secs=%.0f\n",
    m, mean(r$cs_spur), mean(r$cs_true), mean(r$mult), mean(r$min_diag), max(r$rhat_all), max(r$rhat_on, na.rm = TRUE), sum(r$offmode), sum(r$div), mean(r$secs))) }
cat("\nRead: (science) cs_spur separation back to ~0.15? mult~1? (modes) rhat_on clean while rhat_all high -> minor-mode contamination; both high -> systematic. min_diag near 0 -> zero-avoidance insufficient.\n")
