## ex2 SCIENCE check for the Cartesian model (ife_named_cartesian.stan): does it reproduce the two ex2
## signatures now that sampling is clean (0 div, 0 treedepth, 2x faster)?
##  (1) INTS worse than NO_INTS at the (true=0) effect  (2) INTS inflates cor_sq for the SPURIOUS units.
## Cartesian = original unconstrained loadings + error scale anchored to the ESTIMATED signal scale
## scale_est[n] = sigma_data[n]*||Lambda[n,:]||, which is algebraically the unit-norm + est-sigma model.
## fit_scales = 0 (per-unit scale is estimated through ||Lambda||, NOT a separate sigma parameter).
## ad = 0.8 for iteration speed. Master sets ife_mod and exports it (sample_model resolves its MASTER
## closure -- a clusterEvalQ override alone is silently ignored). Crash-safe per-fit rows.
library(foreach); library(doParallel); library(posterior)
source("../sample_model.r"); source("../pathfinder_init.r")
ife_mod <- cmdstan_model("../ife_named_cartesian.stan")
stopifnot(basename(ife_mod$stan_file()) == "ife_named_cartesian.stan", "Lambda" %in% PF_PARAM_BASES)

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
  list(Y = t(lat + rnorm(nrow(lat) * ncol(lat), sd = 0.1 * max(apply(lat, 1, sd)))),
       groups = c("treated", rep("true", N_comp_true), rep("spurious", N_comp_spur), rep("uncorrelated", N_unc)))
}

K_LAT <- 3L; NPT <- 5L; run_seed <- 88213
GRID <- expand.grid(rep = 1:2, N_comp = c(2, 3), sim = c(0.7, 0.9), model = c("no_ints", "ints"), stringsAsFactors = FALSE)
SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
ROWDIR <- file.path(SP, "ex2_cart_rows"); dir.create(ROWDIR, showWarnings = FALSE)

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
      nonstationary = FALSE, num_treated = NPT, fit_scales = 0, alpha_diag = 10, pathfinder_init = TRUE,
      type = "posterior", quiet = TRUE, ad = 0.8, iter = 500, iter_warm = 500, max_treedepth = 12,
      n_chains = 4, seed = 42, return_draws = c("tau", "scale_est", "Lambda")), extra)))["elapsed"]
    d <- f$draws
    absz <- mean(abs(f$effect_means / f$effect_sds)); eabs <- mean(abs(f$effect_means))
    cs <- unpermute_untreated(f$cor_sq, perm); ug <- groups[2:N_]
    cs_by <- sapply(c("true", "spurious", "uncorrelated"), function(g) if (any(ug == g)) mean(cs[ug == g]) else NA)
    tau <- median(as.numeric(extract_variable_matrix(d, "tau")))
    sc <- mean(sapply(1:N_, function(n) mean(as.numeric(extract_variable_matrix(d, sprintf("scale_est[%d]", n))))))
    mult <- sc / mean(os)                         # estimated scale relative to the data scale (~1 wanted)
    yp <- f$y_pred; T_ <- nrow(fit_ys); inc <- 0
    for (n in 1:N_) for (tt in 1:T_) { b <- quantile(yp[, tt, n], c(.005, .995)); inc <- inc + (fit_ys[tt, n] >= b[1] && fit_ys[tt, n] <= b[2]) }
    sd_ <- f$sampler_diag
    wr(data.frame(rep = rp, N_comp = N_comp, sim = sim, model = model, absz = round(absz, 3), eabs = round(eabs, 3),
      cs_true = round(cs_by[1], 3), cs_spur = round(cs_by[2], 3), cs_unc = round(cs_by[3], 3),
      tau = round(tau, 3), scale_est = round(sc, 2), mult = round(mult, 2), cover = round(inc / (T_ * N_), 3),
      div = sd_$n_div, tree = sd_$n_tree, rhat = round(sd_$rhat_max, 2), secs = round(secs), row.names = NULL))
  }, error = function(e) wr(data.frame(rep = rp, N_comp = N_comp, sim = sim, model = model, absz = NA, eabs = NA,
    cs_true = NA, cs_spur = NA, cs_unc = NA, tau = NA, scale_est = NA, mult = NA, cover = NA, div = NA, tree = NA,
    rhat = NA, secs = NA, row.names = NULL)))
}

cl <- makeCluster(4, outfile = "")
registerDoParallel(cl); invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, { library(cmdstanr); library(posterior) }))
clusterExport(cl, c("fit_one", "sample_model", "ife_mod", "pathfinder_inits", "draw_to_init", "PF_PARAM_BASES",
  "ruv", "anchor_order", "unpermute_untreated", "sim_model_intercepts", "GRID", "K_LAT", "NPT", "run_seed", "ROWDIR"))
res <- tryCatch(foreach(row = seq_len(nrow(GRID)), .combine = "rbind", .packages = c("cmdstanr","posterior")) %dopar% fit_one(row), error = function(e) NULL)
tryCatch(stopCluster(cl), error = function(e) NULL)
if (is.null(res) || nrow(res) < nrow(GRID)) { rfs <- list.files(ROWDIR, pattern="^row_.*\\.rds$", full.names=TRUE); res <- do.call(rbind, lapply(rfs, readRDS)) }
save(res, file = file.path(SP, "ex2_cartesian_science.RData"))
res <- res[order(res$model, res$sim, res$N_comp, res$rep), ]
cat("\nex2 CARTESIAN model | science check | ad=0.8, 4 chains | both models\n")
cat("(targets: ints cs_spur > no_ints cs_spur with a clear gap; ints worse on absz/eabs; tau ~0.1 in prior; div ~0)\n\n")
print(res, row.names = FALSE)
cat("\n=== by model ===\n")
for (m in c("no_ints","ints")) { r <- res[res$model==m & is.finite(res$cs_spur), ]
  if (nrow(r)) cat(sprintf("%-8s n=%d: cs_spur=%.2f cs_true=%.2f cs_unc=%.2f | absz=%.2f eabs=%.2f | tau=%.3f mult=%.2f | div=%d tree=%d rhat=%.2f secs=%.0f\n",
    m, nrow(r), mean(r$cs_spur), mean(r$cs_true), mean(r$cs_unc), mean(r$absz), mean(r$eabs),
    mean(r$tau), mean(r$mult), sum(r$div), sum(r$tree), max(r$rhat), mean(r$secs))) }
cat("\nRead: gap (ints cs_spur - no_ints cs_spur) vs polar-unitnorm 0.27 and original-fixed 0.15; and is tau in its prior with div~0?\n")
