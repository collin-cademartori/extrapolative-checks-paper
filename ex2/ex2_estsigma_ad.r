## Probe (NOT shipped): are ex2's estimate-sigma divergences just a low-adapt_delta artifact?
## ex2 fits run ad=0.8; ex1 runs ad=0.98. Test est-sigma + alpha=10 at ad in {0.95, 0.99} on the same
## 4 datasets, both models. Compare div to the ad=0.8 result (no_ints sum=100, ints sum=26).
## If divergences clear at ad~0.98, then "estimate sigma + alpha + ad=0.98" is a clean UNIFIED config.

library(foreach); library(doParallel)
source("../sample_model.r"); source("../pathfinder_init.r")
for (e in parse("ex2_sim_study.r")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
    as.character(e[[2]]) %in% c("ruv", "sim_model_intercepts", "anchor_order", "unpermute_untreated"))
    eval(e, globalenv())
}

N_DATA <- 4L; K <- 3L; IT <- 400L; WM <- 500L; NCH <- 4L
set.seed(20260823)  # same datasets as the other ex2 probes
datasets <- lapply(seq_len(N_DATA), function(i) {
  g <- sim_model_intercepts(N_unc = 2, N_comp_true = 2, N_comp_spur = 3, K_unc = 1, sim = 0.9, T_times = 30)
  list(Y = g$Y, groups = g$groups)
})

model_cfg <- function(fit_ys, model) {
  if (model == "no_ints") {
    list(overall_scales = apply(fit_ys, 2, function(x) sqrt(mean(x^2))), include_factor_means = TRUE)
  } else {
    list(overall_scales = apply(fit_ys, 2, sd), include_ints = TRUE, int_scale = 3, int_loc = 4)
  }
}

fit_cell <- function(di, model, ad) {
  tryCatch({
    ds <- datasets[[di]]; test_ys <- ds$Y; grp <- ds$groups[-1]
    perm <- anchor_order(test_ys, K); fit_ys <- test_ys[, perm]
    args <- c(list(N_units = ncol(fit_ys), T_times = nrow(fit_ys), K_latent = K,
      err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.05, data = fit_ys,
      autocor_a = 90, autocor_b = 10, nonstationary = FALSE, num_treated = 5,
      alpha_diag = 10, pathfinder_init = TRUE, type = "posterior", quiet = TRUE, ad = ad,
      iter = IT, iter_warm = WM, max_treedepth = 12, n_chains = NCH, seed = 42,
      parallel_chains = 1, fit_scales = TRUE, return_draws = "sigma"), model_cfg(fit_ys, model))
    t <- system.time(f <- do.call(sample_model, args))
    cs <- unpermute_untreated(f$cor_sq, perm)
    sd_ <- f$sampler_diag
    data.frame(ds = di, model = model, ad = ad, secs = round(t[["elapsed"]]),
      rhat = round(sd_$rhat_max, 3), div = sd_$n_div, tree = sd_$n_tree,
      ebfmi = round(sd_$ebfmi_min, 2), offmode = sd_$n_offmode,
      cor_true = round(mean(cs[grp == "true"]), 2), cor_spur = round(mean(cs[grp == "spurious"]), 2),
      delta_abs = round(mean(abs(f$effect_means)), 3), row.names = NULL)
  }, error = function(e) data.frame(ds = di, model = model, ad = ad, secs = NA, rhat = NA, div = NA,
    tree = NA, ebfmi = NA, offmode = NA, cor_true = NA, cor_spur = NA, delta_abs = NA, row.names = NULL))
}

combos <- expand.grid(di = seq_len(N_DATA), model = c("no_ints", "ints"), ad = c(0.95, 0.99))
cl <- makeCluster(7, outfile = "")
registerDoParallel(cl); invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, { source("../sample_model.r"); source("../pathfinder_init.r") }))
clusterExport(cl, c("fit_cell", "model_cfg", "datasets", "anchor_order", "unpermute_untreated",
  "K", "IT", "WM", "NCH"))
res <- foreach(k = seq_len(nrow(combos)), .combine = "rbind",
  .packages = c("cmdstanr", "posterior")) %dopar% fit_cell(combos$di[k], combos$model[k], combos$ad[k])
stopCluster(cl)

SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
save(res, file = file.path(SP, "ex2_estsigma_ad.RData"))
res <- res[order(res$model, res$ad, res$ds), ]
cat(sprintf("\nex2 est-sigma + alpha=10, adapt_delta sweep | K=%d, %dch x %d(warm %d) | %d datasets\n", K, NCH, IT, WM, N_DATA))
cat("(ref: ad=0.8 est+a10 gave div sums no_ints=100, ints=26)\n\n")
print(res, row.names = FALSE)
for (m in c("no_ints", "ints")) for (a in c(0.95, 0.99)) {
  r <- res[res$model == m & res$ad == a, ]
  cat(sprintf("%-8s ad=%.2f : div(sum)=%2d  max tree=%d  mean rhat=%.3f  mean secs=%.0f  delta_abs=%.3f\n",
    m, a, sum(r$div, na.rm = TRUE), max(r$tree, na.rm = TRUE), mean(r$rhat, na.rm = TRUE),
    mean(r$secs, na.rm = TRUE), mean(r$delta_abs, na.rm = TRUE)))
}
cat("\nRead: does raising ad drive div -> ~0 (=> ex2 est-sigma was just an adapt_delta artifact; unified config viable)?\n")
