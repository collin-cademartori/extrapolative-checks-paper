## Fail-fast (NOT shipped): is ESTIMATING sigma viable in ex2 now that anchor+pathfinder are in?
## Both models (no_ints, ints) x {fixed sigma (current), estimate sigma} on a few DGP datasets
## (N_comp=3, sim=0.9 -> 8 units), sampling-fixes config (K=3, alpha_diag=0, anchor, pathfinder,
## ad=0.8), SHORT chains. Watch for the sigma-Lambda confound fingerprints (div / treedepth / rhat /
## offmode) and whether the science survives (cor_true>cor_spur, delta~0, S2 loc_cor_pval).

library(foreach); library(doParallel)
source("../sample_model.r"); source("../pathfinder_init.r")
# Reuse the exact DGP + anchor helpers from the shipped study.
for (e in parse("ex2_sim_study.r")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
    as.character(e[[2]]) %in% c("ruv", "sim_model_intercepts", "anchor_order", "unpermute_untreated"))
    eval(e, globalenv())
}

N_DATA <- 4L; K <- 3L; IT <- 400L; WM <- 400L; NCH <- 4L
set.seed(20260823)
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

fit_cell <- function(di, model, est) {
  tryCatch({
    ds <- datasets[[di]]; test_ys <- ds$Y; grp <- ds$groups[-1]
    perm <- anchor_order(test_ys, K); fit_ys <- test_ys[, perm]
    args <- c(list(N_units = ncol(fit_ys), T_times = nrow(fit_ys), K_latent = K,
      err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.05, data = fit_ys,
      autocor_a = 90, autocor_b = 10, nonstationary = FALSE, num_treated = 5,
      alpha_diag = 0, pathfinder_init = TRUE, type = "posterior", quiet = TRUE, ad = 0.8,
      iter = IT, iter_warm = WM, max_treedepth = 12, n_chains = NCH, seed = 42,
      parallel_chains = 1, fit_scales = if (est) TRUE else 0,
      return_draws = "sigma"), model_cfg(fit_ys, model))
    t <- system.time(f <- do.call(sample_model, args))
    cs <- unpermute_untreated(f$cor_sq, perm)
    sig_post <- tryCatch(mean(sapply(seq_len(ncol(fit_ys)), function(n)
      mean(as.numeric(posterior::extract_variable_matrix(f$draws, sprintf("sigma[%d]", n)))))),
      error = function(e) mean(args$overall_scales))
    sd_ <- f$sampler_diag
    data.frame(ds = di, model = model, sigma = if (est) "est" else "fixed", secs = round(t[["elapsed"]]),
      sig = round(sig_post, 1), rhat = round(sd_$rhat_max, 3), div = sd_$n_div, tree = sd_$n_tree,
      ebfmi = round(sd_$ebfmi_min, 2), ess_d1 = round(sd_$ess_delta1), offmode = sd_$n_offmode,
      cor_true = round(mean(cs[grp == "true"]), 2), cor_spur = round(mean(cs[grp == "spurious"]), 2),
      delta_abs = round(mean(abs(f$effect_means)), 3), loc_pval = round(f$loc_cor_pval, 2), row.names = NULL)
  }, error = function(e) data.frame(ds = di, model = model, sigma = if (est) "est" else "fixed",
    secs = NA, sig = NA, rhat = NA, div = NA, tree = NA, ebfmi = NA, ess_d1 = NA, offmode = NA,
    cor_true = NA, cor_spur = NA, delta_abs = NA, loc_pval = NA, row.names = NULL))
}

combos <- expand.grid(di = seq_len(N_DATA), model = c("no_ints", "ints"), est = c(FALSE, TRUE),
  stringsAsFactors = FALSE)
cl <- makeCluster(7, outfile = "")
registerDoParallel(cl); invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, { source("../sample_model.r"); source("../pathfinder_init.r") }))
clusterExport(cl, c("fit_cell", "model_cfg", "datasets", "anchor_order", "unpermute_untreated",
  "K", "IT", "WM", "NCH"))
res <- foreach(k = seq_len(nrow(combos)), .combine = "rbind",
  .packages = c("cmdstanr", "posterior")) %dopar% fit_cell(combos$di[k], combos$model[k], combos$est[k])
stopCluster(cl)

SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
save(res, file = file.path(SP, "ex2_est_sigma_test.RData"))
res <- res[order(res$model, res$sigma, res$ds), ]
cat(sprintf("\nex2 estimate-sigma viability | both models x {fixed, est} | K=%d, %dch x %d/%d | %d datasets\n", K, NCH, IT, WM, N_DATA))
cat("science should survive: cor_true > cor_spur, delta_abs ~ 0. sampling should stay clean: div=0, low tree, rhat<=~1.02, offmode=0\n\n")
print(res, row.names = FALSE)
cat("\nRead: does the 'est' block match 'fixed' on sampling (no new div/tree/rhat/offmode) AND keep the cor_true>cor_spur contrast?\n")
