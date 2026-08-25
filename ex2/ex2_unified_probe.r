## Probe (NOT shipped): can a UNIFIED config (alpha_diag on + estimate sigma, matching ex1) work in
## ex2? And how much does alpha_diag bias delta here? 3 arms x both models x same 4 datasets as the
## est-sigma test (N_comp=3, sim=0.9). SHORT chains.
##   cur   = fixed sigma, alpha=0   (current ex2; delta~0 reference)
##   adiag = fixed sigma, alpha=10  (isolates alpha_diag's delta effect, sampling held clean)
##   unif  = est sigma,   alpha=10  (the unified candidate: does alpha tame the divergences? delta?)
## Read: does 'unif' sample clean (div~0, like ex1) AND keep delta near the 'cur' reference?

library(foreach); library(doParallel)
source("../sample_model.r"); source("../pathfinder_init.r")
for (e in parse("ex2_sim_study.r")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
    as.character(e[[2]]) %in% c("ruv", "sim_model_intercepts", "anchor_order", "unpermute_untreated"))
    eval(e, globalenv())
}

N_DATA <- 4L; K <- 3L; IT <- 400L; WM <- 400L; NCH <- 4L
set.seed(20260823)  # same datasets as ex2_est_sigma_test
datasets <- lapply(seq_len(N_DATA), function(i) {
  g <- sim_model_intercepts(N_unc = 2, N_comp_true = 2, N_comp_spur = 3, K_unc = 1, sim = 0.9, T_times = 30)
  list(Y = g$Y, groups = g$groups)
})

ARMS <- list(cur = list(alpha = 0, est = FALSE), adiag = list(alpha = 10, est = FALSE),
  unif = list(alpha = 10, est = TRUE))

model_cfg <- function(fit_ys, model) {
  if (model == "no_ints") {
    list(overall_scales = apply(fit_ys, 2, function(x) sqrt(mean(x^2))), include_factor_means = TRUE)
  } else {
    list(overall_scales = apply(fit_ys, 2, sd), include_ints = TRUE, int_scale = 3, int_loc = 4)
  }
}

fit_cell <- function(di, model, arm) {
  a <- ARMS[[arm]]
  tryCatch({
    ds <- datasets[[di]]; test_ys <- ds$Y; grp <- ds$groups[-1]
    perm <- anchor_order(test_ys, K); fit_ys <- test_ys[, perm]
    args <- c(list(N_units = ncol(fit_ys), T_times = nrow(fit_ys), K_latent = K,
      err_scale = 0, err_scale_mean = 0.1, err_scale_sd = 0.05, data = fit_ys,
      autocor_a = 90, autocor_b = 10, nonstationary = FALSE, num_treated = 5,
      alpha_diag = a$alpha, pathfinder_init = TRUE, type = "posterior", quiet = TRUE, ad = 0.8,
      iter = IT, iter_warm = WM, max_treedepth = 12, n_chains = NCH, seed = 42,
      parallel_chains = 1, fit_scales = if (a$est) TRUE else 0, return_draws = "sigma"),
      model_cfg(fit_ys, model))
    t <- system.time(f <- do.call(sample_model, args))
    cs <- unpermute_untreated(f$cor_sq, perm)
    sig_post <- tryCatch(mean(sapply(seq_len(ncol(fit_ys)), function(n)
      mean(as.numeric(posterior::extract_variable_matrix(f$draws, sprintf("sigma[%d]", n)))))),
      error = function(e) mean(args$overall_scales))
    sd_ <- f$sampler_diag
    data.frame(ds = di, model = model, arm = arm, secs = round(t[["elapsed"]]),
      sig = round(sig_post, 1), rhat = round(sd_$rhat_max, 3), div = sd_$n_div, tree = sd_$n_tree,
      ebfmi = round(sd_$ebfmi_min, 2), offmode = sd_$n_offmode,
      cor_true = round(mean(cs[grp == "true"]), 2), cor_spur = round(mean(cs[grp == "spurious"]), 2),
      delta_abs = round(mean(abs(f$effect_means)), 3), loc_pval = round(f$loc_cor_pval, 2), row.names = NULL)
  }, error = function(e) data.frame(ds = di, model = model, arm = arm, secs = NA, sig = NA, rhat = NA,
    div = NA, tree = NA, ebfmi = NA, offmode = NA, cor_true = NA, cor_spur = NA, delta_abs = NA,
    loc_pval = NA, row.names = NULL))
}

combos <- expand.grid(di = seq_len(N_DATA), model = c("no_ints", "ints"), arm = names(ARMS),
  stringsAsFactors = FALSE)
cl <- makeCluster(7, outfile = "")
registerDoParallel(cl); invisible(clusterCall(cl, setwd, getwd()))
invisible(clusterEvalQ(cl, { source("../sample_model.r"); source("../pathfinder_init.r") }))
clusterExport(cl, c("fit_cell", "model_cfg", "datasets", "anchor_order", "unpermute_untreated",
  "ARMS", "K", "IT", "WM", "NCH"))
res <- foreach(k = seq_len(nrow(combos)), .combine = "rbind",
  .packages = c("cmdstanr", "posterior")) %dopar% fit_cell(combos$di[k], combos$model[k], combos$arm[k])
stopCluster(cl)

SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
save(res, file = file.path(SP, "ex2_unified_probe.RData"))
res <- res[order(res$model, match(res$arm, names(ARMS)), res$ds), ]
cat(sprintf("\nex2 unified-config probe | cur(fix,a0) / adiag(fix,a10) / unif(est,a10) | K=%d, %dch x %d/%d | %d datasets\n\n", K, NCH, IT, WM, N_DATA))
print(res, row.names = FALSE)
for (m in c("no_ints", "ints")) for (arm in names(ARMS)) {
  r <- res[res$model == m & res$arm == arm, ]
  cat(sprintf("%-8s %-6s : div(sum)=%2d  mean rhat=%.3f  mean delta_abs=%.3f  mean cor_true/spur=%.2f/%.2f\n",
    m, arm, sum(r$div, na.rm = TRUE), mean(r$rhat, na.rm = TRUE), mean(r$delta_abs, na.rm = TRUE),
    mean(r$cor_true, na.rm = TRUE), mean(r$cor_spur, na.rm = TRUE)))
}
cat("\nRead: does 'unif' zero out the divergences 'est+a0' had? Is unif no_ints delta_abs close to 'cur'? (that's the alpha bias)\n")
