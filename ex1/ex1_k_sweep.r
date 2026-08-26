## Does the loadings' weak identification (bimodal Lambda marginals -> intermittent high rhat) get
## better or worse with more factors? Fit K_latent in {3, 4, 5} while holding the DGP at K = 4, so only
## the MODEL's factor count varies -- in ex1_sim_study the two are coupled, which would confound data
## with model. Hypothesis: extra factors beyond what the data support are weakly determined and free to
## re-label, so K = 5 should be worse and K = 3 better.
## Reports rhat split three ways -- ALL params, LAMBDA only, ESTIMANDS only (delta/tau/rho/scale_est) --
## plus cor_sq, the sign/label-invariant functional actually reported from the loadings.
library(posterior)
source("../sample_model.r"); source("../pathfinder_init.r")
ife_mod <- cmdstan_model("../ife_named_cartesian.stan")
rms <- function(x) sqrt(mean(x^2))

GEN_K <- 4L; REPS <- 30L; NPT <- 5L
set.seed(40318)
study_units <- sample.int(2 * REPS, size = REPS, replace = FALSE)
pp_seed <- sample.int(.Machine$integer.max, 1)
td <- sample_model(overall_scales = rep(1, 8), err_scale = 2, alpha_diag = 10, autocor_a = 8,
  autocor_b = 2, nonstationary = TRUE, num_treated = 0, type = "prior_pred", K_latent = GEN_K,
  iter = 2 * REPS, seed = pp_seed)

blk <- function(d, pat) {
  v <- grep(pat, dimnames(d)$variable, value = TRUE)
  if (!length(v)) return(c(rhat = NA, ess = NA))
  ss <- summarise_draws(subset_draws(d, variable = v), "rhat", "ess_bulk")
  c(rhat = suppressWarnings(max(ss$rhat, na.rm = TRUE)), ess = suppressWarnings(min(ss$ess_bulk, na.rm = TRUE)))
}

cat("\nex1 stat_strong | K_latent sweep {3,4,5} | DGP fixed at K=4 | 2x RMS, tau~N(0.05,0.05)+, ad=0.8\n")
cat("rhat split: ALL params | LAMBDA only | ESTIMANDS (delta,tau,rho,scale_est) | cor_sq (sign-invariant)\n\n")
out <- list()
for (ds in c(42, 27)) {
  ys <- td$ys[ds, , ]; N_ <- ncol(ys); T_ <- nrow(ys)
  for (K in c(3L, 4L, 5L)) for (sd_i in c(500, 7, 1234)) {
    r <- tryCatch({
      f <- sample_model(N_units = N_, T_times = T_, K_latent = K,
        overall_scales = 2 * apply(ys, 2, rms),
        err_scale = 0, err_scale_mean = 0.05, err_scale_sd = 0.05, data = ys,
        autocor_a = 97, autocor_b = 3, nonstationary = FALSE, num_treated = NPT,
        fit_scales = FALSE, alpha_diag = 10, pathfinder_init = TRUE, type = "posterior", quiet = TRUE,
        ad = 0.8, iter = 1000, iter_warm = 500, n_chains = 4, seed = sd_i, parallel_chains = 4,
        return_draws = c("lp__", "tau", "rho", "Lambda", "delta", "scale_est", "cor_sq"))
      d <- f$draws
      lam <- blk(d, "^Lambda\\[")
      est <- blk(d, "^(delta|rho|scale_est)\\[|^tau$")
      cs  <- blk(d, "^cor_sq\\[")
      # how bimodal are the loadings: median ESS gain from taking |.|
      lv <- grep("^Lambda\\[", dimnames(d)$variable, value = TRUE)
      gain <- median(sapply(lv, function(v) { m <- extract_variable_matrix(d, v)
        e <- ess_bulk(m); if (is.na(e) || e < 1) NA else ess_bulk(abs(m)) / e }), na.rm = TRUE)
      data.frame(ds = ds, K = K, seed = sd_i, rhat_all = round(f$sampler_diag$rhat_max, 3),
        rhat_lam = round(lam["rhat"], 3), ess_lam = round(lam["ess"]),
        rhat_est = round(est["rhat"], 4), ess_est = round(est["ess"]),
        rhat_cor = round(cs["rhat"], 3), absgain = round(gain, 2),
        div = f$sampler_diag$n_div, offmode = f$sampler_diag$n_offmode, row.names = NULL)
    }, error = function(e) { cat(sprintf("  !! ds%d K%d seed%d FAILED\n", ds, K, sd_i)); NULL })
    if (is.null(r)) next
    out[[length(out) + 1]] <- r
    cat(sprintf("ds%-3d K=%d seed%-6d rhat_all=%-6.3f | Lambda rhat=%-6.3f ess=%-5.0f |.|gain=%.1fx | estimands rhat=%-7.4f ess=%-5.0f | cor_sq rhat=%-6.3f | div=%d\n",
      ds, K, sd_i, r$rhat_all, r$rhat_lam, r$ess_lam, r$absgain, r$rhat_est, r$ess_est, r$rhat_cor, r$div))
  }
}
res <- do.call(rbind, out)
saveRDS(res, "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad/ex1_k_sweep.rds")
cat("\n=== by K (DGP truth is K=4) ===\n")
for (K in c(3L, 4L, 5L)) { s <- res[res$K == K, ]; if (!nrow(s)) next
  cat(sprintf("K=%d (n=%d): rhat_all med=%.3f max=%.3f | Lambda rhat max=%.3f ess min=%.0f | |.|gain med=%.1fx | ESTIMANDS rhat max=%.4f | cor_sq rhat max=%.3f\n",
    K, nrow(s), median(s$rhat_all), max(s$rhat_all), max(s$rhat_lam), min(s$ess_lam),
    median(s$absgain), max(s$rhat_est), max(s$rhat_cor))) }
cat("\nRead: if rhat/ESS degrade from K=3 to K=5, surplus factors are the driver and reducing K is a real\n")
cat("fix. If cor_sq and the estimands stay clean at every K, the instability is confined to the\n")
cat("unidentified loading configuration regardless of K.\n")
