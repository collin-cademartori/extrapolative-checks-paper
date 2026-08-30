## Pathfinder-based initialization for the HMC fits (branch dgp-diag). The interactive-fixed-effects
## posterior has tiny but real collapsed-loading minor modes (a diagonal Lambda[k,k] -> 0). Seeding
## HMC from a coherent point INSIDE the dominant basin keeps chains out of those modes; a short HMC
## pilot works but is slow, so we use multi-path Pathfinder (a fast L-BFGS variational mode-finder)
## instead -- its posterior-approximation error is irrelevant because we only use its draws as HMC
## starting points, not as inference.
##
## pathfinder_inits(): run multi-path Pathfinder, keep only draws in the dominant lp mode (within
## lp_gap of the top log density -- a single, clean criterion), and return one full-parameter init
## per chain, drawn from distinct dominant draws (natural, constraint-valid dispersion -- no manual
## perturbation). Returns NULL on failure so the caller can fall back to Stan's default init.

# The `parameters` block of ife_named.stan (constrained draws of these seed every free parameter;
# zero-length ones -- e.g. sigma_raw when sigma is fixed -- simply won't appear in the draws).
PF_PARAM_BASES <- c("sigma_raw", "tau_param", "omega_sq_param", "Phi_innovations",
  "Phi_means_param", "rho", "Lambda", "gamma_raw", "delta_raw")

# Rebuild one nested Stan init list from a single named draw (names like "Lambda[2,1]", "rho[3]",
# "tau_param[1]"): scalars -> length-1, "v[i]" -> vector, "m[i,j]" -> matrix. Structural zeros of a
# cholesky_factor_cov are filled from the array's zero initialization, so absent upper entries stay 0.
draw_to_init <- function(vals) {
  base <- sub("\\[.*", "", names(vals))
  out <- list()
  for (b in unique(base)) {
    idx <- which(base == b)
    subs <- regmatches(names(vals)[idx], regexpr("\\[.*\\]", names(vals)[idx]))
    if (length(subs) == 0) {                        # scalar (no brackets)
      out[[b]] <- unname(vals[idx])
    } else {
      ij <- do.call(rbind, lapply(strsplit(gsub("\\[|\\]", "", subs), ","), as.integer))
      d <- apply(ij, 2, max)
      arr <- array(0, dim = d)
      for (r in seq_along(idx)) arr[matrix(ij[r, , drop = FALSE], nrow = 1)] <- vals[idx[r]]
      out[[b]] <- if (length(d) == 1) as.vector(arr) else arr
    }
  }
  out
}

pathfinder_inits <- function(mod, stat_data, n_chains, num_paths = NULL, draws = NULL,
                             lp_gap = 10, seed = NULL, quiet = TRUE, output_dir = NULL) {
  num_paths <- if (is.null(num_paths)) max(8L, 2L * n_chains) else num_paths
  draws <- if (is.null(draws)) max(400L, 40L * n_chains) else draws

  pf <- tryCatch(
    mod$pathfinder(data = stat_data, num_paths = num_paths, draws = draws, seed = seed,
      refresh = 0, show_messages = !quiet, show_exceptions = !quiet, output_dir = output_dir),
    error = function(e) { warning("pathfinder failed: ", conditionMessage(e)); NULL })
  if (is.null(pf)) return(NULL)

  # pf$draws() must be guarded too, not just $pathfinder(). Pathfinder can return an object whose
  # CSV was never written -- when every path fails, as happens if the initialization lands somewhere
  # the likelihood is degenerate -- and then reading it throws from read_cmdstan_csv rather than
  # returning an error object. Unguarded, that propagated out of sample_model and killed the fit
  # instead of falling back to Stan's default init, which is the whole point of returning NULL here.
  present <- tryCatch(
    intersect(PF_PARAM_BASES, sub("\\[.*", "", posterior::variables(pf$draws()))),
    error = function(e) { warning("pathfinder draws unreadable: ", conditionMessage(e)); NULL })
  if (is.null(present)) return(NULL)
  dd <- tryCatch(as.data.frame(posterior::as_draws_df(pf$draws(c("lp__", present)))),
    error = function(e) { warning("pathfinder draw extraction failed: ", conditionMessage(e)); NULL })
  if (is.null(dd) || nrow(dd) == 0) return(NULL)
  dd <- dd[is.finite(dd$lp__), , drop = FALSE]
  if (nrow(dd) == 0) return(NULL)

  # dominant lp mode: draws within lp_gap of the top log density (a single, clean criterion)
  dom <- which(dd$lp__ >= max(dd$lp__) - lp_gap)
  if (length(dom) == 0) return(NULL)

  pcols <- setdiff(names(dd), c("lp__", "lp_approx__", ".chain", ".iteration", ".draw"))
  if (!is.null(seed)) set.seed(seed)
  # distinct dominant draws when there are enough; fall back to with-replacement only if the
  # dominant set is smaller than the chain count (else sample() errors).
  pick <- sample(dom, n_chains, replace = length(dom) < n_chains)
  inits <- lapply(pick, function(i) draw_to_init(unlist(dd[i, pcols])))
  attr(inits, "n_dominant") <- length(dom)
  attr(inits, "n_pf_draws") <- nrow(dd)
  # Pathfinder writes its own CSV to the temp directory (~4 MB here) and, like $sample(), cmdstanr
  # only removes it on garbage collection -- which does not reliably happen inside a long study loop.
  # The draws are already extracted into `dd`, so drop the file now. See sample_model.r for the
  # full rationale.
  try(unlink(pf$output_files(), force = TRUE), silent = TRUE)
  inits
}
