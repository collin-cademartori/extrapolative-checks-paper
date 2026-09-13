## Pathfinder initialization for the HMC fits. The posterior has minor modes with negligibly small
## posterior density. Initializing chains inside the dominant mode keeps them
## out. pathfinder_inits() runs multi-path Pathfinder, keeps the draws with log density within
## lp_gap of the maximum, and initializes each chain from one of those draws.

# Parameters of ife_named.stan for pathfinder to initialize.
PF_PARAM_BASES <- c("tau_param", "omega_sq_param", "Phi_innovations",
  "Phi_means_param", "rho", "Lambda", "gamma_raw", "delta_raw")

# Convert a named draw (like "Lambda[2,1]" or "rho[3]") into a nested Stan init list.
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

  present <- tryCatch(
    intersect(PF_PARAM_BASES, sub("\\[.*", "", posterior::variables(pf$draws()))),
    error = function(e) { warning("pathfinder draws unreadable: ", conditionMessage(e)); NULL })
  if (is.null(present)) return(NULL)
  dd <- tryCatch(as.data.frame(posterior::as_draws_df(pf$draws(c("lp__", present)))),
    error = function(e) { warning("pathfinder draw extraction failed: ", conditionMessage(e)); NULL })
  if (is.null(dd) || nrow(dd) == 0) return(NULL)
  dd <- dd[is.finite(dd$lp__), , drop = FALSE]
  if (nrow(dd) == 0) return(NULL)

  # Dominant mode defined as draws within lp_gap of the maximum log density.
  dom <- which(dd$lp__ >= max(dd$lp__) - lp_gap)
  if (length(dom) == 0) return(NULL)

  pcols <- setdiff(names(dd), c("lp__", "lp_approx__", ".chain", ".iteration", ".draw"))
  # One distinct draw per chain (with replacement only if there are fewer draws than chains).
  if (!is.null(seed)) {
    have_state <- exists(".Random.seed", envir = globalenv())
    old_state <- if (have_state) get(".Random.seed", envir = globalenv()) else NULL
    set.seed(seed)
    pick <- sample(dom, n_chains, replace = length(dom) < n_chains)
    if (have_state) assign(".Random.seed", old_state, envir = globalenv())
  } else {
    pick <- sample(dom, n_chains, replace = length(dom) < n_chains)
  }
  inits <- lapply(pick, function(i) draw_to_init(unlist(dd[i, pcols])))
  attr(inits, "n_dominant") <- length(dom)
  attr(inits, "n_pf_draws") <- nrow(dd)
  try(unlink(pf$output_files(), force = TRUE), silent = TRUE)
  inits
}
