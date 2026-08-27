## Verify the memory mechanism blamed for the OOM kill, BEFORE committing to another long run.
##
## THE HYPOTHESIS
## The second long ex1 run was killed by the OS at ~70 of 200 reps, exactly as the first round-3
## fits were starting. The claim is that sample_model() called model_sample$draws() with NO
## arguments to extract delta (5 columns) and cor_sq (7), which materializes all ~1270 of the
## model's columns; cmdstanr then caches that array for the life of the fit, so it stays resident.
## Measured on a laptop at 10.1 KB per stored draw, that is 61 MB at the base iter = 2000 x 3 chains
## but 607 MB at iter = 20000 -- in every worker that reaches round 3, concurrently.
##
## WHAT WOULD FALSIFY IT
##   * If "bare" and "named" come out roughly equal, the bare $draws() call is not the driver and
##     the OOM has another cause -- in which case the round-3 reduction to 12000 was the wrong
##     remedy and should be reconsidered.
##   * If the bare/named gap does not grow with iterations, the caching story is wrong.
##
## WHAT IT MEASURES
## The primary metric is R's own heap high-water mark, via gc()'s "max used" over the extraction
## window. It is deterministic and behaves identically on Linux and macOS. Process RSS is reported
## alongside but should NOT be judged by on macOS: `ps` gives whatever is resident at the moment of
## the call, so whether a garbage collection happened to run first dominates it. Two macOS trial
## runs of an earlier version of this script disagreed by 380 MB at the same rung for that reason,
## once even showing the bare arm as CHEAPER. On Linux the RSS column reads VmHWM, a genuine
## process peak, and is what the worker-count projection should be based on.
##
## Each configuration runs in its own fresh Rscript subprocess, so one cannot contaminate the next.
##
## LAPTOP REFERENCE (macOS, heap metric, for comparison against what you get):
##   iter=1000   bare 254 MB  named 197 MB   saved  56 MB (1.3x)
##   iter=2500   bare 457 MB  named 271 MB   saved 187 MB (1.7x)
## i.e. the saving grows at roughly 87 MB per 1000 sampling iterations, which extrapolates to about
## 1 GB per fit saved at iter = 12000. The post-fix floor there was ~148 MB + 50 KB per iteration.
##
## The fit is built by calling ife_mod$sample() directly with ex1's stat_weak configuration, rather
## than through sample_model(), so the result does not depend on which version of sample_model.r is
## checked out -- both arms are measured against the same fit code, differing only in the one line
## under test.
##
## USAGE
##   Rscript ex1_memory_check.r                 # rungs 2000, 8000, 12000
##   Rscript ex1_memory_check.r 2000,8000,20000 # also reproduce the rung that died
##   Rscript ex1_memory_check.r 2000,8000 24     # ... and project for 24 workers
## Runtime is roughly 20-40 min for the default rungs (six fits, run sequentially and deliberately
## NOT in parallel, so the measurements stay clean and the test cannot itself exhaust memory).

suppressPackageStartupMessages({library(cmdstanr); library(posterior)})

peak_rss_mb <- function() {
  if (file.exists("/proc/self/status")) {
    v <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
    if (length(v)) return(as.numeric(gsub("[^0-9]", "", v)) / 1024)
  }
  as.numeric(system(sprintf("ps -o rss= -p %d", Sys.getpid()), intern = TRUE)) / 1024
}
is_true_peak <- file.exists("/proc/self/status")

total_ram_gb <- function() {
  if (file.exists("/proc/meminfo")) {
    v <- grep("^MemTotal:", readLines("/proc/meminfo"), value = TRUE)
    if (length(v)) return(as.numeric(gsub("[^0-9]", "", v)) / 1024 / 1024)
  }
  out <- suppressWarnings(try(system("sysctl -n hw.memsize", intern = TRUE), silent = TRUE))
  if (!inherits(out, "try-error") && length(out)) return(as.numeric(out) / 1024^3)
  NA_real_
}

rms <- function(x) sqrt(mean(x^2))

# ---- the fit under test: ex1's stat_weak configuration ----------------------------------------
build_data <- function(unit = 90L) {
  # Same prior-predictive DGP as the study, so the fit sees realistic data.
  set.seed(40318)
  study_units <- sample.int(400L, size = 200L, replace = FALSE)
  pp_seed <- sample.int(.Machine$integer.max, 1)
  ife_mod <- cmdstan_model("../ife_named.stan")
  pp <- ife_mod$sample(
    data = list(M_units = 8, T_times = 20, K_latent = 4L, Y = matrix(0, 20, 8),
      a_rho = 8, b_rho = 2, tau_val = 2, m_tau = 0, s_tau = 0, sigma_data = rep(1, 8),
      fit_overall_scales = 0, nonstationary = 1, unit_intercepts = 0, factor_means = 0,
      sample_posterior = 0, num_treated = 0, gamma_scale = 1, gamma_loc = 0, alpha_diag = 10),
    chains = 1, iter_warmup = 200, iter_sampling = 400, refresh = 0,
    show_messages = FALSE, show_exceptions = FALSE, seed = pp_seed)
  ys <- extract_variable_array(pp$draws("Y_prior"), "Y_prior")[unit, 1, , ]
  try(unlink(pp$output_files(), force = TRUE), silent = TRUE)
  list(mod = ife_mod, ys = ys)
}

run_one <- function(mode, iter) {
  bd <- build_data()
  ys <- bd$ys
  stat_data <- list(
    M_units = ncol(ys), T_times = nrow(ys), K_latent = 5L, Y = ys,
    a_rho = 97, b_rho = 3, tau_val = 0, m_tau = 0.1, s_tau = 0.1,
    sigma_data = 2 * apply(ys, 2, rms), fit_overall_scales = 0,
    nonstationary = 0, unit_intercepts = 0, factor_means = 0,
    sample_posterior = 1, num_treated = 5,
    gamma_scale = 1, gamma_loc = 0, alpha_diag = 20)

  # Pathfinder init is included because the study uses it and it contributes to the peak; leaving
  # it out would make the projected safe-worker-count optimistic.
  source("../pathfinder_init.r")
  pfi <- tryCatch(pathfinder_inits(bd$mod, stat_data, 3, seed = 4242, quiet = TRUE),
    error = function(e) NULL)

  f <- bd$mod$sample(data = stat_data, chains = 3, parallel_chains = 3,
    iter_warmup = 500, iter_sampling = iter, adapt_delta = 0.8,
    init = if (is.null(pfi)) 2 else pfi,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE, seed = 4242)

  # --- measurement window opens -----------------------------------------------------------------
  # R's own heap high-water mark is the primary metric, NOT process RSS. gc()'s "max used" is a
  # true high-water mark since the last reset, is deterministic, and works identically on Linux and
  # macOS. Process RSS sampled after the fact is not: it reports whatever is resident at the moment
  # of the call, so whether a garbage collection happened to run first dominates the reading. Two
  # macOS trial runs of this script disagreed by 380 MB at the same rung for exactly that reason.
  invisible(gc(reset = TRUE, full = TRUE))

  # Everything sample_model() does downstream, in the same order. The ONLY difference between the
  # two arms is the delta / cor_sq extraction.
  sample_index <- sample.int(iter, size = iter)
  y_means <- extract_variable_array(f$draws("Y_latent"), "Y_latent")[sample_index, 1, , ]
  y_pred  <- extract_variable_array(f$draws("Y_pred"), "Y_pred")[sample_index, 1, , ]

  if (mode == "bare") {
    effects <- extract_variable_array(f$draws(), "delta")[, 1, ]          # the suspect line
    cor_sq  <- extract_variable_array(f$draws(), "cor_sq")[, 1, ]         # the suspect line
  } else {
    effects <- extract_variable_array(f$draws("delta"), "delta")[, 1, ]   # the fix
    cor_sq  <- extract_variable_array(f$draws("cor_sq"), "cor_sq")[, 1, ] # the fix
  }
  tau <- as.numeric(f$draws("tau"))
  mad <- mean(as.numeric(f$draws("mean_abs_diffs")))
  sv <- f$metadata()$stan_variables
  vars <- intersect(c("Lambda", "Phi_innovations", "rho", "delta_raw", "tau_param"), sv)
  ps <- f$summary(vars, "rhat", "ess_bulk")
  rM <- suppressWarnings(max(f$summary("Lambda_Phi", "rhat")$rhat, na.rm = TRUE))

  g <- gc(full = TRUE)
  heap_peak <- sum(g[, "max used"] * c(56, 8)) / 1e6   # R heap high-water since the reset, MB
  # --- measurement window closes ----------------------------------------------------------------

  n_draws <- iter * 3
  cat(sprintf("RESULT %s %d %.1f %.1f %.3f %.1f\n", mode, iter, heap_peak,
    1024 * heap_peak / n_draws, rM, peak_rss_mb()))
  try(unlink(f$output_files(), force = TRUE), silent = TRUE)
  invisible(NULL)
}

# ---- worker branch --------------------------------------------------------------------------
# Resolve this file's own path so the driver can re-invoke it as a worker. --file= is what Rscript
# passes; the literal is only a fallback for interactive sourcing.
THIS_SCRIPT <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (!length(THIS_SCRIPT)) THIS_SCRIPT <- "ex1_memory_check.r"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && args[1] == "--worker") {
  run_one(args[2], as.integer(args[3]))
  quit(save = "no")
}

# ---- driver ----------------------------------------------------------------------------------
rungs <- if (length(args) >= 1) as.integer(strsplit(args[1], ",")[[1]]) else c(2000L, 8000L, 12000L)
n_workers <- if (length(args) >= 2) as.integer(args[2]) else 20L
ram <- total_ram_gb()

cat("\n=== ex1 memory-mechanism check ===\n")
cat(sprintf("  host RAM      : %s\n", if (is.na(ram)) "unknown" else sprintf("%.1f GB", ram)))
cat(sprintf("  RSS measure   : %s\n",
  if (is_true_peak) "VmHWM (true peak)" else "ps RSS (CURRENT, understates the peak -- prefer Linux)"))
cat(sprintf("  rungs         : %s\n", paste(rungs, collapse = ", ")))
cat(sprintf("  projecting for: %d concurrent workers\n\n", n_workers))
cat("  'bare'  = extract delta/cor_sq via $draws()          <- the pre-fix code\n")
cat("  'named' = extract delta/cor_sq via $draws('delta')   <- the committed fix\n\n")

rows <- list()
for (it in rungs) for (md in c("bare", "named")) {
  cat(sprintf("  running %-5s iter=%-6d ... ", md, it)); flush.console()
  t0 <- Sys.time()
  out <- suppressWarnings(system2("Rscript", c(THIS_SCRIPT, "--worker", md, it),
    stdout = TRUE, stderr = FALSE))
  line <- grep("^RESULT ", out, value = TRUE)
  if (!length(line)) { cat("FAILED (no RESULT line)\n"); next }
  p <- strsplit(line[1], " +")[[1]]
  rows[[length(rows) + 1]] <- data.frame(mode = p[2], iter = as.integer(p[3]),
    heap_mb = as.numeric(p[4]), b_per_draw = as.numeric(p[5]), rhat_M = as.numeric(p[6]),
    rss_mb = as.numeric(p[7]), mins = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  cat(sprintf("heap peak %6.0f MB   (proc RSS %6.0f MB)  (%5.1f min)\n",
    as.numeric(p[4]), as.numeric(p[7]), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

r <- do.call(rbind, rows)
cat("\n\n=== RESULTS ===\n")
cat(sprintf("%-8s %-8s %14s %14s %12s %10s\n",
  "iter", "mode", "R heap MB", "bytes/draw", "proc RSS MB", "rhat_M"))
for (i in seq_len(nrow(r)))
  cat(sprintf("%-8d %-8s %14.0f %14.0f %12.0f %10.3f\n",
    r$iter[i], r$mode[i], r$heap_mb[i], r$b_per_draw[i], r$rss_mb[i], r$rhat_M[i]))
cat("\n  R heap  = gc() high-water mark over the extraction. Deterministic; THE metric to judge by.\n")
cat(sprintf("  proc RSS = %s\n",
  if (is_true_peak) "VmHWM, a true process peak. Use for the worker-count projection."
  else "CURRENT RSS (macOS): noisy, GC-timing dependent. IGNORE on this host."))

cat("\n=== VERDICT ===\n")
for (it in unique(r$iter)) {
  b <- r$heap_mb[r$iter == it & r$mode == "bare"]
  n <- r$heap_mb[r$iter == it & r$mode == "named"]
  if (!length(b) || !length(n)) next
  cat(sprintf("  iter=%-6d bare %5.0f MB vs named %5.0f MB  -> saved %5.0f MB (%.1fx)\n",
    it, b, n, b - n, b / n))
}
big <- r[r$iter == max(r$iter), ]
if (nrow(big) == 2) {
  bb <- big$heap_mb[big$mode == "bare"]; nn <- big$heap_mb[big$mode == "named"]
  cat(sprintf("\n  At the top rung (iter=%d), %d workers would need:\n", max(r$iter), n_workers))
  cat(sprintf("    pre-fix : %6.1f GB   %s\n", n_workers * bb / 1024,
    if (!is.na(ram) && n_workers * bb / 1024 > 0.8 * ram) "<-- EXCEEDS 80% OF RAM (this is the OOM)" else ""))
  cat(sprintf("    post-fix: %6.1f GB   %s\n", n_workers * nn / 1024,
    if (!is.na(ram) && n_workers * nn / 1024 > 0.8 * ram) "<-- STILL TOO HIGH, reduce workers or rungs" else "(fits)"))
  if (!is.na(ram))
    cat(sprintf("\n  Safe worker count post-fix at this rung (80%% of RAM): %d\n",
      floor(0.8 * ram * 1024 / nn)))
}
cat("\n  Hypothesis is SUPPORTED if bare/named grows with iterations and the ratio is large at the\n")
cat("  top rung. It is REFUTED if the two arms are close, in which case the OOM lies elsewhere and\n")
cat("  the ladder change was the wrong remedy.\n")

# The second question this answers, which matters just as much as the first: the bare/named gap is
# the BUG, but the `named` column is the FLOOR -- what a fit costs even with the bug fixed. On a
# macOS laptop trial (rungs 1000 and 3000) the named arm still grew from 355 MB to 756 MB, i.e.
# roughly 200 MB fixed plus ~65 KB per stored draw, far more than the ~10 KB/draw the cached draws
# array alone accounts for. If that slope holds on the desktop, round 3 at 12000 is still large
# enough that the worker count, not just the ladder, needs setting deliberately. Read the
# "safe worker count" line above as the operational answer.
if (nrow(big) == 2) {
  nmd <- r[r$mode == "named", ]
  if (nrow(nmd) >= 2) {
    fit <- lm(heap_mb ~ iter, data = nmd)
    cat(sprintf("\n  Post-fix per-fit cost scales as %.0f MB + %.1f KB per sampling iteration.\n",
      coef(fit)[1], 1024 * coef(fit)[2]))
    for (target in c(8000, 12000, 20000)) {
      pred <- coef(fit)[1] + coef(fit)[2] * target
      cat(sprintf("    extrapolated to iter=%-6d: %6.0f MB/fit -> %5.1f GB for %d workers%s\n",
        target, pred, n_workers * pred / 1024, n_workers,
        if (!is.na(ram) && n_workers * pred / 1024 > 0.8 * ram) "   <-- TOO HIGH" else ""))
    }
    cat("  (extrapolation from the rungs actually measured; run the real rungs for real numbers)\n")
  }
}
cat("\n")
