## Rebuild a simulation study's saved result object from its per-task checkpoint files.
##
## Both studies write one .rds per task as it completes (see the checkpointing block in
## ex1_sim_study.r / ex2_sim_study.r) and only assemble the final .RData at the very end. If a run
## dies before finishing -- disk, memory, anything -- the checkpoints are the surviving record, and
## this reassembles them into exactly what the driver would have written.
##
## Why this is exact rather than approximate: the driver's result is
## `bind_rows` over the per-task rows, and foreach returns results in TASK order regardless of the
## order workers finish. The checkpoint filenames lead with a zero-padded task index, so sorting
## them lexically reproduces that same order. Verified by identical() against a completed run.
##
## USAGE
##   Rscript rebuild_from_checkpoints.r <ckpt_dir> [out.RData]
##
##   # ex1
##   Rscript rebuild_from_checkpoints.r ex1/ckpt_ns_seed40318_reps200 ex1/sim_study_ns.RData
##   # ex2
##   Rscript rebuild_from_checkpoints.r ex2/ckpt_ints_seed52918_n800 ex2/sim_study_ints.RData
##
## The object name written into the .RData is inferred from the checkpoint prefix
## (`rep_` -> sim_study_stat, `task_` -> sim_study_ints), so the existing summary scripts load it
## unchanged. A partial rebuild is written normally and simply has fewer rows; the report below says
## which tasks are missing so you know what the summaries are conditioning on.

suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript rebuild_from_checkpoints.r <ckpt_dir> [out.RData]")
ckpt_dir <- args[1]
if (!dir.exists(ckpt_dir)) stop("no such directory: ", ckpt_dir)

files <- sort(list.files(ckpt_dir, pattern = "\\.rds$", full.names = TRUE))
if (!length(files)) stop("no .rds checkpoints in ", ckpt_dir)

kind <- if (any(grepl("^task_", basename(files)))) "ints" else "stat"
obj_name <- if (kind == "ints") "sim_study_ints" else "sim_study_stat"
out_file <- if (length(args) >= 2) args[2] else
  file.path(dirname(ckpt_dir), if (kind == "ints") "sim_study_ints.RData" else "sim_study_ns.RData")

rows <- lapply(files, function(f)
  tryCatch(readRDS(f), error = function(e) {
    warning("unreadable checkpoint (skipped): ", basename(f), " -- ", conditionMessage(e))
    NULL
  }))
bad <- vapply(rows, is.null, logical(1))
res <- bind_rows(rows[!bad])

cat(sprintf("\n=== rebuild from %s ===\n", ckpt_dir))
cat(sprintf("  checkpoints found : %d\n", length(files)))
if (any(bad)) cat(sprintf("  UNREADABLE        : %d (skipped)\n", sum(bad)))
cat(sprintf("  rows assembled    : %d\n", nrow(res)))
cat(sprintf("  columns           : %d\n", ncol(res)))

# Which tasks are missing: the index is the first number in the filename, and the driver numbers
# tasks 1..N contiguously, so gaps are exactly the tasks that never completed.
idx <- as.integer(sub("^[a-z]+_0*(\\d+)_.*$", "\\1", basename(files)))
idx <- sort(idx[!is.na(idx)])
if (length(idx)) {
  full <- seq_len(max(idx))
  missing <- setdiff(full, idx)
  cat(sprintf("  task index range  : %d..%d\n", min(idx), max(idx)))
  if (length(missing)) {
    cat(sprintf("  MISSING           : %d of %d in that range\n", length(missing), length(full)))
    cat(sprintf("                      %s%s\n", paste(head(missing, 25), collapse = ", "),
      if (length(missing) > 25) sprintf(" ... (+%d more)", length(missing) - 25) else ""))
    cat("  NOTE: a run that died partway also never reached its later tasks, so the true\n")
    cat("        denominator is the intended rep count, not max(task index).\n")
  } else {
    cat("  MISSING           : none in that range\n")
  }
}

if ("failed" %in% names(res)) {
  nf <- sum(res$failed, na.rm = TRUE)
  cat(sprintf("  failed tasks      : %d\n", nf))
  if (nf > 0) {
    cat("  failure messages:\n")
    for (m in unique(res$error[which(res$failed)])) cat(sprintf("    - %s\n", m))
  }
}

assign(obj_name, res)
save(list = obj_name, file = out_file)
cat(sprintf("\n  wrote %s  (object `%s`) -- the summary script loads this unchanged.\n\n",
  out_file, obj_name))
