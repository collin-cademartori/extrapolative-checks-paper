## Confirm the mechanism: do divergences localize to SMALL row-norm ||Lambda_raw[n,:]|| (the funnel neck)?
## Uses the saved unit-norm est-sigma fit (no re-run). Per draw compute each row's norm; test whether
## divergent draws have systematically smaller norms (min over rows, and per row). Mechanism predicts
## divergent draws sit at small ||Lambda_raw|| for some row.
library(posterior)
SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
z <- readRDS(file.path(SP,"div_subspace_fit.rds")); d <- z$d; sdg <- z$sdg; M <- 8L; K <- z$K
divflag <- as.logical(as.numeric(extract_variable_matrix(sdg,"divergent__")))
vmat <- function(v) as.numeric(extract_variable_matrix(d, v))

# per-draw row norms ||Lambda_raw[n,1:min(K,n)]||
rn <- sapply(1:M, function(n){ ks <- 1:min(K,n); m <- sapply(ks, function(k) vmat(sprintf("Lambda_raw[%d,%d]", n, k))); sqrt(rowSums(as.matrix(m)^2)) })
colnames(rn) <- paste0("row", 1:M)
minrn <- apply(rn, 1, min)                 # smallest row norm per draw (the most funnel-prone)
logminrn <- log(minrn)

cat(sprintf("divergences: %d of %d (%.1f%%)\n\n", sum(divflag), length(divflag), 100*mean(divflag)))
cat("=== does SMALL ||Lambda_raw|| localize divergences? ===\n")
cat(sprintf("min-row-norm:  bulk mean=%.3f  DIV mean=%.3f  | z=%.2f\n", mean(minrn[!divflag]), mean(minrn[divflag]),
  (mean(minrn[divflag])-mean(minrn[!divflag]))/sd(minrn)))
cat(sprintf("log(min-norm): bulk mean=%.3f  DIV mean=%.3f  | z=%.2f  (log sharpens the funnel neck)\n",
  mean(logminrn[!divflag]), mean(logminrn[divflag]), (mean(logminrn[divflag])-mean(logminrn[!divflag]))/sd(logminrn)))
cat(sprintf("P(min-norm < 0.5): bulk=%.3f  DIV=%.3f   |   P(min-norm<0.3): bulk=%.3f DIV=%.3f\n",
  mean(minrn[!divflag]<0.5), mean(minrn[divflag]<0.5), mean(minrn[!divflag]<0.3), mean(minrn[divflag]<0.3)))

cat("\n=== per-row: is a small norm in row n associated with divergences? (z = div-mean vs bulk, in sd) ===\n")
for(n in 1:M){ x<-rn[,n]; cat(sprintf("  row%d norm: bulk=%.2f div=%.2f z=%+.2f | P(<0.5) bulk=%.3f div=%.3f\n",
  n, mean(x[!divflag]), mean(x[divflag]), (mean(x[divflag])-mean(x[!divflag]))/sd(x), mean(x[!divflag]<0.5), mean(x[divflag]<0.5))) }

# contrast: does log(min-norm) beat the earlier whitened linear separation (0.61)?
cat(sprintf("\nfor scale: earlier whitened LINEAR separation was 0.61 sd; log(min-norm) alone gives z=%.2f\n",
  (mean(logminrn[divflag])-mean(logminrn[!divflag]))/sd(logminrn)))
cat("Read: if DIV draws sit at markedly smaller min-row-norm (negative z, higher P(<0.5)) -> funnel mechanism confirmed.\n")
