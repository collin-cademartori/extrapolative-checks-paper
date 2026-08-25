## Test H1 (factor-loading funnel / near-degenerate factor) + H2 (sigma coupling) on the saved unit-norm
## est-sigma fit. In the SAMPLED coordinates: does a factor's innovations (Phi_innovations[:,k]) inflate
## toward the prior as its loadings shrink, and do divergences localize at the neck?
## Measures per draw: (a) smallest singular value of the loading matrix (near-degeneracy),
## (b) per-factor raw loading energy Ek=sum_n Lambda_raw[n,k]^2 and innovation energy Ik=mean_t innov[t,k]^2,
## (c) the "most-freed factor" = max_k Ik. Then localize divergences against each + funnel correlations.
library(posterior)
SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
z <- readRDS(file.path(SP,"div_subspace_fit.rds")); d <- z$d; sdg <- z$sdg; M<-8L; K<-z$K; T_<-z$T_
divflag <- as.logical(as.numeric(extract_variable_matrix(sdg,"divergent__")))
V <- function(v) as.numeric(extract_variable_matrix(d, v)); nd <- length(divflag)
zsh <- function(x) (mean(x[divflag])-mean(x[!divflag]))/sd(x)

# raw loading matrix per draw (lower-triangular cholesky_factor_cov): Lambda_raw[n, 1:min(K,n)]
Lraw <- array(0, c(nd, M, K)); for(n in 1:M) for(k in 1:min(K,n)) Lraw[,n,k] <- V(sprintf("Lambda_raw[%d,%d]",n,k))
Lnrm <- array(0, c(nd, M, K)); for(n in 1:M) for(k in 1:min(K,n)) Lnrm[,n,k] <- V(sprintf("Lambda[%d,%d]",n,k))
Innov <- array(0, c(nd, T_, K)); for(t in 1:T_) for(k in 1:K) Innov[,t,k] <- V(sprintf("Phi_innovations[%d,%d]",t,k))

Ek_raw <- t(sapply(1:nd, function(i) colSums(Lraw[i,,]^2)))       # [nd x K] loadings energy per factor (raw)
Ek_nrm <- t(sapply(1:nd, function(i) colSums(Lnrm[i,,]^2)))       # normalized loadings energy per factor (identified use)
Ik     <- t(sapply(1:nd, function(i) colMeans(Innov[i,,]^2)))     # [nd x K] innovation mean-square per factor (~1 if freed)
smin   <- sapply(1:nd, function(i) min(svd(Lnrm[i,,])$d))          # smallest singular value of the loading matrix
minEk  <- apply(Ek_nrm, 1, min)                                    # least-used factor's identified energy
maxIk  <- apply(Ik, 1, max)                                        # most-freed factor's innovation energy

cat(sprintf("divergences: %d of %d (%.1f%%)\n\n", sum(divflag), nd, 100*mean(divflag)))
cat("=== degeneracy / freed-factor localization of divergences (z = div vs bulk, in sd) ===\n")
cat(sprintf("  smallest singular value of Lambda : z=%+.2f | bulk=%.3f div=%.3f  (H1: near-degenerate -> small)\n", zsh(smin), mean(smin[!divflag]), mean(smin[divflag])))
cat(sprintf("  log(smin)                          : z=%+.2f\n", zsh(log(smin))))
cat(sprintf("  min factor identified-energy       : z=%+.2f | bulk=%.3f div=%.3f\n", zsh(minEk), mean(minEk[!divflag]), mean(minEk[divflag])))
cat(sprintf("  max factor innovation-energy (freed): z=%+.2f | bulk=%.3f div=%.3f  (H1: freed -> inflates ~1)\n", zsh(maxIk), mean(maxIk[!divflag]), mean(maxIk[divflag])))

cat("\n=== per-factor: identified loading energy, innovation energy, funnel corr, div localization ===\n")
for(k in 1:K){
  cat(sprintf("  factor%d: E_nrm bulk=%.2f div=%.2f (z=%+.2f) | innov E bulk=%.2f div=%.2f (z=%+.2f) | corr(log E_nrm, innovE)=%+.2f\n",
    k, mean(Ek_nrm[!divflag,k]), mean(Ek_nrm[divflag,k]), zsh(Ek_nrm[,k]),
    mean(Ik[!divflag,k]), mean(Ik[divflag,k]), zsh(Ik[,k]), cor(log(Ek_nrm[,k]), Ik[,k])))
}

# funnel test: does the LEAST-loaded factor per draw have inflated innovations, esp. at divergences?
freed_k <- apply(Ek_nrm, 1, which.min)   # which factor is least identified in each draw
Ik_of_freed <- Ik[cbind(1:nd, freed_k)]
cat(sprintf("\ninnovation energy OF the least-loaded factor: bulk=%.2f div=%.2f (z=%+.2f)\n",
  mean(Ik_of_freed[!divflag]), mean(Ik_of_freed[divflag]), zsh(Ik_of_freed)))
cat(sprintf("corr( min-loading-energy , its innovation-energy ) = %+.2f   (strong negative = the funnel)\n", cor(minEk, Ik_of_freed)))
cat("\nRead: H1 confirmed if divergences sit at small smin / small min loading-energy AND inflated freed-factor innovations, with a strong negative funnel corr.\n")
