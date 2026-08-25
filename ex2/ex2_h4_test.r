## Test H4: divergences concentrate in the treated unit's counterfactual-extrapolation region --
## factor-1 innovations over the (late/post) treatment window, confounded with delta (and sigma_1).
## On the saved fit: (1) whitened separation of divergent vs bulk restricted to candidate SUBSETS --
## if a small treatment-window-factor1 + delta subset matches the full-param 0.61, the divergence
## direction lives there. (2) the delta <-> Phi[post,1] trade-off (the counterfactual ridge) and its
## posterior spread. Treatment window = t 26..30 (T=30, num_treated=5).
library(posterior)
SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
z <- readRDS(file.path(SP,"div_subspace_fit.rds")); d <- z$d; sdg <- z$sdg; T_<-z$T_; K<-z$K
divflag <- as.logical(as.numeric(extract_variable_matrix(sdg,"divergent__")))
V <- function(v) as.numeric(extract_variable_matrix(d, v))
sep <- function(vars){ vars <- intersect(vars, dimnames(d)$variable); if(length(vars)<1) return(NA)
  X <- sapply(vars, V); ctr <- scale(X); keep<-is.finite(colSums(ctr)); ctr<-ctr[,keep,drop=FALSE]
  dm <- colMeans(ctr[divflag,,drop=FALSE])-colMeans(ctr[!divflag,,drop=FALSE])
  S<-cor(ctr); lam<-0.05*mean(diag(S)); dir<-solve(S+lam*diag(ncol(S)),dm); dir<-dir/sqrt(sum(dir^2))
  proj<-as.numeric(ctr%*%dir); (mean(proj[divflag])-mean(proj[!divflag]))/sd(proj) }

tw <- 26:30; late <- 20:30; early <- 1:10
cat(sprintf("divergences: %.1f%%\n\n=== whitened divergence separation by parameter SUBSET (full-param ref ~0.61) ===\n", 100*mean(divflag)))
sets <- list(
  "ALL params (ref)"                 = grep("^(Phi_innovations|sigma_raw|rho|tau_param|Lambda_raw|delta_raw|omega_sq_param|Phi_means_param)\\[", dimnames(d)$variable, value=TRUE),
  "tw factor1 innov + delta + sig1"  = c(sprintf("Phi_innovations[%d,1]",tw), sprintf("delta_raw[%d]",1:5), "sigma_raw[1]"),
  "tw factor1 innov ONLY (26-30)"    = sprintf("Phi_innovations[%d,1]",tw),
  "late factor1 innov (20-30)"       = sprintf("Phi_innovations[%d,1]",late),
  "early factor1 innov (1-10)"       = sprintf("Phi_innovations[%d,1]",early),
  "ALL factor1 innov"                = sprintf("Phi_innovations[%d,1]",1:T_),
  "factor2+3 innov (all t)"          = c(sprintf("Phi_innovations[%d,2]",1:T_), sprintf("Phi_innovations[%d,3]",1:T_)),
  "delta_raw only"                   = sprintf("delta_raw[%d]",1:5),
  "sigma_raw only"                   = sprintf("sigma_raw[%d]",1:8))
for(nm in names(sets)) cat(sprintf("  %-32s sep=%.2f  (%d params)\n", nm, sep(sets[[nm]]), length(intersect(sets[[nm]],dimnames(d)$variable))))

cat("\n=== counterfactual ridge: delta_raw[j] vs factor-1 path Phi[post,1] (trade-off?) ===\n")
Phi_post <- tryCatch(sapply(tw, function(t) V(sprintf("Phi[%d,1]",t))), error=function(e) NULL)
if(!is.null(Phi_post)){
  for(j in 1:5) cat(sprintf("  delta_raw[%d] <-> Phi[%d,1]: corr=%+.2f | Phi[%d,1] post-sd=%.2f\n",
    j, tw[j], cor(V(sprintf("delta_raw[%d]",j)), Phi_post[,j]), tw[j], sd(Phi_post[,j])))
  # posterior sd of factor-1 path: is the POST window much wider (less identified) than PRE?
  Phi_pre <- sapply(early, function(t) V(sprintf("Phi[%d,1]",t)))
  cat(sprintf("\n  factor-1 path posterior sd: PRE(t1-10) mean=%.2f  POST(t26-30) mean=%.2f  (wider post = weaker counterfactual ID)\n",
    mean(apply(Phi_pre,2,sd)), mean(apply(Phi_post,2,sd))))
} else cat("  (Phi not in saved draws)\n")
cat("\nRead: H4 confirmed if the small tw-factor1+delta subset separation ~ the full 0.61 (direction lives there), delta<->Phi anti-correlate, and post factor path is wide.\n")
