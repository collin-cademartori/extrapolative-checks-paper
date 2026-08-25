## Divergence geometry, done right (Q1). Marginal shifts were tiny; look at the SUBSPACE and test the
## innovations claim. Two analyses on one correctly-seeded est-sigma no_ints normal(0,1) fit:
## (A) Are the SAMPLED innovations Phi_innovations posterior-autocorrelated across time? (Non-centering
##     should keep them ~iid; the REALIZATIONS Phi are autocorrelated by construction -- shown as a
##     reference, NOT the sampled coordinates.) posterior sd of each innovation too.
## (B) Whitened/Mahalanobis direction separating DIVERGENT from bulk draws (ridge-LDA over all params),
##     to find a non-axis-aligned culprit + its top-contributing parameters.
library(posterior); source("../sample_model.r"); source("../pathfinder_init.r")
ife_mod <- cmdstan_model("../ife_named_unitnorm_ig_s1.stan"); stopifnot("Lambda_raw" %in% PF_PARAM_BASES)
SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"

ruv <- function(d){v<-rnorm(d); v/sqrt(sum(v*v))}
anchor_order <- function(y,K){N<-ncol(y);yc<-scale(y,center=TRUE,scale=FALSE);sel<-1L;rem<-setdiff(seq_len(N),sel)
  while(length(sel)<K&&length(rem)>0){Q<-qr.Q(qr(yc[,sel,drop=FALSE]));rs<-yc[,rem,drop=FALSE]-Q%*%crossprod(Q,yc[,rem,drop=FALSE]);pk<-rem[which.max(colSums(rs^2))];sel<-c(sel,pk);rem<-setdiff(rem,pk)};c(sel,rem)}
simmod<-function(N_unc,N_spur,sim=0.9,T_=30,Tt=5){Nu<-1+2+N_spur+N_unc;Kg<-3
  ft<-6+arima.sim(list(ar=0.9),n=T_);fa<-(ft-6)+c(rep(0,T_-Tt),rep(-1.9,Tt));fu<-matrix(nrow=1,ncol=T_);cu<-Inf
  while(cu>0.01){fu[1,]<-rnorm(1,0,2)+arima.sim(list(ar=0.9),n=T_);cu<-max(abs(cor(t(fu),t(t(ft)))))}
  fc<-rbind(ft,fa,fu);ld<-matrix(nrow=Nu,ncol=Kg);ld[1,]<-c(1,0,0)
  for(n in 1:2)ld[1+n,]<-c(sqrt(sim),0,sqrt(1-sim)*ruv(1));for(n in seq_len(N_spur))ld[3+n,]<-c(0,sqrt(sim),sqrt(1-sim)*ruv(1));for(n in seq_len(N_unc))ld[3+N_spur+n,]<-c(0,0,ruv(1))
  lat<-ld%*%fc;t(lat+rnorm(nrow(lat)*ncol(lat),sd=0.1*max(apply(lat,1,sd))))}

set.seed(89333); ys<-simmod(2,3); perm<-anchor_order(ys,3); fy<-ys[,perm]; N_<-ncol(fy); T_<-nrow(fy)
sd_ <- list(M_units=N_, T_times=T_, K_latent=3L, Y=fy, a_rho=90, b_rho=10, tau_val=0, m_tau=0.1, s_tau=0.05,
  sigma_data=apply(fy,2,function(x)sqrt(mean(x^2))), fit_overall_scales=1L, nonstationary=0L, unit_intercepts=0L,
  factor_means=1L, sample_posterior=1L, num_treated=5L, gamma_scale=3, gamma_loc=4, alpha_diag=10)
inits <- pathfinder_inits(ife_mod, sd_, 6, seed=42, quiet=TRUE)
f <- ife_mod$sample(data=sd_, chains=6, parallel_chains=6, iter_warmup=1000, iter_sampling=800, adapt_delta=0.9,
  max_treedepth=12, init=inits, seed=42, refresh=0, show_messages=FALSE, show_exceptions=FALSE)
d <- f$draws(); sdg <- f$sampler_diagnostics()
saveRDS(list(d=d, sdg=sdg, T_=T_, K=3L), file.path(SP,"div_subspace_fit.rds"))
divflag <- as.logical(as.numeric(extract_variable_matrix(sdg,"divergent__")))
cat(sprintf("divergences: %d of %d (%.1f%%)\n\n", sum(divflag), length(divflag), 100*mean(divflag)))

vmat <- function(v) as.numeric(extract_variable_matrix(d, v))

cat("=== (A) innovation vs realization posterior autocorrelation (per factor k) ===\n")
cat("innov = SAMPLED Phi_innovations; Phi = transformed AR realization (reference). lag-j = mean over t of cor across draws.\n")
for(k in 1:3){
  I <- sapply(1:T_, function(t) vmat(sprintf("Phi_innovations[%d,%d]", t, k)))   # [draws x T]
  P <- tryCatch(sapply(1:T_, function(t) vmat(sprintf("Phi[%d,%d]", t, k))), error=function(e) NULL)
  ac <- function(M, lag) mean(sapply(1:(T_-lag), function(t) cor(M[,t], M[,t+lag])))
  sd_innov <- mean(apply(I,2,sd))
  line <- sprintf("  k=%d innov: sd=%.2f | lag1=%.2f lag2=%.2f lag3=%.2f", k, sd_innov, ac(I,1), ac(I,2), ac(I,3))
  if(!is.null(P)) line <- paste0(line, sprintf("   || Phi realiz: lag1=%.2f lag2=%.2f", ac(P,1), ac(P,2)))
  cat(line, "\n")
}

cat("\n=== (B) whitened separating direction (ridge-LDA): divergent vs bulk over ALL params ===\n")
allv <- grep("^(Phi_innovations|sigma_raw|rho|tau_param|Lambda_raw|delta_raw|omega_sq_param|Phi_means_param)\\[", dimnames(d)$variable, value=TRUE)
X <- sapply(allv, vmat); ctr <- scale(X)  # standardize
keep <- is.finite(colSums(ctr)); ctr <- ctr[,keep]; allv <- allv[keep]
dmean <- colMeans(ctr[divflag,,drop=FALSE]) - colMeans(ctr[!divflag,,drop=FALSE])
S <- cor(ctr); lam <- 0.05*mean(diag(S)); dir <- solve(S + lam*diag(ncol(S)), dmean); dir <- dir/sqrt(sum(dir^2))
proj <- as.numeric(ctr %*% dir)
sep <- (mean(proj[divflag]) - mean(proj[!divflag]))/sd(proj)
cat(sprintf("MARGINAL best |z| shift = %.2f   vs   WHITENED-direction separation = %.2f sd\n", max(abs(dmean)), sep))
contrib <- dir * dmean   # contribution of each param to the separation
ord <- order(-abs(contrib))
cat("top parameters loading the divergence direction (contribution = dir*dmean):\n")
for(j in head(ord,18)) cat(sprintf("  %-20s dir=%+.2f dmean=%+.2f contrib=%+.3f\n", allv[j], dir[j], dmean[j], contrib[j]))
cat(sprintf("\nsaved fit to div_subspace_fit.rds for re-analysis. (marginal was ~0.09; whitened sep=%.2f tells the real story)\n", sep))
