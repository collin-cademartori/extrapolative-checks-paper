## Localize the residual divergences (Q1): correctly-seeded est-sigma no_ints normal(0,1) fit still
## diverges (~509 at ad=0.9, a handful at ad=0.99). Which parameters go extreme at divergent transitions?
## Also check curvature suspects: rho near 1 (near-unit-root AR funnel), tau near 0, Lambda_raw diagonal
## near 0 (despite inv_gamma), omega_sq near bounds. Single-process; PF now seeds Lambda_raw (fixed).
library(posterior); source("../sample_model.r"); source("../pathfinder_init.r")
ife_mod <- cmdstan_model("../ife_named_unitnorm_ig_s1.stan")
stopifnot("Lambda_raw" %in% PF_PARAM_BASES)

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
divflag <- as.logical(as.numeric(extract_variable_matrix(sdg, "divergent__")))
cat(sprintf("divergences: %d of %d (%.1f%%) | treedepth mean=%.1f\n", sum(divflag), length(divflag),
  100*mean(divflag), mean(as.numeric(extract_variable_matrix(sdg,"treedepth__")))))

# which params are extreme at divergent transitions (standardized mean shift)
vars <- grep("^(sigma\\[|sigma_raw|tau|delta_raw|rho|Phi_means_param|omega_sq_param)", dimnames(d)$variable, value=TRUE)
# add the Lambda_raw diagonal entries (collapse suspects)
vars <- c(vars, sprintf("Lambda_raw[%d,%d]", 1:3, 1:3))
zex <- sapply(vars, function(v){x<-as.numeric(extract_variable_matrix(d,v)); (mean(x[divflag])-mean(x[!divflag]))/sd(x)})
cat("\n=== params most shifted at divergent draws (|z|, div-mean vs bulk) ===\n")
for(j in head(order(-abs(zex)),16)) cat(sprintf("  %-18s z=%+.2f | bulk_mean=%.3f div_mean=%.3f\n",
  vars[j], zex[j], mean(as.numeric(extract_variable_matrix(d,vars[j]))[!divflag]),
  mean(as.numeric(extract_variable_matrix(d,vars[j]))[divflag])))

cat("\n=== curvature suspects (posterior ranges) ===\n")
for(v in c("rho[1]","rho[2]","rho[3]","tau","omega_sq_param[1]","Lambda_raw[1,1]","Lambda_raw[2,2]","Lambda_raw[3,3]")) {
  x <- as.numeric(extract_variable_matrix(d,v)); cat(sprintf("  %-18s mean=%.3f q01=%.3f q99=%.3f min=%.3f max=%.3f\n",
    v, mean(x), quantile(x,.01), quantile(x,.99), min(x), max(x)))
}
