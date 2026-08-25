## Isolate the divergence cause (user's hint: fixed-sigma + NON-normalized Lambda had few divergences).
## 2x2 on the SAME dataset: {unit-norm (normalized), original (non-normalized)} x {est sigma, fixed sigma}.
## Same seed/chains/ad; pathfinder per-model (Lambda_raw for unit-norm, Lambda for original). Report
## divergences + treedepth + rhat. If normalization is the culprit: unit-norm rows high regardless of
## est/fixed; if est-sigma: est rows high; if interaction: only unit-norm+est high.
library(posterior); source("../sample_model.r"); source("../pathfinder_init.r")
uni <- cmdstan_model("../ife_named_unitnorm_ig_s1.stan")   # normalized, sigma_raw~normal(0,1)
orig <- cmdstan_model("../ife_named.stan")                  # original (non-normalized) loadings
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
base <- list(M_units=N_, T_times=T_, K_latent=3L, Y=fy, a_rho=90, b_rho=10, tau_val=0, m_tau=0.1, s_tau=0.05,
  sigma_data=apply(fy,2,function(x)sqrt(mean(x^2))), nonstationary=0L, unit_intercepts=0L, factor_means=1L,
  sample_posterior=1L, num_treated=5L, gamma_scale=3, gamma_loc=4, alpha_diag=10)
run <- function(mod, est, lam_base) {
  assign("PF_PARAM_BASES", c("sigma_raw","tau_param","omega_sq_param","Phi_innovations","Phi_means_param","rho",lam_base,"gamma_raw","delta_raw"), envir=globalenv())
  sdta <- c(base, list(fit_overall_scales=as.integer(est)))
  inits <- pathfinder_inits(mod, sdta, 6, seed=42, quiet=TRUE)
  f <- mod$sample(data=sdta, chains=6, parallel_chains=6, iter_warmup=1000, iter_sampling=800, adapt_delta=0.9,
    max_treedepth=12, init=if(is.null(inits)) 2 else inits, seed=42, refresh=0, show_messages=FALSE, show_exceptions=FALSE)
  dg <- f$diagnostic_summary(quiet=TRUE); td <- as.numeric(extract_variable_matrix(f$sampler_diagnostics(),"treedepth__"))
  rh <- tryCatch(max(f$summary(c("tau","delta","sigma"),"rhat")$rhat,na.rm=TRUE), error=function(e) NA)
  c(div=sum(dg$num_divergent), tree=sum(dg$num_max_treedepth), td_mean=round(mean(td),1), rhat=round(rh,2))
}
cat("\n2x2 divergence isolation (same dataset, 6 chains, ad=0.9, no_ints, est/fixed sigma):\n\n")
cat(sprintf("%-22s %s\n", "config", "div  tree  td_mean  rhat"))
for(cfg in list(c("unitnorm","est"), c("unitnorm","fix"), c("original","est"), c("original","fix"))) {
  mod <- if(cfg[1]=="unitnorm") uni else orig; lb <- if(cfg[1]=="unitnorm") "Lambda_raw" else "Lambda"
  r <- run(mod, cfg[2]=="est", lb)
  cat(sprintf("%-22s %-4d %-5d %-8.1f %s\n", paste(cfg,collapse="+"), r["div"], r["tree"], r["td_mean"], r["rhat"]))
}
cat("\nRead: normalized rows (unit-norm) high regardless of est/fix -> normalization; est rows high -> est-sigma; only unitnorm+est high -> interaction.\n")
