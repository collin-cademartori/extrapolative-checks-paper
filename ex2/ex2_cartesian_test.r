## INTERVENTION (not correlation): is tau's variability the amplifier behind the unit-norm + est-sigma
## divergences? Mechanism: logL = -0.5*sum[Y/(tau*sigma) - (u.Phi)/tau]^2 - T*log(tau*sigma), so the
## curvature of the whole (loadings, factors) block carries 1/tau^2 (~190x at tau~0.073). With sigma
## FIXED, tau is pinned by the residuals; with sigma FREE, tau trades off against sigma (measured
## corr(log tau, log sigma) = -0.23), so the block's curvature MOVES between iterations -- which one
## adapted step size + diagonal metric cannot track.
## Test: unit-norm, sigma FREE throughout, tau free vs tau FIXED at its posterior mean (tau_val>0).
## If fixing tau collapses divergences (~408 -> small) with sigma still free, the amplifier is tau.
library(posterior); source("../sample_model.r"); source("../pathfinder_init.r")
ife_mod <- cmdstan_model("../ife_named_cartesian.stan")
# this model's loadings parameter is Lambda (Cartesian), not Lambda_raw
PF_PARAM_BASES <- c("sigma_raw","tau_param","omega_sq_param","Phi_innovations","Phi_means_param","rho","Lambda","gamma_raw","delta_raw")
SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
tau_hat <- mean(as.numeric(extract_variable_matrix(readRDS(file.path(SP,"div_subspace_fit.rds"))$d, "tau")))
cat(sprintf("posterior-mean tau from the problematic fit = %.4f\n", tau_hat))

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
base <- list(M_units=N_, T_times=T_, K_latent=3L, Y=fy, a_rho=90, b_rho=10, m_tau=0.1, s_tau=0.05,
  sigma_data=apply(fy,2,function(x)sqrt(mean(x^2))), fit_overall_scales=0L, nonstationary=0L, unit_intercepts=0L,
  factor_means=1L, sample_posterior=1L, num_treated=5L, gamma_scale=3, gamma_loc=4, alpha_diag=10)

run <- function(tau_val, label) {
  sdta <- c(base, list(tau_val=tau_val))
  inits <- pathfinder_inits(ife_mod, sdta, 6, seed=42, quiet=TRUE)
  t0 <- Sys.time()
  f <- ife_mod$sample(data=sdta, chains=6, parallel_chains=6, iter_warmup=1000, iter_sampling=800,
    adapt_delta=0.9, max_treedepth=12, init=if(is.null(inits)) 2 else inits, seed=42, refresh=0,
    show_messages=FALSE, show_exceptions=FALSE)
  dg <- f$diagnostic_summary(quiet=TRUE); sdg <- f$sampler_diagnostics()
  td <- as.numeric(extract_variable_matrix(sdg,"treedepth__"))
  d <- f$draws()
  sd_logtau <- if(tau_val==0) sd(log(as.numeric(extract_variable_matrix(d,"tau")))) else NA
  sig <- sapply(1:N_, function(n) as.numeric(extract_variable_matrix(d, sprintf("scale_est[%d]",n))))
  tau_mean <- mean(as.numeric(extract_variable_matrix(d,"tau")))
  rh <- tryCatch(max(f$summary(c("delta","sigma"),"rhat")$rhat,na.rm=TRUE), error=function(e) NA)
  cat(sprintf("%-26s div=%-5d tree=%-4d td_mean=%.1f rhat=%.2f | %s | mean scale_est=%.2f | %.0f s\n",
    label, sum(dg$num_divergent), sum(dg$num_max_treedepth), mean(td), rh,
    sprintf("tau=%.3f",tau_mean), mean(colMeans(sig)),
    as.numeric(difftime(Sys.time(), t0, units="secs"))))
}
cat("\n=== tau intervention (unit-norm, sigma FREE in both arms, ad=0.9, 6 chains) ===\n")
cat("baseline for reference: unitnorm+est tau-free previously gave div=408, tree=138\n\n")
run(0, "CARTESIAN (fixed sigma_data)")
cat("\nRead: if fixing tau collapses divergences while sigma stays free -> tau's variability (1/tau^2 curvature scaling) is the amplifier.\n")
