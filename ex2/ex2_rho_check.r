## Quick ex2 check: does the POSTERIOR rho match its Beta(90,10) prior mean (0.90)? The derived sigma
## multiple is 1/(E[realized AR sd] x pass-through) and E[realized sd] is very sensitive to rho near 1
## (T=30: rho .90 -> 0.723, .95 -> ~0.60, .97 -> 0.473), so the multiple must be derived at the rho the
## POSTERIOR actually uses, not the prior mean. Also reports the realized Phi amplitude directly (the
## quantity the multiple is really about) and tau/scale_est. Current config, Cartesian model, ad=0.8.
library(posterior)
source("../sample_model.r"); source("../pathfinder_init.r")
ife_mod <- cmdstan_model("../ife_named_cartesian.stan")
stopifnot(basename(ife_mod$stan_file()) == "ife_named_cartesian.stan")

ruv <- function(d){v<-rnorm(d); v/sqrt(sum(v*v))}
anchor_order <- function(y,K){N<-ncol(y);yc<-scale(y,center=TRUE,scale=FALSE);sel<-1L;rem<-setdiff(seq_len(N),sel)
  while(length(sel)<K&&length(rem)>0){Q<-qr.Q(qr(yc[,sel,drop=FALSE]));rs<-yc[,rem,drop=FALSE]-Q%*%crossprod(Q,yc[,rem,drop=FALSE]);pk<-rem[which.max(colSums(rs^2))];sel<-c(sel,pk);rem<-setdiff(rem,pk)};c(sel,rem)}
simmod<-function(N_unc,N_spur,sim=0.9,T_=30,Tt=5){Nu<-1+2+N_spur+N_unc;Kg<-3
  ft<-6+arima.sim(list(ar=0.9),n=T_);fa<-(ft-6)+c(rep(0,T_-Tt),rep(-1.9,Tt));fu<-matrix(nrow=1,ncol=T_);cu<-Inf
  while(cu>0.01){fu[1,]<-rnorm(1,0,2)+arima.sim(list(ar=0.9),n=T_);cu<-max(abs(cor(t(fu),t(t(ft)))))}
  fc<-rbind(ft,fa,fu);ld<-matrix(nrow=Nu,ncol=Kg);ld[1,]<-c(1,0,0)
  for(n in 1:2)ld[1+n,]<-c(sqrt(sim),0,sqrt(1-sim)*ruv(1));for(n in seq_len(N_spur))ld[3+n,]<-c(0,sqrt(sim),sqrt(1-sim)*ruv(1));for(n in seq_len(N_unc))ld[3+N_spur+n,]<-c(0,0,ruv(1))
  lat<-ld%*%fc;t(lat+rnorm(nrow(lat)*ncol(lat),sd=0.1*max(apply(lat,1,sd))))}
rms <- function(x) sqrt(mean(x^2))

K <- 3L; NPT <- 5L; run_seed <- 88213
cat("\nex2 posterior rho check | prior rho ~ Beta(90,10), mean 0.900, sd 0.030\n")
cat("prior-predictive realized amplitude at rho=0.90, T=30: sd 0.723, RMS 0.947\n\n")
for (rp in 1:2) for (N_comp in c(2,3)) {
  set.seed(run_seed + 1000*rp + 10*N_comp + 90)
  ys <- simmod(5-N_comp, N_comp, sim=0.9, T_=30)
  perm <- anchor_order(ys, K); fy <- ys[,perm]; N_ <- ncol(fy); T_ <- nrow(fy)
  for (arm in c("no_ints","ints")) {
    ints <- arm=="ints"
    os <- if (ints) apply(fy,2,sd) else apply(fy,2,rms)
    extra <- if (ints) list(include_ints=TRUE,int_scale=3,int_loc=4) else list(include_factor_means=TRUE)
    f <- tryCatch(do.call(sample_model, c(list(N_units=N_, T_times=T_, K_latent=K, overall_scales=os,
      err_scale=0, err_scale_mean=0.1, err_scale_sd=0.05, data=fy, autocor_a=90, autocor_b=10,
      nonstationary=FALSE, num_treated=NPT, fit_scales=0, alpha_diag=10, pathfinder_init=TRUE,
      type="posterior", quiet=TRUE, ad=0.8, iter=500, iter_warm=500, n_chains=4, seed=42,
      parallel_chains=4, return_draws=c("tau","rho","Phi","scale_est")), extra)),
      error=function(e){cat(sprintf("  !! rep%d nc%d %s FAILED\n",rp,N_comp,arm)); NULL})
    if (is.null(f)) next
    d <- f$draws
    rho_m <- sapply(1:K, function(k) mean(as.numeric(extract_variable_matrix(d, sprintf("rho[%d]",k)))))
    amp <- t(sapply(1:K, function(k){ P <- sapply(1:T_, function(t) as.numeric(extract_variable_matrix(d, sprintf("Phi[%d,%d]",t,k))))
      c(sd=mean(apply(P,1,sd)), rms=mean(apply(P,1,rms))) }))
    tau <- median(as.numeric(extract_variable_matrix(d,"tau")))
    sc <- mean(sapply(1:N_, function(n) mean(as.numeric(extract_variable_matrix(d, sprintf("scale_est[%d]",n))))))
    cat(sprintf("rep%d nc%d %-8s: rho = %s (mean %.3f) | realized sd = %s (mean %.3f -> mult %.2fx)\n",
      rp, N_comp, arm, paste(sprintf("%.3f",rho_m),collapse=","), mean(rho_m),
      paste(sprintf("%.2f",amp[,"sd"]),collapse=","), mean(amp[,"sd"]), 1/mean(amp[,"sd"])))
    cat(sprintf("            tau=%.3f | scale_est/os = %.2f | realized RMS mean=%.2f | div=%d rhat=%.2f\n",
      tau, sc/mean(os), mean(amp[,"rms"]), f$sampler_diag$n_div, f$sampler_diag$rhat_max))
  }
}
cat("\nRead: if posterior rho >> 0.90 the ex2 derived multiple is larger than 1.4x; the realized-sd column\n")
cat("gives the multiple DIRECTLY (1/realized sd), bypassing the rho -> amplitude step entirely.\n")
