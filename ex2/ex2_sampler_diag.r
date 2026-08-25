## Sampler deep-dive (NOT shipped): trap-prone est-sigma no_ints fit, normal(0,1).
## Q3: does pathfinder seed Lambda_raw? A/B: PF_PARAM_BASES with "Lambda_raw" (correct) vs "Lambda"
##     (buggy, leaves loadings at random init). Inspect the actual inits + missing-init, compare
##     offmode/div/rhat. Q2: for the correct arm, localize modes -- which params differ most between
##     off-mode and on-mode chains. Q1: localize divergences -- which params are extreme at divergent
##     transitions. Single-process (no cluster env bug); 8 chains run in parallel.
library(posterior)
source("../sample_model.r"); source("../pathfinder_init.r")
ife_mod <- cmdstan_model("../ife_named_unitnorm_ig_s1.stan")   # normal(0,1)

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
stat_data <- list(M_units=N_, T_times=T_, K_latent=3L, Y=fy, a_rho=90, b_rho=10, tau_val=0, m_tau=0.1, s_tau=0.05,
  sigma_data=apply(fy,2,function(x)sqrt(mean(x^2))), fit_overall_scales=1L, nonstationary=0L, unit_intercepts=0L,
  factor_means=1L, sample_posterior=1L, num_treated=5L, gamma_scale=3, gamma_loc=4, alpha_diag=10)

run_arm <- function(pf_bases) {
  assign("PF_PARAM_BASES", pf_bases, envir=globalenv())
  inits <- pathfinder_inits(ife_mod, stat_data, 8, seed=42, quiet=TRUE)
  init_names <- if(is.null(inits)) "PF_FAILED" else paste(sort(names(inits[[1]])), collapse=",")
  has_Lraw <- grepl("Lambda_raw", init_names); has_L <- grepl("(^|,)Lambda(,|$)", init_names)
  f <- ife_mod$sample(data=stat_data, chains=8, parallel_chains=8, iter_warmup=1000, iter_sampling=800,
    adapt_delta=0.9, max_treedepth=12, init=if(is.null(inits)) 2 else inits, seed=42, refresh=0, show_messages=FALSE, show_exceptions=FALSE)
  list(f=f, init_names=init_names, has_Lraw=has_Lraw, has_L=has_L)
}

summ <- function(f) {
  d <- f$draws()
  lp <- colMeans(extract_variable_matrix(d, "lp__")); gap <- max(lp)-lp; on <- which(gap<=5)
  idv <- grep("^(sigma|tau|delta|Lambda|rho|Phi_means|omega)", dimnames(d)$variable, value=TRUE)
  idv <- idv[!grepl("Lambda_raw", idv)]  # identified params
  ra <- suppressWarnings(max(summarise_draws(subset_draws(d, variable=idv), "rhat")$rhat, na.rm=TRUE))
  ron <- if(length(on)>=2) suppressWarnings(max(summarise_draws(subset_draws(d, variable=idv, chain=on), "rhat")$rhat, na.rm=TRUE)) else NA
  dv <- f$diagnostic_summary(quiet=TRUE)
  list(lp=lp, on=on, off=setdiff(seq_along(lp), on), rhat_all=ra, rhat_on=ron, ndiv=sum(dv$num_divergent), d=d)
}

cat("=== Q3: pathfinder init A/B (same dataset/seed) ===\n")
for(nm in c("correct(Lambda_raw)","buggy(Lambda)")) {
  bases <- if(grepl("raw", nm)) c("sigma_raw","tau_param","omega_sq_param","Phi_innovations","Phi_means_param","rho","Lambda_raw","gamma_raw","delta_raw") else c("sigma_raw","tau_param","omega_sq_param","Phi_innovations","Phi_means_param","rho","Lambda","gamma_raw","delta_raw")
  a <- run_arm(bases); s <- summ(a$f)
  cat(sprintf("%-20s: init has Lambda_raw=%s Lambda=%s | offmode=%d/8 rhat_all=%.2f rhat_on=%.2f div=%d\n",
    nm, a$has_Lraw, a$has_L, length(s$off), s$rhat_all, s$rhat_on, s$ndiv))
  cat("   init params:", a$init_names, "\n")
  if(grepl("raw", nm)) COR <- s  # keep the correct arm for deep-dive
}

cat("\n=== Q2: localize modes (correct arm) -- params most different between off-mode and on-mode chains ===\n")
if(length(COR$off) >= 1 && length(COR$on) >= 1) {
  d <- COR$d; vars <- grep("^(sigma_raw|tau_param|Lambda_raw|delta_raw|rho|Phi_means_param|omega_sq_param|Phi_innovations|Lambda\\[|sigma\\[)", dimnames(d)$variable, value=TRUE)
  cm <- sapply(vars, function(v){m<-extract_variable_matrix(d,v); colMeans(m)})  # [chain x var]
  onm <- colMeans(cm[COR$on,,drop=FALSE]); offm <- colMeans(cm[COR$off,,drop=FALSE])
  pooled_sd <- apply(cm, 2, sd); pooled_sd[pooled_sd<1e-8] <- NA
  zdiff <- abs(offm - onm)/pooled_sd
  ord <- order(-zdiff); top <- head(ord, 15)
  cat(sprintf("on-mode chains: %s | off-mode chains: %s\n", paste(COR$on,collapse=","), paste(COR$off,collapse=",")))
  for(j in top) cat(sprintf("  %-18s zdiff=%.1f | on=%.2f off=%.2f\n", vars[j], zdiff[j], onm[j], offm[j]))
} else cat("  (no off-mode chains in correct arm -- seeding fixed it!)\n")

cat("\n=== Q1: localize divergences (correct arm) -- params most extreme at divergent transitions ===\n")
dv <- extract_variable_matrix(COR$f$sampler_diagnostics(), "divergent__")
d <- COR$d; vars <- grep("^(sigma_raw|tau_param|Lambda_raw|delta_raw|rho|Phi_means_param|omega_sq_param)", dimnames(d)$variable, value=TRUE)
divflag <- as.logical(dv)
if(sum(divflag) >= 3) {
  zex <- sapply(vars, function(v){x<-as.numeric(extract_variable_matrix(d,v)); (mean(x[divflag])-mean(x[!divflag]))/sd(x)})
  ord <- order(-abs(zex)); for(j in head(ord,15)) cat(sprintf("  %-18s z=%.2f (div mean vs bulk)\n", vars[j], zex[j]))
  cat(sprintf("(total divergent draws: %d of %d)\n", sum(divflag), length(divflag)))
} else cat(sprintf("  (only %d divergent draws -- too few to localize)\n", sum(divflag)))
