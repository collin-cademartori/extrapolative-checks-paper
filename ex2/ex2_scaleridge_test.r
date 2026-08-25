## Test the GLOBAL SCALE (product) ridge hypothesis. Exact likelihood invariance in unit-norm + est sigma:
##   innovations * c  ->  Phi * c  ->  signal * c ;  sigma_n / c  restores signal ;  tau * c restores error sd
## => a 1-D non-identified direction spanning ~90 innovations + 8 sigmas + tau, broken only by priors.
## A product ridge between a high-dim block and a scale is a FUNNEL whose neck is a LINEAR COMBINATION,
## which is why no single parameter / row-norm / factor / window localized the divergences.
## Predictions: corr(log innovScale, log sigma) < 0 ; corr(log tau, log sigma) < 0 ;
##              corr(log innovScale, log tau) > 0 ; divergences concentrate along that ridge coordinate,
##              separating BETTER than the full-parameter whitened LDA (0.61).
library(posterior)
SP <- "/private/tmp/claude-502/-Users-collin-Documents-Wake-Research-PriorPredictivePaper-reprod/0edcfb95-808f-4867-8f30-15ec1e2aec59/scratchpad"
z <- readRDS(file.path(SP,"div_subspace_fit.rds")); d <- z$d; sdg <- z$sdg; T_ <- z$T_; K <- z$K; M <- 8L
divflag <- as.logical(as.numeric(extract_variable_matrix(sdg,"divergent__"))); nd <- length(divflag)
V <- function(v) as.numeric(extract_variable_matrix(d, v))
zsh <- function(x) (mean(x[divflag])-mean(x[!divflag]))/sd(x)

# --- the three blocks of the invariance, on the log scale ---
Innov <- sapply(1:T_, function(t) sapply(1:K, function(k) NULL))  # placeholder
IM <- matrix(0, nd, T_*K); j <- 0
for(t in 1:T_) for(k in 1:K){ j <- j+1; IM[,j] <- V(sprintf("Phi_innovations[%d,%d]",t,k)) }
g_innov <- log(sqrt(rowMeans(IM^2)))                       # global innovation scale (log c)
PM <- matrix(0, nd, T_*K); j <- 0
for(t in 1:T_) for(k in 1:K){ j <- j+1; PM[,j] <- V(sprintf("Phi[%d,%d]",t,k)) }
g_phi <- log(sqrt(rowMeans(PM^2)))                         # realized factor amplitude (log)
SG <- sapply(1:M, function(n) V(sprintf("sigma[%d]",n)))
s_sig <- rowMeans(log(SG))                                 # global log sigma
l_tau <- log(V("tau"))

cat(sprintf("divergences: %d of %d (%.1f%%)\n\n", sum(divflag), nd, 100*mean(divflag)))
cat("=== (1) is the invariance visible as a posterior ridge? (predicted signs) ===\n")
cat(sprintf("  corr(log innovScale, log sigma) = %+.2f   [predict < 0]\n", cor(g_innov, s_sig)))
cat(sprintf("  corr(log PhiAmp   , log sigma) = %+.2f   [predict < 0]\n", cor(g_phi, s_sig)))
cat(sprintf("  corr(log tau      , log sigma) = %+.2f   [predict < 0]\n", cor(l_tau, s_sig)))
cat(sprintf("  corr(log innovScale, log tau ) = %+.2f   [predict > 0]\n", cor(g_innov, l_tau)))
cat(sprintf("  corr(log PhiAmp   , log tau ) = %+.2f   [predict > 0]\n", cor(g_phi, l_tau)))

cat("\n=== (2) do divergences localize along the ridge coordinate? (z = div vs bulk, in sd) ===\n")
for(nm in c("g_innov","g_phi","s_sig","l_tau")) cat(sprintf("  %-9s z=%+.2f | bulk=%.3f div=%.3f\n", nm,
  zsh(get(nm)), mean(get(nm)[!divflag]), mean(get(nm)[divflag])))
# the ridge coordinate itself: the invariance direction (innov up, sigma down, tau up), standardized
ridge <- scale(g_innov)[,1] - scale(s_sig)[,1] + scale(l_tau)[,1]
cat(sprintf("  RIDGE (g_innov - s_sig + l_tau): z=%+.2f\n", zsh(ridge)))
# and the orthogonal "overall size" combination for contrast
orth <- scale(g_innov)[,1] + scale(s_sig)[,1]
cat(sprintf("  contrast (g_innov + s_sig)     : z=%+.2f\n", zsh(orth)))

cat("\n=== (3) whitened separation of the 3-block scale summary vs full-param reference (0.61) ===\n")
sep <- function(X){ ctr <- scale(X); keep<-is.finite(colSums(ctr)); ctr<-ctr[,keep,drop=FALSE]
  dm <- colMeans(ctr[divflag,,drop=FALSE])-colMeans(ctr[!divflag,,drop=FALSE])
  S<-cor(ctr); lam<-0.05*mean(diag(S)); dir<-solve(S+lam*diag(ncol(S)),dm); dir<-dir/sqrt(sum(dir^2))
  p<-as.numeric(ctr%*%dir); (mean(p[divflag])-mean(p[!divflag]))/sd(p) }
cat(sprintf("  3 scale summaries (g_innov, s_sig, l_tau)      sep=%.2f\n", sep(cbind(g_innov,s_sig,l_tau))))
cat(sprintf("  + g_phi and per-unit log sigma (11 summaries)  sep=%.2f\n", sep(cbind(g_innov,g_phi,l_tau,log(SG)))))

cat("\n=== (4) funnel shape: is the CONDITIONAL spread of innovations tied to the scale? ===\n")
qs <- quantile(g_innov, seq(0,1,0.2)); bin <- cut(g_innov, qs, include.lowest=TRUE, labels=FALSE)
cat("  g_innov quintile | mean innov |sd| across coords | divergence RATE\n")
for(b in 1:5){ ix <- bin==b
  cat(sprintf("    Q%d  g=%+.2f   coordspread=%.3f   divrate=%.3f\n", b, mean(g_innov[ix]),
    mean(apply(IM[ix,,drop=FALSE],2,sd)), mean(divflag[ix]))) }
cat("\nRead: predicted signs + a ridge z / sep clearly beating the marginals (and a monotone divergence rate across g_innov quintiles) => global product-scale funnel.\n")
