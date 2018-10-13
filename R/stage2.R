
## recover joint effect of case-control model

# gam - gamma, marginal logOR, nominal SE of gam is not needed
# n1 - matrix of sharing sample sizes of cases among gam
# n0 - matrix of sharing sample sizes of controls among gam
# ref - reference panel
stage2 <- function(outcome, n1, n0, ref){
  
  gam <- outcome$beta
  
  L <- length(gam) # number of markers
  N <- nrow(ref) # number of reference samples
  X <- as.matrix(ref)
  rm(ref)
  
  parms <- rep(0, 2*L+1)
  
  fit <- runNR(parms, X, n1, n0, gam)
  
  Sigma.pi <- fit$Sigma.the[2:(L+1), 2:(L+1)]
  se.pi <- sqrt(diag(Sigma.pi))
  
  rs <- outcome$iv
  names(fit$pi) <- rs
  names(se.pi) <- rs
  rownames(Sigma.pi) <- rs
  colnames(Sigma.pi) <- rs
  
  list(a = fit$a, pi = fit$pi, mu = fit$mu, gam = gam, 
       cov.pi = Sigma.pi, se.pi = se.pi, 
       max.abs.diff=fit$max.abs.diff, max.rel.diff=fit$max.rel.diff, 
       score = fit$score, niter = fit$niter, converged=fit$converged)
  
}





