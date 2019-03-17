
## recover joint effect of case-control model
## modified from Bill's code

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
  
  ## the = (a, pi, mu)
  ##     true model: Delta = exp(a + X * pi)
  ## marginal model: delta_l = exp(mu_l + X_l * gam_l), l = 1, ..., L
  ## pi is joint effect of interest
  Sigma.pi <- fit$Sigma.the[2:(L+1), 2:(L+1)]
  se.pi <- sqrt(diag(Sigma.pi))
  
  rs <- outcome$iv
  names(fit$pi) <- rs
  names(se.pi) <- rs
  rownames(Sigma.pi) <- rs
  colnames(Sigma.pi) <- rs
  
  list(a = fit$a, pi = fit$pi, mu = fit$mu, 
       gam = gam, cov.gam = fit$Sigma.gam, 
       cov.pi = Sigma.pi, se.pi = se.pi, 
       g2 = fit$g2, H2 = fit$H2, 
       max.abs.diff=fit$max.abs.diff, max.rel.diff=fit$max.rel.diff, 
       score = fit$score, niter = fit$niter, converged=fit$converged)
  
}





