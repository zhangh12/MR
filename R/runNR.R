
## Modified from Wheeler Bill's code

# X - reference panel
# n1Mat - matrix of sharing sample sizes of cases among gam
# n0Mat - matrix of sharing sample sizes of controls among gam
# gamVec - gamma, marginal logOR, nominal SE of gam is not needed
runNR <- function(parms, X, n1Mat, n0Mat, gamVec, op=NULL) {
  
  op <- default.list(op, 
                     c("tol.rel", "tol.diff", "niter", "eps", "print"), 
                     list(1e-4, 1e-4, 100, 1e-8, 0))
  tol.rel  <- op$tol.rel
  tol.diff <- op$tol.diff 
  niter    <- op$niter
  eps      <- op$eps
  print    <- op$print
  conv     <- FALSE
  
  stop <- 0
  x0   <- parms
  for (i in 1:niter) {
    x1 <- NewtonRaphson(x0, X, n1Mat, n0Mat, gamVec) 
    v1 <- abs(x1 - x0)
    v2 <- v1/(abs(x0) + eps)
    d1 <- max(v1)
    d2 <- max(v2)
    if (print) {
      print(paste("Iteration: ", i, ", max abs diff=", d1, ", max rel diff=", d2, sep=""))
    }
    
    psi <- PSI0(x1, X, n1Mat, n0Mat, gamVec, computeJ = FALSE)
    
    # if ((d1 < tol.diff) || (d2 < tol.rel)) {
    if (max(abs(psi)) < 1e-9) {
      conv <- TRUE
      break
    }
    x0 <- x1
  }
  
  dev <- PSI0(x1, X, n1Mat, n0Mat, gamVec, computeSigma.theta = TRUE)
  if(max(abs(dev$PSI)) > 1e-9){
    msg <- 'stage 2 may not converge'
    print(dev$PSI)
    stop(msg)
  }
  
  ##     true model: Delta = exp(a + X * pi)
  ## marginal model: delta_l = exp(mu_l + X_l * gam_l), l = 1, ..., L
  ## pi is joint effect of interest
  L <- length(gamVec)
  a <- x1[1]
  pi <- x1[2:(L+1)]
  mu <- x1[(L+2):(2*L+1)]
  gam <- gamVec
  
  list(parms=x1, score = dev$PSI, 
       a = a, pi = pi, mu = mu, gam = gam, 
       Sigma.gam = dev$Sigma.gam, Sigma.the = dev$Sigma.the, 
       g2 = dev$g2, H2 = dev$H2, 
       niter=i, max.abs.diff=d1, max.rel.diff=d2, converged=conv)
  
} # END: runNR


