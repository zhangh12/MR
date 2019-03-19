
gfunction <- function(para, map, n1, n0, ref, pi){
  
  L <- (length(para) - 3) / 7
  n <- nrow(ref)
  
  bet <- para[map$bet]
  a <- para[map$a]
  alp <- para[map$alp]
  the <- para[map$the]
  mu <- para[map$mu]
  gam <- para[map$gam]
  lam <- para[map$lam]
  
  ref <- as.matrix(ref)
  rho <- diag(n1/n0)
  
  rcpp_gfunction(bet, a, alp, the, mu, gam, lam, ref, rho, pi)
  
}


gfunction.0 <- function(para, map, n1, n0, ref){
  
  L <- (length(para) - 3) / 7
  n <- nrow(ref)
  
  bet <- para[map$bet]
  a <- para[map$a]
  alp <- para[map$alp]
  the <- para[map$the]
  mu <- para[map$mu]
  gam <- para[map$gam]
  lam <- para[map$lam]
  
  ref <- as.matrix(ref)
  
  #tilt(para, map, ref)
  lin <- as.vector(ref %*% alp)
  Delta <- as.vector(exp(a + bet * lin))
  delta <- exp(t(t(ref) * gam + mu))
  
  if(0){ # do not delete, showing broadcasting in R
    tmp1 <- ref * lin
    tmp2 <- t(t(ref) * the) * ref
    tmp3 <- tmp1 - tmp2
  }
  
  g <- matrix(NA, nrow = n, ncol = 3 * L + 1)
  g[, 1:L] <- ref * (lin - t(t(ref) * the))
  
  rho <- diag(n1/n0)
  
  # Define PSI in 3 parts, columns 1:L, (L+1):2, and 2L+1
  onePlusRhoDelta <- 1 + t(t(delta) * rho)
  deltaOverOnePlusRhoDelta <- delta/onePlusRhoDelta
  PSI1 <- (Delta - delta)/onePlusRhoDelta
  PSI2 <- PSI1 * ref
  PSI3 <- Delta - 1
  
  g[, -c(1:L)] <- cbind(PSI1, PSI2, PSI3)
  
  g1 <- rcpp_gfunction(bet, a, alp, the, mu, gam, lam, ref, rho)
  
  g
  
}
