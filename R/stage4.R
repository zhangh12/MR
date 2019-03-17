
stage4 <- function(alp, cov.alp, pi, cov.pi){
  
  inv.alp <- solve(cov.alp)
  inv.pi <- solve(cov.pi)
  
  L <- length(alp)
  Omega <- matrix(0, 2*L, 2*L)
  Omega[1:L, 1:L] <- inv.alp
  Omega[(L+1):(2*L), (L+1):(2*L)] <- inv.pi
  
  par <- NR.stage3(alp, pi, Omega)
  cov3 <- stderr(par, alp, pi, Omega)
  
  fit <- list(exposure = list(alp = alp, cov.alp = cov.alp, 
                              alp2 = par[-1], cov.alp2 = cov3$cov.alp), 
              outcome = list(pi = pi, cov.pi = cov.pi), 
              Omega = Omega, 
              coefficients = par[1], 
              se = cov3$se)
  
  fit
  
}
