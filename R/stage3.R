

stage3 <- function(alp, cov.alp, pi, cov.pi){
  
  inv.alp <- solve(cov.alp)
  inv.pi <- solve(cov.pi)
  
  par <- NR.stage3(c(0, alp), alp, inv.alp, pi, inv.pi)
  cov3 <- stderr(par, alp, inv.alp, pi, inv.pi)
  
  fit <- list(exposure = list(alp = alp, cov.alp = cov.alp, 
                              alp2 = par[-1], cov.alp2 = cov3$cov.alp), 
              outcome = list(pi = pi, cov.pi = cov.pi), 
              coefficients = c(Causal=par[1]), 
              se = cov3$se)
  
  fit
  
}
