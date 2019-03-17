
## score function in stage 3

# alp - joint effect of exposure model
# pi - joint effect of case-control model (logOR)
dev.stage3 <- function(par, alp, pi, Omega11, Omega12, Omega21, Omega22){
  
  bet <- par[1]
  alp1 <- par[-1]
  
  d1 <- t(alp1) %*% Omega21 %*% (alp1 - alp) - t(alp1) %*% Omega22 %*% pi + bet * t(alp1) %*% Omega22 %*% alp1
  d2 <- (Omega11 + bet * (Omega12 + Omega21) + bet^2 * Omega22) %*% alp1 - ((Omega11 + bet * Omega21) %*% alp + (Omega12 + bet * Omega22) %*% pi)
  
  d1 <- as.vector(d1)
  d2 <- as.vector(d2)
  d <- c(d1, d2)
  
  d
  
}
