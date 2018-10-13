
## score function in stage 3

# alp - joint effect of exposure model
# inv.alp - precise matrix of alp
# pi - joint effect of case-control model (logOR)
# inv.pi - precise matrix of pi
dev.stage3 <- function(par, alp, inv.alp, pi, inv.pi){
  
  bet <- par[1]
  alp1 <- par[-1]
  
  d1 <- t(alp1) %*% inv.pi %*% (pi - bet*alp1)
  d2 <- inv.alp %*% (alp - alp1) + bet * inv.pi %*% (pi - bet * alp1)
  d2 <- as.vector(d2)
  d <- -c(d1, d2)
  
  d
  
}
