
## object function in stage 3

# alp - joint effect of exposure model
# inv.alp - precise matrix of alp
# pi - joint effect of case-control model (logOR)
# inv.pi - precise matrix of pi
obj.stage3 <- function(par, alp, inv.alp, pi, inv.pi){
  
  bet <- par[1]
  alp1 <- par[-1]
  
  t(alp - alp1) %*% inv.alp %*% (alp - alp1) + t(pi - bet*alp1) %*% inv.pi %*% (pi - bet*alp1)
  
}
