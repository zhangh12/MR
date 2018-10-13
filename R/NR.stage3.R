## object function in stage 3

# alp - joint effect of exposure model
# inv.alp - precise matrix of alp
# pi - joint effect of case-control model (logOR)
# inv.pi - precise matrix of pi
NR.stage3 <- function(par, alp, inv.alp, pi, inv.pi){
  
  bet0 <- par[1]
  alp0 <- par[-1]
  
  k <- 500
  conv <- FALSE
  while(k > 0){
    tmp <- inv.pi %*% pi
    bet1 <- as.vector(1/(t(alp0) %*% inv.pi %*% alp0) * t(alp0) %*% tmp)
    alp1 <- as.vector(solve(inv.alp + bet1^2 * inv.pi) %*% (inv.alp %*% alp + bet1 * tmp))
    bet0 <- bet1
    alp0 <- alp1
    k <- k-1
    sc <- dev.stage3(c(bet0, alp0), alp, inv.alp, pi, inv.pi)
    if(max(abs(sc)) < 1e-9){
      conv <- TRUE
      break
    }
  }
  
  if(!conv){
    msg <- 'stage 3 may not converge'
    print(sc)
    stop(msg)
  }
  
  par <- c(bet0, alp0)
  names(par) <- c('Causal', names(alp))
  par
  
}


