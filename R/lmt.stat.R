
## LM tests for H0: bet = bet0
# set level = 0 so that fac = 0, i.e. returning score statistics
lmt.stat <- function(bet0, alp, inv.alp, pi, inv.pi, level){
  
  eta <- solve(inv.alp + bet0^2 * inv.pi) %*% (inv.alp %*% alp + bet0 * inv.pi %*% pi)
  sc <- t(eta) %*% inv.pi %*% (bet0 * eta - pi)
  v <- -t(2*bet0 * eta-pi) %*% inv.pi %*% solve(inv.alp + bet0^2 * inv.pi) %*% inv.pi %*% pi + 
    t(eta) %*% inv.pi %*% eta
  
  fac <- qchisq(level, df = 1)
  as.vector(sc^2/v) - fac
  
}
