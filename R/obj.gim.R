
obj.gim <- function(para, map, the0, inv.the, pi, inv.pi, ref){
  
  L <- (length(para) - 3) / 7
  n <- nrow(ref)
  
  bet <- para[map$bet]
  alp <- para[map$alp]
  the <- para[map$the]
  lam <- para[map$lam]
  
  ref <- as.matrix(ref)
  
  # h <- rcpp_hess_lambda(bet, a, alp, the, mu, gam, lam, ref, rho)
  
  # o <- -sum(log(1 + g %*% lam)) - .5 * t(the0-the) %*% inv.the %*% (the0-the) - .5 * t(pi-bet*alp) %*% inv.pi %*% (pi-bet*alp)
  o <- rcpp_obj(bet, alp, the, lam, the0, inv.the, pi, inv.pi, ref)
  
  o
  
}

