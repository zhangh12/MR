
obj.gim <- function(para, map, lam, the0, inv.the, pi, inv.pi, ref){
  
  bet <- para[map$bet]
  alp <- para[map$alp]
  the <- para[map$the]
  if(is.null(lam)){
    lam <- para[map$lam]
  }
  # h <- rcpp_hess_lambda(bet, a, alp, the, mu, gam, lam, ref, rho)
  
  # o <- -sum(log(1 + g %*% lam)) - .5 * t(the0-the) %*% inv.the %*% (the0-the) - .5 * t(pi-bet*alp) %*% inv.pi %*% (pi-bet*alp)
  rcpp_obj(bet, alp, the, lam, the0, inv.the, pi, inv.pi, ref)
  
}

