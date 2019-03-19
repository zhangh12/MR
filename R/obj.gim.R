
obj.gim <- function(para, map, the0, inv.the, gam0, inv.gam, 
                    n1, n0, ref){
  
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
  
  # g <- gfunction(para, map, n1, n0, ref)
  # o1 <- -sum(log(1 + g %*% lam)) - .5 * t(the0-the) %*% inv.the %*% (the0-the) - .5 * t(gam0-gam) %*% inv.gam %*% (gam0-gam)
  
  # h <- rcpp_hess_lambda(bet, a, alp, the, mu, gam, lam, ref, rho)
  o <- rcpp_obj(bet, a, alp, the, mu, gam, lam, ref, rho, the0, inv.the, gam0, inv.gam)
  
  o
  
}

