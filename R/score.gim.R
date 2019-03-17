
score.gim <- function(para, map, the0, inv.the, gam0, inv.gam, 
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
  
  d.the <- as.vector(inv.the %*% (the - the0))
  d.gam <- as.vector(inv.gam %*% (gam - gam0))
  rho <- diag(n1/n0)
  
  sc <- rcpp_score(bet, a, alp, the, mu, gam, lam, ref, rho, d.the, d.gam)
  # sc <- numDeriv::grad(obj.gim, para, 
  #                      map = map, the0 = the0, inv.the = inv.the, 
  #                      gam0 = gam0, inv.gam = inv.gam, 
  #                      n1 = n1, n0 = n0, ref = ref)
  
  sc
  
}

