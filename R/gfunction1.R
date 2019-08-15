
gfunction1 <- function(para, map, rho, ref, pi){
  
  bet <- para[map$bet]
  a <- para[map$a]
  alp <- para[map$alp]
  the <- para[map$the]
  mu <- para[map$mu]
  gam <- para[map$gam]
  
  #lin <- as.vector(ref %*% alp)
  #ref * (lin - t(t(ref) * the))
  
  g <- rcpp_gfunction1(bet, a, alp, the, mu, gam, rho, ref, pi)
  g
  
}

