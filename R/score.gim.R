
score.gim <- function(para, map, the0, inv.the, pi, inv.pi, ref){
  
  L <- (length(para) - 3) / 7
  n <- nrow(ref)
  
  bet <- para[map$bet]
  alp <- para[map$alp]
  the <- para[map$the]
  lam <- para[map$lam]
  
  ref <- as.matrix(ref)
  
  sc <- rcpp_score(bet, alp, the, lam, the0, inv.the, pi, inv.pi, ref)
  # sc1 <- numDeriv::grad(obj.gim, para,
  #                      map = map, the0 = the0, inv.the = inv.the,
  #                      pi = pi, inv.pi = inv.pi, ref = ref)
  
  sc
  
}

