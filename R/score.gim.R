
score.gim <- function(para, map, lam, the0, inv.the, pi, inv.pi, ref){
  
  bet <- para[map$bet]
  alp <- para[map$alp]
  the <- para[map$the]
  if(is.null(lam)){
    lam <- para[map$lam]
  }
  
  ref <- as.matrix(ref)
  
  sc <- rcpp_score(bet, alp, the, lam, the0, inv.the, pi, inv.pi, ref)
  # sc1 <- numDeriv::grad(obj.gim, para,
  #                      map = map, the0 = the0, inv.the = inv.the,
  #                      pi = pi, inv.pi = inv.pi, ref = ref)
  
  sc
  
}

