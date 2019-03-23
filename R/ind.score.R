

ind.score <- function(para, map, the0, inv.the, pi, inv.pi, ref){
  
  bet <- para[map$bet]
  alp <- para[map$alp]
  the <- para[map$the]
  lam <- para[map$lam]
  
  rcpp_ind_score(bet, alp, the, lam, the0, inv.the, pi, inv.pi, ref)
  
}
