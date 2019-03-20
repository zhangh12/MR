
gfunction <- function(para, map, ref){
  
  alp <- para[map$alp]
  the <- para[map$the]
  
  #lin <- as.vector(ref %*% alp)
  #ref * (lin - t(t(ref) * the))
  
  rcpp_gfunction(alp, the, ref)
  
}

