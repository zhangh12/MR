
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
  
  g <- gfunction(para, map, n1, n0, ref)
  
  o <- sum(log(1 + g %*% lam)) + 1/2 * as.vector(t(the0 - the) %*% inv.the %*% (the0 - the)) + 
    1/2 * as.vector(t(gam0 - gam) %*% inv.gam %*% (gam0 - gam))
  
  o
  
}

