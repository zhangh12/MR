
NR <- function(para, map, the0, inv.the, gam0, inv.gam, 
               n1, n0, ref){
  
  f <- score.gim(para, map, the0, inv.the, gam0, inv.gam, 
                  n1, n0, ref)
  
  hess <- numDeriv::jacobian(score.gim, para, 
                             map = map, the0 = the0, inv.the = inv.the, 
                             gam0 = gam0, inv.gam = inv.gam, 
                             n1 = n1, n0 = n0, ref = ref)
  
  Jinv <- try(solve(hess), silent=FALSE)
  if (!("try-error" %in% class(Jinv))) {
    ret  <- para - Jinv %*% f
    ret  <- as.vector(ret)
  } else {
    # For now, add random small numbers to parms
    print(Jinv)
    ret <- as.vector(para) + 0.01*runif(length(para))
  } 
  
  ret
  
}
