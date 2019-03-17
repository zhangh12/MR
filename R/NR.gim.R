
NR.gim <- function(para, map, the0, inv.the, gam0, inv.gam, 
                   n1, n0, ref, op = NULL){
  
  op <- default.list(op, 
                     c("tol.rel", "tol.diff", "niter", "eps", "print"), 
                     list(1e-4, 1e-4, 100, 1e-8, 1))
  tol.rel  <- op$tol.rel
  tol.diff <- op$tol.diff 
  niter    <- op$niter
  eps      <- op$eps
  print    <- op$print
  conv     <- FALSE
  
  stop <- 0
  x0   <- para
  for (i in 1:niter) {
    x1 <- NR(x0, map, the0, inv.the, gam0, inv.gam, 
             n1, n0, ref) 
    v1 <- abs(x1 - x0)
    v2 <- v1/(abs(x0) + eps)
    d1 <- max(v1)
    d2 <- max(v2)
    o <- obj.gim(x1, map, the0, inv.the, gam0, inv.gam, n1, n0, ref)
    if (print) {
      print(paste("Iteration: ", i, ", max abs diff=", d1, ", max rel diff=", d2, ", obj=", o, sep=""))
    }
    
    sc <- score.gim(x1, map, the0, inv.the, gam0, inv.gam, 
                    n1, n0, ref)
    
    # if ((d1 < tol.diff) || (d2 < tol.rel)) {
    if (max(abs(sc)) < 1e-9) {
      conv <- TRUE
      break
    }
    x0 <- x1
  }
  
  sc <- score.gim(x1, map, the0, inv.the, gam0, inv.gam, 
                  n1, n0, ref)
  if(max(abs(sc)) > 1e-9){
    msg <- 'NR.gim may not converge'
    print(sc)
    stop(msg)
  }
  
  list(para = x1, score = sc, 
       niter=i, max.abs.diff=d1, max.rel.diff=d2, converged=conv)
  
}
