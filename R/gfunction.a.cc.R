
gfunction.a.cc <- function(ref, delta, dem1){

  n <- nrow(ref)
  np <- ncol(ref)
  g.a <- matrix(0, n, 3*np+1)
  for(j in 1:np){
    g.a[, j+np] <- delta/dem1[, j]
    g.a[, j+2*np] <- g.a[, j+np] * ref[, j]
  }
  g.a[, 3*np+1] <- delta

  g.a

}
