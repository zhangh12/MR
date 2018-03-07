
gfunction.cc <- function(ref, the, zlin, delta, delta1, dem1){

  n <- nrow(ref)
  np <- ncol(ref)
  g <- matrix(NA, n, 3*np+1)
  for(j in 1:np){
    g[, j] <- ref[, j] * (zlin - ref[, j] * the[j])
    g[, j+np] <- (delta - delta1[, j])/dem1[, j]
    g[, j+2*np] <- (delta - delta1[, j])/dem1[, j] * ref[, j]
  }
  g[, 3*np+1] <- delta-1

  g

}
