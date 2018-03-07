
gfunction.bet.cc <- function(ref, zlin, delta, dem1){

  n <- nrow(ref)
  np <- ncol(ref)
  g.bet <- matrix(0, n, 3*np+1)
  for(j in 1:np){
    g.bet[, j+np] <- delta/dem1[, j] * zlin
    g.bet[, j+2*np] <- g.bet[, j+np] * ref[, j]
  }
  g.bet[, 3*np+1] <- delta * zlin

  g.bet

}
