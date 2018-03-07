

gfunction.alp.cc <- function(ref, bet, delta, dem1){

  n <- nrow(ref)
  np <- ncol(ref)
  g.alp <- list()
  for(k in 1:np){
    m <- matrix(0, n, 3*np+1)
    m[, 1:np] <- ref * ref[, k]
    for(j in 1:np){
      m[, j+np] <- bet * delta/dem1[, j] * ref[, k]
      m[, j+2*np] <- m[, j+np] * ref[, j]
    }
    m[, 3*np+1] <- bet * delta * ref[, k]
    g.alp[[k]] <- m
    rm(m)
  }

  g.alp

}
