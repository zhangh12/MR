

gfunction.gam.cc <- function(ref, delta1, dem, dem1){

  n <- nrow(ref)
  np <- ncol(ref)
  g.gam <- list()
  for(k in 1:np){
    m <- matrix(0, n, 3*np+1)
    m[, k+np] <- -dem[, k]/dem1[, k]^2 * delta1[, k] * ref[, k]
    m[, k+2*np] <- m[, k+np] * ref[, k]
    g.gam[[k]] <- m
    rm(m)
  }

  g.gam

}
