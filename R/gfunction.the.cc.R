

gfunction.the.cc <- function(ref){

  n <- nrow(ref)
  np <- ncol(ref)
  g.the <- list()
  for(k in 1:np){
    m <- matrix(0, n, 3*np+1)
    m[, k] <- -ref[, k]^2
    g.the[[k]] <- m
    rm(m)
  }

  g.the

}
