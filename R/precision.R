
# compute covariance of alp and the = (a, pi, mu)
# then return precision matrix of alp and pi

precision <- function(fit1, fit2){
  
  g1 <- fit1$g1
  H1 <- fit1$H1
  g2 <- fit2$g2
  H2 <- fit2$H2
  
  if(nrow(g1) != nrow(g2)){
    msg <- 'debug cov.joint.effect (err 1)'
    stop(msg)
  }
  
  if(ncol(g2) != ncol(g1) * 2 + 1){
    msg <- 'debug cov.joint.effect (err 2)'
    stop(msg)
  }
  
  N <- nrow(g1)
  L <- ncol(g1)
  c <- solve(H1) %*% (cov(g1, g2)/N) %*% solve(H2)
  
  V <- c[, 2:(L+1)]
  #V[] <- 0
  
  Omega <- matrix(NA, 2*L, 2*L)
  Omega[1:L, 1:L] <- fit1$cov.alp
  Omega[(L+1):(2*L), (L+1):(2*L)] <- fit2$cov.pi
  Omega[1:L, (L+1):(2*L)] <- V
  Omega[(L+1):(2*L), 1:L] <- t(V)
  
  Omega <- solve(Omega)
  
  Omega
  
}
