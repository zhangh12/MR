
alp.stage3 <- function(bet, alp, pi, Omega){
  
  L <- length(alp)
  if(nrow(Omega) != 2*L){
    msg <- 'debug NR.stage3'
    stop(msg)
  }
  
  Omega11 <- Omega[1:L, 1:L]
  Omega12 <- Omega[1:L, (L+1):(2*L)]
  Omega21 <- t(Omega12)
  Omega22 <- Omega[(L+1):(2*L), (L+1):(2*L)]
  
  mat1 <- Omega11 + bet * (Omega12 + Omega21) + bet^2 * Omega22
  mat2 <- (Omega11 + bet * Omega21) %*% alp + (Omega12 + bet * Omega22) %*% pi
  
  as.vector(solve(mat1) %*% mat2)
  
}
