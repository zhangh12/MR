
## LM tests for H0: bet = bet0
# set level = 0 so that fac = 0, i.e. returning score statistics
lmt.stat <- function(bet0, alp, pi, V, Omega11, Omega12, Omega21, Omega22, level){
  
  mat1 <- solve(Omega11 + bet0 * (Omega12 + Omega21) + bet0^2 * Omega22)
  G1 <- mat1 %*% (Omega11 + bet0 * Omega21)
  G2 <- mat1 %*% (Omega12 + bet0 * Omega22)
  eta <- G1 %*% alp + G2 %*% pi
  sc <- as.vector(t(eta) %*% Omega21 %*% (eta - alp) - t(eta) %*% Omega22 %*% pi + bet0 * t(eta) %*% Omega22 %*% eta)
  
  H0 <- t(eta - alp) %*% Omega12 + t(eta) %*% Omega21 - t(pi) %*% Omega22 + 2 * bet0 * t(eta) %*% Omega22
  H1 <- -t(eta) %*% Omega21
  H2 <- -t(eta) %*% Omega22
  L <- length(alp)
  mat2 <- matrix(NA, nrow = 1, ncol = 2*L)
  mat2[1, 1:L] <- H0 %*% G1 + H1
  mat2[1, (L+1):(2*L)] <- H0 %*% G2 + H2
  
  v <- as.vector(mat2 %*% V %*% t(mat2))
  
  fac <- qchisq(level, df = 1)
  as.vector(sc^2/v) - fac
  
}
