
## compute SE of causal effect

# alp - joint effect of exposure model
# pi - joint effect of case-control model (logOR)
stderr <- function(bet, alp1, alp, pi, Omega){
  
  L <- length(alp)
  
  Omega11 <- Omega[1:L, 1:L]
  Omega12 <- Omega[1:L, (L+1):(2*L)]
  Omega21 <- t(Omega12)
  Omega22 <- Omega[(L+1):(2*L), (L+1):(2*L)]
  
  hess <- matrix(0, L+1, L+1)
  hess[1, 1] <- t(alp1) %*% Omega22 %*% alp1
  hess[-1, 1] <- Omega12 %*% alp1 + bet * Omega22 %*% alp1
  hess[1, -1] <- t(hess[-1, 1])
  hess[-1, -1] <- Omega11 + bet * (Omega12 + Omega21) + bet^2 * Omega22
  
  cov <- solve(hess)
  rownames(cov) <- names(par)
  colnames(cov) <- names(par)
  se <- sqrt(cov[1,1])
  
  list(cov.alp = cov[-1,-1], se = se)
  
}




stderr0 <- function(par, alp, pi, Omega){
  
  L <- length(alp)
  
  Omega11 <- Omega[1:L, 1:L]
  Omega12 <- Omega[1:L, (L+1):(2*L)]
  Omega21 <- t(Omega12)
  Omega22 <- Omega[(L+1):(2*L), (L+1):(2*L)]
  
  bet <- par[1]
  alp1 <- par[-1]
  
  info <- matrix(0, L+1, L+1)
  info[1, 1] <- t(alp1) %*% Omega22 %*% alp1
  info[-1, 1] <- (Omega12 + bet * Omega22) %*% alp1
  info[1, -1] <- t(info[-1, 1])
  info[-1, -1] <- Omega11 + bet * (Omega12 + Omega21) + bet^2 * Omega22
  
  hess <- matrix(0, L+1, L+1)
  hess[1, 1] <- info[1, 1]
  hess[-1, 1] <- (Omega12 + Omega21 + 2 * bet * Omega22) %*% alp1 - (Omega21 %*% alp + Omega22 %*% pi)
  hess[1, -1] <- t(hess[-1, 1])
  hess[-1, -1] <- info[-1, -1]
  
  cov2 <- solve(hess)
  
  cov <- solve(hess) %*% info %*% solve(hess)
  rownames(cov) <- names(par)
  colnames(cov) <- names(par)
  se <- sqrt(cov[1,1])
  
  list(cov.alp = cov[-1,-1], se = se)
  
}

