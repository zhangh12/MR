
# alp - joint effect of exposure model
# pi - joint effect of case-control model (logOR)
NR.stage3 <- function(alp, pi, Omega){
  
  L <- length(alp)
  if(nrow(Omega) != 2*L){
    msg <- 'debug NR.stage3'
    stop(msg)
  }
  
  Omega11 <- Omega[1:L, 1:L]
  Omega12 <- Omega[1:L, (L+1):(2*L)]
  Omega21 <- t(Omega12)
  Omega22 <- Omega[(L+1):(2*L), (L+1):(2*L)]
  
  bet0 <- 0
  alp0 <- alp
  
  k <- 20000
  conv <- FALSE
  while(k > 0){
    bet1 <- (t(alp0) %*% Omega22 %*% pi + t(alp0) %*% Omega21 %*% (alp - alp0)) / (t(alp0) %*% Omega22 %*% alp0)
    bet1 <- as.vector(bet1)
    
    alp1 <- solve(Omega11 + bet1 * (Omega12 + Omega21) + bet1^2 * Omega22) %*% ((Omega11 + bet1 * Omega21) %*% alp + (Omega12 + bet1 * Omega22) %*% pi)
    alp1 <- as.vector(alp1)
    
    bet0 <- bet1
    alp0 <- alp1
    k <- k-1
    
    sc <- dev.stage3(c(bet0, alp0), alp, pi, Omega11, Omega12, Omega21, Omega22)
    
    if(max(abs(sc)) < 1e-5){
      conv <- TRUE
      break
    }
  }
  
  if(!conv){
    msg <- 'stage 3 may not converge'
    print(sc)
    stop(msg)
  }
  
  par <- c(bet0, alp0)
  names(par) <- c('Causal', names(alp))
  par
  
}


