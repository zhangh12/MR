
NR.gim <- function(para, map, the0, inv.the, pi, inv.pi, ref){
  
  L <- (length(para) - 3) / 7
  n <- nrow(ref)
  
  bet <- para[map$bet]
  alp <- para[map$alp]
  the <- para[map$the]
  lam <- para[map$lam]
  
  x0 <- para
  fn <- obj.gim(x0, map, the0, inv.the, pi, inv.pi, ref)
  niter <- 20
  tol <- 1e-6
  for(i in seq(niter)){
    s <- score.gim(x0, map, the0, inv.the, pi, inv.pi, ref)
    
    cat('bet = ', x0[map$bet], ',\ts = ', max(abs(s)), '\n')
    h <- numDeriv::jacobian(score.gim, x0, 
                            map = map, the0 = the0, inv.the = inv.the, 
                            pi = pi, inv.pi = inv.pi, ref = ref)
    
    if(max(abs(s)) < tol){
      print('converged')
      break
    }
    
    d <- solve(h, s)
    
    for(j in 1:10){
      tau <- .9^(j-1)
      x1 <- x0 - tau * d
      g <- gfunction(x1, map, ref)
      tmp <- as.vector(1 + g %*% x1[map$lam])
      fn1 <- obj.gim(x1, map, the0, inv.the, pi, inv.pi, ref)
      if(all(tmp > 1/n) && fn1 > tail(fn, 1)){
        break
      }
    }
    
    x0 <- x1
    fn <- c(fn, obj.gim(x0, map, the0, inv.the, pi, inv.pi, ref))
    plot(fn, pch=20)
  }
  
  s <- score.gim(x0, map, the0, inv.the, pi, inv.pi, ref)
  h <- numDeriv::jacobian(score.gim, x0, 
                          map = map, the0 = the0, inv.the = inv.the, 
                          pi = pi, inv.pi = inv.pi, ref = ref)
  
  max.ev <- max(eigen(h+t(h))$value)
  
  cat('max ev = ', max.ev, ', s = ', max(abs(s)), '\n')
  
  if(max(abs(s)) > tol){
    if(i == niter){
      warning('Iteration limit had been reached')
    }else{
      warning('Procedure does not stop at a stationary point')
    }
    bet <- NA
    se <- NA
  }else{
    if(max.ev > 0){
      warning('Procedure may have stopped at a saddle point')
    }
    bet <- x0[map$bet]
    se <- sqrt(-solve(h)[map$bet, map$bet])
  }
  
  list(coefficients = bet, se = se)
  
}
