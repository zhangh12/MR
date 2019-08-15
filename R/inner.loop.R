
inner.loop <- function(para, map, rho, ref, pi){
  
  #message('Inner loop...')
  
  L <- (length(para) - 3) / 7
  n <- nrow(ref)
  
  lam0 <- rep(0, 3 * L + 1)
  
  g <- gfunction1(para, map, rho, ref, pi)
  
  lam0 <- lam0[-c(1:2)]
  g <- g[, -c(1:2)]
  
  M <- chol(t(g) %*% g)
  g <- g %*% solve(M)
  lam0 <- as.vector(M %*% lam0)
  
  foo <- function(lam, g, deriv){
    n <- nrow(g)
    pr <- 1/as.vector(1 + g %*% lam)
    
    if(deriv == 0){
      return( sum(log(pr)) )
    }
    
    if(deriv == 1){
      return( -as.vector(t(g) %*% pr) )
    }
    
    if(deriv == 2){
      return( t(g) %*% (g * pr^2) )
    }
  }
  
  niter <- 100
  tol <- 1e-6
  fn <- foo(lam0, g, 0)
  for(i in seq(niter)){
    
    s <- foo(lam0, g, 1)
    h <- foo(lam0, g, 2)
    
    min.ev <- min(eigen(h)$value)
    
    if(max(abs(s)) < tol){
      message('inner loop converges')
      break
    }
    
    d <- solve(h, s)
    
    for(j in 1:10){
      tau <- .5^(j-1)
      lam1 <- lam0 - tau * d
      if(all(1 + g %*% lam1 > 1/n)){
        break
      }
    }
    
    lam0 <- lam1
    fn0 <- foo(lam0, g, 0)
    fn <- c(fn, fn0)
    cat('\tgrad=', max(abs(s)), '\tstep=', tau, '\tmin.ev=', min.ev, '\tfn=', fn0, '\n')
    plot(fn, pch=20, main='Inner Loop')
  }
  
  if(min.ev < 0){
    warning(paste0('Inner loop may stop at a saddle point: min ev = ', min.ev))
  }
  
  if(max(abs(s)) > tol){
    warning(paste0('Inner loop does not stop at a stationary points: max grad = ', max(abs(s))))
  }
  
  print(s)
  print(lam0)
  
  lam0
  
}
