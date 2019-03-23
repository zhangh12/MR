
NR.lambda <- function(para, map, ref){
  
  #message('Inner loop...')
  
  bet <- para[map$bet]
  alp <- para[map$alp]
  the <- para[map$the]
  
  L <- length(alp)
  n <- nrow(ref)
  
  lam0 <- rep(0, L)
  
  g <- gfunction(para, map, ref)
  # M <- chol(t(g) %*% g)
  # g.star <- g %*% solve(M)
  # t(g.star) %*% g.star # close to a identity matrix
  # lam.star <- as.vector(M %*% lam)
  
  g.star <- g
  foo <- function(lam, g, deriv){
    n <- nrow(g)
    m <- ncol(g)
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
  fn <- foo(lam0, g.star, 0)
  for(i in seq(niter)){
    
    s <- foo(lam0, g.star, 1)
    h <- foo(lam0, g.star, 2)
    if(max(abs(s)) < tol){
      #print('converged')
      break
    }

    d <- solve(h, s)
    
    for(j in 1:10){
      tau <- .5^(j-1)
      lam1 <- lam0 - tau * d
      tmp <- as.vector(1 + g.star %*% lam1)
      if(all(tmp > 1/n)){
        break
      }
    }
    
    lam0 <- lam1
    fn <- c(fn, foo(lam0, g.star, 0))
    #plot(fn, pch=20, main='Inner Loop')
  }
  
  min.ev <- min(eigen(h)$value)
  if(min.ev < 0){
    warning(paste0('Inner loop may stop at a saddle point: min ev = ', min.ev))
  }
  
  if(max(abs(s)) > tol){
    warning(paste0('Inner loop does not stop at a stationary points: max grad = ', max(abs(s))))
  }
  
  lam0
  
}
