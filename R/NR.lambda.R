
NR.lambda <- function(para, map, ref){
  
  L <- (length(para) - 3) / 7
  n <- nrow(ref)
  
  bet <- para[map$bet]
  alp <- para[map$alp]
  the <- para[map$the]
  lam <- para[map$lam]
  
  g <- gfunction(alp, the, ref)
  # M <- chol(t(g) %*% g)
  # g.star <- g %*% solve(M)
  # t(g.star) %*% g.star # close to a identity matrix
  # lam.star <- as.vector(M %*% lam)
  
  g.star <- g
  lam.star <- lam
  
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
  lam0 <- lam.star
  fn <- foo(lam0, g.star, 0)
  for(i in seq(niter)){
    print(lam0)
    s <- foo(lam0, g.star, 1)
    h <- foo(lam0, g.star, 2)
    if(max(abs(s)) < 1e-6){
      print('converged')
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
    plot(fn, pch=20)
  }
  
  cat('min ev = ', min(eigen(h)$value), ', s = ', max(abs(s)), '\n')
  
  NULL
  
}
