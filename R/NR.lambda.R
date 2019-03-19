
NR.lambda <- function(para, map, n1, n0, ref, pi){
  
  L <- (length(para) - 3) / 7
  n <- nrow(ref)
  
  bet <- para[map$bet]
  a <- para[map$a]
  alp <- para[map$alp]
  the <- para[map$the]
  mu <- para[map$mu]
  gam <- para[map$gam]
  lam <- para[map$lam]
  
  ref <- as.matrix(ref)
  ref <- ref[sample(1:n, n, TRUE), ]
  ref <- scale(ref, scale = FALSE)
  rho <- diag(n1/n0)
  
  g <- gfunction(para, map, n1, n0, ref, pi)
  #g <- g[, -c(1:L)]
  M <- chol(t(g) %*% g)
  g.star <- g %*% solve(M)
  # t(g.star) %*% g.star # close to a identity matrix
  lam.star <- as.vector(M %*% lam)
  
  foo <- function(lam, g, deriv){
    n <- nrow(g)
    m <- ncol(g)
    pr <- 1/as.vector(1 + g %*% lam)
    
    if(deriv == 0){
      return( sum(log(pr)) )
    }
    
    if(deriv == 1){
      return( as.vector(t(g) %*% pr) )
    }
    
    if(deriv == 2){
      return( -t(g) %*% (g * pr^2) )
    }
  }
  
  niter <- 100
  lam0 <- lam.star
  fn <- foo(lam0, g.star, 0)
  for(i in seq(niter)){
    s0 <- foo(lam0, g.star, 1)
    if(max(abs(s0)) < 1e-6){
      print('converged')
      break
    }
    h <- foo(lam0, g.star, 2)
    d <- solve(h, s0)
    
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
  
  NULL
  
}

# 
# if(0){
#   obj.lambda(lam, bet, a, alp, the, mu, gam, ref, rho)
#   
#   niter <- 100
#   for(i in seq(niter)){
#     
#     s <- score.lambda(lam, bet, a, alp, the, mu, gam, ref, rho)
#     
#     if(max(abs(s)) < 1e-6){
#       break
#     }
#     
#     h <- hess.lambda(lam, bet, a, alp, the, mu, gam, ref, rho)
#     o <- obj.lambda(lam, bet, a, alp, the, mu, gam, ref, rho)
#     
#     d <- as.vector(solve(h) %*% s)
#     cat(i, '\tobj=', o, ', dev1=', max(abs(s)), ', d=', range(d), '\n')
#     d <- d/sqrt(sum(d^2))
#     r <- 1
#     for(i in seq(50)){
#       r <- .9^i
#       lam1 <- lam - r * d
#       tmp <- 1/n/as.vector(1 + g %*% lam)
#       if(all(tmp > 0) && all(tmp < 1)){
#         lam <- lam1
#         break
#       }
#     }
#     
#     if(any(tmp <= 0)){
#       stop()
#     }
#     
#   }
#   
#   s <- score.lambda(lam, bet, a, alp, the, mu, gam, ref, rho)
#   if(max(abs(s)) > 1e-6){
#     stop('NR.lambda does not converge')
#   }
#   
#   lam
# }
