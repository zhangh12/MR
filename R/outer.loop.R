

outer.loop <- function(para, map, the0, inv.the, pi, inv.pi, ref){
  
  L <- (length(para) - 1) / 3
  
  n <- nrow(ref)
  ref <- ref[sample(n, n, TRUE), ]
  ref <- scale(ref, scale = FALSE)
  
  x0 <- para[c(map$bet, map$alp, map$the)]
  if(abs(x0[1]) > 1){
    x0[1] <- 0
  }
  niter <- 50
  tol <- 1e-6
  
  lam <- NR.lambda(x0, map, ref)
  fn <- obj.gim(x0, map, lam = lam, the0, inv.the, pi, inv.pi, ref)
  
  for(i in seq(niter)){
    bet <- x0[map$bet]
    alp <- x0[map$alp]
    the <- x0[map$the]
    
    ## inner loop
    lam <- NR.lambda(x0, map, ref)
    s0 <- rcpp_deriv1(bet, alp, the, lam, the0, inv.the, pi, inv.pi, ref)
    
    R11 <- rcpp_R11(bet, alp, the, lam, the0, inv.the, pi, inv.pi, ref)
    R12 <- rcpp_R12(bet, alp, the, lam, the0, inv.the, pi, inv.pi, ref)
    R22 <- rcpp_R22(bet, alp, the, lam, the0, inv.the, pi, inv.pi, ref)
    
    h <- R11 - R12 %*% solve(R22) %*% t(R12)
    #h1 <- hessian(obj.gim, c(x0,lam), map=map,lam=NULL,the0=the0,inv.the=inv.the,pi=pi,inv.pi=inv.pi,ref=ref)
    #h2 <- jacobian(score.gim, c(x0,lam), map=map,lam=NULL,the0=the0,inv.the=inv.the,pi=pi,inv.pi=inv.pi,ref=ref)
    #R11 <- h2[1:(2 * L + 1), 1:(2 * L + 1)]
    #R12 <- h2[1:(2 * L + 1), (2 * L + 2):(3 * L + 1)]
    #R22 <- h2[(2 * L + 2):(3 * L + 1), (2 * L + 2):(3 * L + 1)]
    
    ei <- eigen(h)
    if(any(ei$value >= -1e-2)){
      message('modifying hess')
      ev <- ei$values
      ev <- ev * (ev < -1e-2) - .01 * (ev >= -1e-2)
      h <- ei$vectors %*% diag(ev) %*% t(ei$vectors)
    }
    
    max.ev <- max(eigen(h)$value)
    
    deriv1 <- max(abs(s0))
    max.fn <- tail(fn, 1)
    
    if(deriv1 < tol && max.ev < 0){
      ## converge
      break
    }
    
    #s1 <- grad(obj.gim, c(x0,lam), map=map,lam=NULL,the0=the0,inv.the=inv.the,pi=pi,inv.pi=inv.pi,ref=ref)
    #s2 <- score.gim(c(x0,lam), map=map,lam=NULL,the0=the0,inv.the=inv.the,pi=pi,inv.pi=inv.pi,ref=ref)
    
    d <- solve(h, s0)
    
    for(j in 1:20){
      tau <- .5^(j-1)
      x1 <- x0 - tau * d
      g <- gfunction(x1, map, ref)
      lam1 <- NR.lambda(x1, map, ref)
      tmp <- as.vector(1 + g %*% lam1)
      fn1 <- obj.gim(x1, map, lam1, the0, inv.the, pi, inv.pi, ref)
      if(all(tmp > 1/n) && fn1 > max.fn){
        break
      }
    }
    
    x0 <- x1
    fn <- c(fn, obj.gim(x1, map, lam1, the0, inv.the, pi, inv.pi, ref))
    #print(max(abs(s)))
    cat('bet=', x0[map$bet], '\tgrad=', deriv1, '\tstep=', tau, '\tmax.ev=', max.ev, '\tfn=', max(fn), '\n')
    plot(fn, pch=20, main = 'Outer Loop')
    
  }
  
  deriv1 <- formatC(max(abs(s0)), digits = 1, format = 'e')
  max.ev <- formatC(max.ev, digits = 1, format = 'e')
  max.fn <- formatC(tail(fn, 1), digits = 1, format = 'e')
  
  if(max(abs(s0)) < tol && max.ev < 0){
    code <- 0
    message(paste0('outer loop converges. grad=', deriv1, ', max.ev=', max.ev, ', max.fn=', max.fn))
  }else{
    if(max(abs(s0)) >= tol){
      code <- 1
      message(paste0('outer loop cannot find a stationary point. grad = ', deriv1, ', max.ev = ', max.ev, '\n'))
    }else{
      if(max.ev >= 0){
        code <- 2
        message(paste0('outer loop stops at a saddle point. grad = ', deriv1, ', max.ev = ', max.ev, '\n'))
      }
    }
  }
  
  lam <- NR.lambda(x0, map, ref)
  para <- c(x0, lam)
  V <- cov.gim(para, map, the0, inv.the, pi, inv.pi, ref)
  
  list(para = para, V = V, 
       bet.gim = para[map$bet], 
       se.gim = sqrt(diag(V))[map$bet], 
       code = code)
  
}

