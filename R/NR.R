
# Newton-Raphson algorithm
NR <- function(para, para.id, ref, the0, gam0, V, nsample){

  inv.V <- solve(V)
  np <- length(para)
  para.null <- rep(NA, np)

  i <- 0
  while(i<100){

    s0 <- score.cc(para, para.id, ref, the0, gam0, inv.V, nsample)
    s1 <- grad(obj.cc, para, para.id = para.id, ref = ref, the0 = the0, gam0 = gam0, inv.V = inv.V, nsample = nsample)

    if(all(abs(s0) < 1e-6)){
      break
    }

    t0 <- try(inv.h0 <- solve(hess.cc(para, para.id, ref, the0, gam0, inv.V, nsample)), silent = TRUE)


    if('try-error' %in% class(t0)){
      return(list(coefficients = para.null, score = para.null, conv =0))
    }
    d0 <- as.vector(inv.h0 %*% s0)
    if(max(abs(d0)) > 1){
      #d0 <- d0/max(abs(d0))
    }
    para <- para - d0

    i <- i + 1
    print(s0)
  }

  #svd(inv.h0)$d

  #print(s0)

  list(coefficients = para, score = s0, conv = ifelse(all(abs(s0) < 1e-6), 1, 0))

}
