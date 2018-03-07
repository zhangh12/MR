
score.cc <- function(para, para.id, ref, the0, gam0, inv.V, nsample){

  id.alp <- para.id$id.alp
  id.mu <- para.id$id.mu
  id.the <- para.id$id.the
  id.gam <- para.id$id.gam
  id.lam <- para.id$id.lam

  bet <- para[1]
  a <- para[2]
  alp <- para[min(id.alp):max(id.alp)]
  mu <- para[min(id.mu):max(id.mu)]
  the <- para[min(id.the):max(id.the)]
  gam <- para[min(id.gam):max(id.gam)]
  lam <- para[min(id.lam):max(id.lam)]

  n <- nrow(ref)
  np <- ncol(ref)

  n0 <- nsample$n0
  n1 <- nsample$n1
  nz <- nsample$nz

  zlin <- as.vector(ref %*% alp)
  delta <- as.vector(exp(a + bet * ref %*% alp))
  delta1 <- matrix(NA, n, np)
  dem <- matrix(NA, n, np)
  dem1 <- matrix(NA, n, np)
  for(j in 1:np){
    delta1[, j] <- exp(mu[j] + ref[, j] * gam[j])
    dem1[, j] <- n0[j] + n1[j] * delta1[, j]
    dem[, j] <- n0[j] + n1[j] * delta
  }

  g <- gfunction.cc(ref, the, zlin, delta, delta1, dem1)
  g.bet <- gfunction.bet.cc(ref, zlin, delta, dem1)
  g.a <- gfunction.a.cc(ref, delta, dem1)
  g.alp <- gfunction.alp.cc(ref, bet, delta, dem1)
  g.mu <- gfunction.mu.cc(ref, delta1, dem, dem1)
  g.the <- gfunction.the.cc(ref)
  g.gam <- gfunction.gam.cc(ref, delta1, dem, dem1)

  pr <- as.vector(1/(1+g %*% lam))

  npara <- length(para)
  sc <- rep(NA, npara)
  names(sc) <- names(para)

  tmp <- as.vector(g.bet %*% lam)
  sc[1] <- -sum(tmp * pr)

  tmp <- as.vector(g.a %*% lam)
  sc[2] <- -sum(tmp * pr)

  k <- 2
  for(i in 1:np){
    tmp <- as.vector(g.alp[[i]] %*% lam)
    k <- k+1
    sc[k] <- -sum(tmp * pr)
  }

  k <- max(id.alp)
  for(i in 1:np){
    tmp <- as.vector(g.mu[[i]] %*% lam)
    k <- k+1
    sc[k] <- -sum(tmp * pr)
  }

  pi <- c(the, gam)
  pi0 <- c(the0, gam0)
  dqf <- as.vector(inv.V %*% (pi - pi0))
  k <- max(id.mu)
  for(i in 1:np){
    tmp <- as.vector(g.the[[i]] %*% lam)
    k <- k+1
    sc[k] <- -sum(tmp * pr) - dqf[i]
  }

  k <- max(id.the)
  for(i in 1:np){
    tmp <- as.vector(g.gam[[i]] %*% lam)
    k <- k+1
    sc[k] <- -sum(tmp * pr) - dqf[i + np]
  }

  sc[min(id.lam):max(id.lam)] <- -as.vector(t(g) %*% pr)

  sc

}
