
obj.cc <- function(para, para.id, ref, the0, gam0, inv.V, nsample){

  id.alp <- para.id$id.alp
  id.mu <- para.id$id.mu
  id.the <- para.id$id.the
  id.gam <- para.id$id.gam
  id.lam <- para.id$id.lam

  n <- nrow(ref)
  np <- ncol(ref)

  n0 <- nsample$n0
  n1 <- nsample$n1
  nz <- nsample$nz

  bet <- para[1]
  a <- para[2]
  alp <- para[min(id.alp):max(id.alp)]
  mu <- para[min(id.mu):max(id.mu)]
  the <- para[min(id.the):max(id.the)]
  gam <- para[min(id.gam):max(id.gam)]
  lam <- para[min(id.lam):max(id.lam)]

  zlin <- as.vector(ref %*% alp)
  delta <- as.vector(exp(a + bet * ref %*% alp))
  delta1 <- matrix(NA, n, np)
  dem1 <- matrix(NA, n, np)
  for(j in 1:np){
    delta1[, j] <- exp(mu[j] + ref[, j] * gam[j])
    dem1[, j] <- n0[j] + n1[j] * delta1[, j]
  }

  g <- gfunction.cc(ref, the, zlin, delta, delta1, dem1)

  pi <- c(the, gam)
  pi0 <- c(the0, gam0)

  -sum(log(1+g %*% lam)) - .5 * as.vector(t(pi-pi0) %*% inv.V %*% (pi-pi0))

}
