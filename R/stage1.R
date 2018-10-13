
## recover joint effect of exposure model

# the - theta, marginal ols
# se - nominal SE of marginal ols, provided by lm
# ns - matrix of sharing sample sizes among marginal ols
# ref - reference panel
stage1 <- function(exposure, ns, ref){
  
  the <- exposure$beta
  se <- exposure$se
  
  L <- length(the) # number of markers
  N <- nrow(ref) # number of reference samples
  X <- scale(as.matrix(ref), scale = FALSE) # centerized reference
  rm(ref)
  
  # redefine covariance
  cov <- function(x){
    m <- nrow(x)
    stats::cov(x) * (m-1)/m
  }
  
  C <- cov(X)
  H <- diag(diag(C))
  Vx <- diag(C) # var(x)
  n <- diag(ns) # diagonal elements are sample size for marginal model
  Vy <- mean(((n-2)*se^2+the^2) * Vx) # var(y), averaged over markers
  
  alp <- as.vector(solve(C) %*% H %*% the) # join effect, true model y ~ x * alp + N(0, tau)
  tau <- as.vector(Vy - t(alp) %*% C %*% alp)*N/(N-L-1) # variance of error term in true model
  
  Xalp <- as.vector(X %*% alp)
  EU <- matrix(NA, N, L)
  for(l in 1:L){
    EU[, l] <- X[, l] * (Xalp - X[, l] * the[l])
  }
  
  # cov(U)
  D <- cov(EU)
  B <- D + C * tau # X has been centerized, otherwise use the following one
  # B <- (t(EU) %*% EU)/N + (t(X) %*% X)/N * tau
  
  cov.the <- solve(H * ns) %*% (B * ns) %*% solve(H * ns)
  se.the <- sqrt(diag(cov.the))
  
  cov.alp <- solve(C) %*% (D/N + H %*% cov.the %*% H) %*% solve(C)
  se.alp <- sqrt(diag(cov.alp))
  
  rs <- exposure$iv
  names(alp) <- rs
  names(se.alp) <- rs
  names(the) <- rs
  names(se.the) <- rs
  rownames(cov.alp) <- rs
  colnames(cov.alp) <- rs
  rownames(cov.the) <- rs
  colnames(cov.the) <- rs
  
  list(alp = alp, se.alp = se.alp, cov.alp = cov.alp, tau=tau, 
       the = the, se.the = se.the, cov.the = cov.the)
  
}
