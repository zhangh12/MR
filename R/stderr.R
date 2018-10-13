
## compute SE of causal effect

# alp - joint effect of exposure model
# inv.alp - precise matrix of alp
# pi - joint effect of case-control model (logOR)
# inv.pi - precise matrix of pi
stderr <- function(par, alp, inv.alp, pi, inv.pi){
  
  bet <- par[1]
  alp1 <- par[-1]
  
  L <- length(alp)
  info <- matrix(0, L+1, L+1)
  info[1, 1] <- t(alp1) %*% inv.pi %*% alp1
  info[-1, 1] <- bet * inv.pi %*% alp1
  info[1, -1] <- t(info[-1, 1])
  info[-1, -1] <- inv.alp + bet^2 * inv.pi
  
  hess <- matrix(0, L+1, L+1)
  hess[1, 1] <- info[1, 1]
  hess[-1, 1] <- inv.pi %*% (2 * bet * alp1 - pi)
  hess[1, -1] <- t(hess[-1, 1])
  hess[-1, -1] <- info[-1, -1]
  
  ## slightly conservative
  # se <- sqrt((solve(hess) %*% info %*% solve(hess))[1, 1])
  
  
  # fit <- optim(par, obj.stage3, dev.stage3, method = 'BFGS',
  #              alp = alp, inv.alp = inv.alp,
  #              pi = pi, inv.pi = inv.pi, hessian = TRUE)
  ## equivalent to sqrt(solve(fit$hessian)[1,1])
  cov <- solve(hess)
  rownames(cov) <- names(par)
  colnames(cov) <- names(par)
  se <- sqrt(cov[1,1])
  
  list(cov.alp = cov[-1,-1], se = se)
  
}

