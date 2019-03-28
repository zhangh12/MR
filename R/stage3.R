

stage3 <- function(fit1, fit2){
  
  Omega <- precision(fit1, fit2)
  
  fit <- optim(0, obj.stage3, method = 'Brent', lower = -10, upper = 10, 
               alp = fit1$alp, pi = fit2$pi, Omega = Omega, 
               hessian = TRUE)
  
  bet <- fit$par
  alp1 <- alp.stage3(bet, fit1$alp, fit2$pi, Omega)
  cov3 <- stderr(bet, alp1, fit1$alp, fit2$pi, Omega)
  
  fit <- list(exposure = list(alp = fit1$alp, cov.alp = fit1$cov.alp, 
                              the = fit1$the, cov.the = fit1$cov.the), 
              outcome = list(pi = fit2$pi, cov.pi = fit2$cov.pi, 
                             gam = fit2$gam, cov.gam = fit2$cov.gam), 
              Omega = Omega, 
              Sigma = solve(Omega), 
              coefficients = bet, 
              se = cov3$se)
  
  L <- length(fit1$alp)
  map <- list(bet = 1, 
              alp = 1 + (1:L), 
              the = 1 + L + (1:L),
              lam = 1 + 2 * L + (1:L))
  
  para <- c(bet = bet, 
            alp = fit1$alp, 
            the = fit1$the, 
            lam = runif(L, .0, .0))
  
  fit$ini <- list(L = L, 
                  map = map, 
                  para = para, 
                  inv.the = solve(fit1$cov.the), 
                  inv.pi = solve(fit2$cov.pi), 
                  the0 = fit1$the, 
                  pi = fit2$pi, 
                  cov.pi = fit2$cov.pi, 
                  alp = fit1$alp, 
                  cov.alp = fit1$cov.alp)
  
  fit
  
  
}


stage3.0 <- function(fit1, fit2){
  
  Omega <- precision(fit1, fit2)
  par <- NR.stage3(fit1$alp, fit2$pi, Omega)
  cov3 <- stderr(par, fit1$alp, fit2$pi, Omega)
  
  fit <- list(exposure = list(alp = fit1$alp, cov.alp = fit1$cov.alp, 
                              alp2 = par[-1], cov.alp2 = cov3$cov.alp), 
              outcome = list(pi = fit2$pi, cov.pi = fit2$cov.pi), 
              Omega = Omega, 
              coefficients = par[1], 
              se = cov3$se)
  
  fit
  
}
