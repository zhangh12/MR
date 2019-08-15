
test <- function(fit1, fit2, ref, bet){
  
  n <- nrow(ref)
  ref <- ref[sample(n, n, TRUE), ]
  ref <- scale(ref, scale = FALSE)
  
  L <- length(fit1$alp)
  map <- list(bet = 1, 
              a = 2, 
              alp = 2 + (1:L), 
              the = 2 + L + (1:L), 
              mu = 2 + 2 * L + (1:L), 
              gam = 2 + 3 * L + (1:L), 
              lam = 2 + 4 * L + (1:(3 * L + 1)))
  
  a <- log(n / sum(exp(bet * ref %*% fit1$alp)))
  para <- c(bet = bet, 
            #a = fit2$a, 
            a = a, 
            alp = fit1$alp, 
            the = fit1$the, 
            mu = fit2$mu, 
            gam = fit2$gam, 
            lam = runif(3 * L + 1, .0, .0))
  
  n1 <- fit2$n1
  n0 <- fit2$n0
  rho <- diag(n1/n0)
  
  lam <- inner.loop(para, map, rho, ref, fit2$pi)
  
  #lam1 <- NR.lambda(para, map, ref)
  
  NULL
  
}
