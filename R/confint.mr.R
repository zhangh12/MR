
confint.mr <- function(object, parm, level = 0.95, ...){
  
  cf <- coef(object)
  se <- sqrt(vcov(object))
  fac <- qnorm((level+1)/2)
  wald.test <- data.frame(LCL = cf - fac * se, 
                          UCL = cf + fac * se, 
                          P = summary(object)$coef$'Pr(>|z|)')
  colnames(wald.test)[1:2] <- paste((1+c(-1,1)*level)/2*100, '%')
  rownames(wald.test) <- 'Wald'
  
  alp <- object$exposure$alp
  inv.alp <- solve(object$exposure$cov.alp)
  pi <- object$outcome$pi
  inv.pi <- solve(object$outcome$cov.pi)
  
  lm.test <- lmt(bet=cf, se, alp, inv.alp, pi, inv.pi, level)
  ci <- rbind(lm.test, wald.test)
  
  ci
  
}
