
confint.mr <- function(object, parm, level = 0.95, ...){
  
  cf <- coef(object)
  se <- sqrt(vcov(object))
  fac <- qnorm((level+1)/2)
  wald.test <- data.frame(LCL = cf - fac * se, 
                          UCL = cf + fac * se, 
                          P = summary(object)$coef$'Pr(>|z|)')
  colnames(wald.test)[1:2] <- paste((1+c(-1,1)*level)/2*100, '%')
  rownames(wald.test) <- 'Wald'
  
  lm.test <- lmt(object, level, plot = isTRUE(list(...)$plot))
  ci <- rbind(lm.test, wald.test)
  
  ci
  
}
