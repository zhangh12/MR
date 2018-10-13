
mr <- function(exposure, n, outcome, n1, n0, ref){
  
  dat <- reformat(exposure, n, outcome, n1, n0, ref)
  exposure <- dat$exposure
  n <- dat$n
  outcome <- dat$outcome
  n1 <- dat$n1
  n0 <- dat$n0
  ref <- dat$ref
  
  fit1 <- stage1(exposure, n, ref)
  fit2 <- stage2(outcome, n1, n0, ref)
  fit3 <- stage3(fit1$alp, fit1$cov.alp, fit2$pi, fit2$cov.pi)
  
  fit3$call <- match.call()
  class(fit3) <- 'mr'
  
  fit3
  
}
