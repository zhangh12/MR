
## exposure: data frame with columns iv, beta, se, and study (optional)
##        n: sample size in each study of exposure
##  outcome: data frame with columns iv, beta, and study (optional)
##       n1: number of cases in each study of outcome
##       n0: number of controls in each study of outcome
##      ref: a data frame of reference panel of iv
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
  fit3 <- stage3(fit1, fit2)
  
  fit3$call <- match.call()
  class(fit3) <- 'mr'
  
  ini <- fit3$ini
  
  para <- ini$para
  map <- ini$map
  inv.the <- ini$inv.the
  inv.pi <- ini$inv.pi
  the0 <- ini$the0
  pi <- ini$pi
  
  nr <- outer.loop(para, map, the0, inv.the, pi, inv.pi, ref)
  fit3$bet.gim <- nr$bet.gim
  fit3$se.gim <- nr$se.gim
  fit3$code <- nr$code
  
  #fit3$coefficients <- nr$coefficients
  #fit3$se <- nr$se
  
  fit3
  
}
