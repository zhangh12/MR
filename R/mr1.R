
## exposure: data frame with columns iv, beta, se, and study (optional)
##        n: sample size in each study of exposure
##  outcome: data frame with columns iv, beta, and study (optional)
##       n1: number of cases in each study of outcome
##       n0: number of controls in each study of outcome
##      ref: a data frame of reference panel of iv
mr1 <- function(exposure, n, outcome, n1, n0, ref){
  
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
  
  test(fit1, fit2, ref, fit3$ini$para[1])
  
  NULL
  
}
