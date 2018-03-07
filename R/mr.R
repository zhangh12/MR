
mr <- function(){

  dat <- gen.data()
  para0 <- dat$para0
  ref <- dat$ref
  the0 <- dat$the0
  gam0 <- dat$gam0
  V <- dat$V
  nsample <- dat$nsample

  ini <- init.cc(para0, ref, the0, gam0, V)

  para <- ini$para
  para.id <- ini$para.id
  ref <- ini$ref
  the0 <- ini$the0
  gam0 <- ini$gam0
  V <- ini$V
  inv.V <- solve(V)

  fit <- NR(para, para.id, ref, the0, gam0, V, nsample)

  fit

}
