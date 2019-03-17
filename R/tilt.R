
tilt <- function(para, map, ref){
  
  L <- (length(para) - 3) / 7
  n <- nrow(ref)
  
  bet <- para[map$bet]
  a <- para[map$a]
  alp <- para[map$alp]
  mu <- para[map$mu]
  gam <- para[map$gam]
  
  ref <- as.matrix(ref)
  
  Delta <- as.vector(exp(a + bet * ref %*% alp))
  delta <- exp(t(t(ref) * gam + mu))
  
  list(Delta = Delta, delta = delta)
  
}
