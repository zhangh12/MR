
init.cc <- function(para0, ref, the0, gam0, V){

  ref <- scale(ref, scale = FALSE)
  np <- ncol(ref)
  lam <- rep(.01, 3*np+1)
  para <- c(0, para0$a, para0$alp, para0$mu, para0$the, para0$gam, lam)
  names(para) <- c('bet', 'a',
                   paste0('alp', 1:np),
                   paste0('mu', 1:np),
                   paste0('the', 1:np),
                   paste0('gam', 1:np),
                   paste0('lam', 1:(3*np+1)))

  id.alp <- data.frame(start = 3, end = 3+np-1)
  id.mu <- data.frame(start = 3+np, end = 3+2*np-1)
  id.the <- data.frame(start = 3+2*np, end = 3+3*np-1)
  id.gam <- data.frame(start = 3+3*np, end = 3+4*np-1)
  id.lam <- data.frame(start = 3+4*np, end = 3+7*np)
  para.id <- list(id.alp = id.alp, id.mu = id.mu, id.the = id.the, id.gam = id.gam, id.lam = id.lam)

  list(para = para, para.id = para.id, ref = ref, the0 = the0, gam0 = gam0, V = V)

}
