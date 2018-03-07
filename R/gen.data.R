
gen.data <- function(){

  set.seed(1)
  N <- 1e5
  np <- 4

  maf <- .4
  pr <- c((1-maf)^2, 2*maf*(1-maf), maf^2)
  x0 <- matrix(sample(0:2, N * np, TRUE, pr), nrow = N)

  alp <- seq(.1, .2, length.out = np)
  s2 <- .1
  z0 <- x0 %*% alp + rnorm(N, sd = sqrt(s2))
  bet <- 1.5

  pr0 <- rep(1/N, N)
  pr1 <- pr0 * exp(bet * z0)
  pr1 <- pr1/sum(pr1)

  n0 <- 1000
  n1 <- 1000
  nz <- 500
  nx <- 2000

  idy <- c(sample(1:N, n0, TRUE, pr0), sample(1:N, n1, TRUE, pr1))
  x <- x0[idy, , drop = FALSE]
  colnames(x) <- paste0('x', 1:np)
  y <- c(rep(0, n0), rep(1, n1))
  cc <- data.frame(y, x)
  gam <- rep(NA, np)
  ysd <- rep(NA, np)
  for(j in 1:np){
    form <- paste0('y~x', j)
    fit <- glm(form, data = cc, family = 'binomial')
    gam[j] <- coef(fit)[paste0('x',j)]
    ysd[j] <- summary(fit)$coefficients[paste0('x',j), 'Std. Error']
  }

  idz <- sample(1:N, nz, TRUE, pr0)
  x <- x0[idz, , drop = FALSE]
  colnames(x) <- paste0('x', 1:np)
  z <- x %*% alp + rnorm(nz, sd = sqrt(s2))
  ex <- data.frame(z, x)
  the <- rep(NA, np)
  zsd <- rep(NA, np)
  for(j in 1:np){
    form <- paste0('z~x',j)
    fit <- glm(form, data = ex, family = 'gaussian')
    the[j] <- coef(fit)[paste0('x',j)]
    zsd[j] <- summary(fit)$coefficients[paste0('x',j), 'Std. Error']
  }

  idx <- sample(1:N, nx, TRUE, pr0)
  x <- x0[idx, , drop = FALSE]
  x <- scale(x, scale = FALSE)

  para0 <-list()
  para0$the <- the
  para0$alp <- as.vector(solve(t(x) %*% x) %*% diag(diag(t(x) %*% x)) %*% the)
  para0$bet <- bet
  para0$a <- log(nx/sum(exp(para0$bet * x %*% para0$alp)))
  para0$gam <- gam
  para0$mu <- rep(NA, np)
  for(j in 1:np){
    para0$mu[j] <- log(nx/sum(exp(x[, j] * para0$gam[j])))
  }

  V <- diag(c(zsd, ysd))
  n0 <- rep(n0, np)
  n1 <- rep(n1, np)
  nz <- rep(nz, np)
  nsample <- list(n0 = n0, n1 = n1, nz = nz)

  list(ref = x, the0 = the, gam0 = gam, ysd = ysd, zsd = zsd, para0 = para0, V = V,
       nsample = nsample)

}
