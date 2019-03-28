
## score test for causal effect

lmt.gim <- function(object, level, plot = FALSE){
  
  bet <- coef(object)
  se <- sqrt(vcov(object))
  alp <- object$exposure$alp
  pi <- object$outcome$pi
  cov.alp <- object$exposure$cov.alp
  cov.pi <- object$outcome$cov.pi
  
  para <- object$para
  map <- object$map
  the0 <- object$the0
  inv.the <- object$inv.the
  pi <- object$pi
  inv.pi <- object$inv.pi
  ref <- object$ref
  
  sc <- as.vector(t(alp) %*% inv.pi %*% pi)
  x <- as.vector(rmvnorm(1e5, sigma = cov.pi) %*% inv.pi %*% alp)
  p.sc <- mean(x > abs(sc) | x < -abs(sc))
  
  fac <- qnorm((level+1)/2)
  
  stat0 <- outer.loop.confint(bet, para, map, the0, inv.the, pi, inv.pi, ref, level)
  b <- bet
  stat <- stat0
  iter <- 0
  while(TRUE){
    iter <- iter + 1
    bet1 <- bet - iter * se/2
    if(bet1 < bet - 5*fac*se){
      rt1 <- list(root = bet-qnorm((level+1)/2) *se)
      msg <- 'lower confidence limit might be inaccurate'
      warning(msg)
      break
    }
    stat1 <- outer.loop.confint(bet1, para, map, the0, inv.the, pi, inv.pi, ref, level)
    b <- c(bet1, b)
    stat <- c(stat1, stat)
    if(stat0 * stat1 < 0){
      rt1 <- uniroot(outer.loop.confint, lower = bet1, upper = bet, 
                     extendInt = 'yes', check.conv = TRUE, trace = 10, 
                     para = para, map = map, the0 = the0, inv.the = inv.the, 
                     pi = pi, inv.pi = inv.pi, ref = ref, level = level)
      break
    }
  }
  
  iter <- 0
  while(TRUE){
    iter <- iter + 1
    bet1 <- bet + iter * se/2
    if(bet1 > bet + 5*fac*se){
      rt2 <- list(root = bet+qnorm((level+1)/2) *se)
      msg <- 'upper confidence limit might be inaccurate'
      warning(msg)
      break
    }
    stat1 <- outer.loop.confint(bet1, para, map, the0, inv.the, pi, inv.pi, ref, level)
    stat <- c(stat, stat1)
    b <- c(b, bet1)
    if(stat0 * stat1 < 0){
      rt2 <- uniroot(outer.loop.confint, lower = bet, upper = bet1, 
                     extendInt = 'yes', check.conv = TRUE, trace = 10, 
                     para = para, map = map, the0 = the0, inv.the = inv.the, 
                     pi = pi, inv.pi = inv.pi, ref = ref, level = level)
      break
    }
  }
  
  if(plot){
    gap <- qchisq(level, df = 1)
    width1 <- max(b)-min(b)
    width2 <- max(stat)-min(stat)
    plot(b, stat+gap, type = 'l', pch = 20, main = 'LM-based Confidence Interval', 
         xlab = expression(paste('Causal Effect ', beta)), ylab = 'LM Statistics', 
         xlim = c(min(b)-width1*0.1, max(b)+width1*0.2), 
         ylim = c(min(stat)-width2*0.1+gap, max(stat)+width2*0.2+gap))
    abline(v = bet, col = 'blue')
    abline(h = gap, col = 'red')
    legend(x = 'topright', 
           lty = c(1, 1), 
           bty = 'n', 
           col = c('red', 'blue'), 
           legend = c('Threshold', 'Causal Effect'))
  }
  
  ci <- data.frame(LCL = rt1$root, UCL = rt2$root, P = p.sc)
  colnames(ci)[1:2] <- paste((1+c(-1,1)*level)/2*100, '%')
  rownames(ci) <- 'LM.gim'
  
  ci
  
}

