
## score test for causal effect

lmt <- function(bet, se, alp, inv.alp, pi, inv.pi, level, plot = FALSE){
  
  stat.null <- lmt.stat(0, alp, inv.alp, pi, inv.pi, level = 0.0)
    
  p <- pchisq(stat.null, df = 1, lower.tail = FALSE)
  
  fac <- qnorm((level+1)/2)
  
  stat0 <- lmt.stat(bet, alp, inv.alp, pi, inv.pi, level)
  b <- bet
  stat <- stat0
  iter <- 0
  while(TRUE){
    iter <- iter + 1
    bet1 <- bet - iter * se/2
    if(bet1 < bet - 5*fac*se){
      rt1 <- list(root = bet1)
      msg <- 'lower confidence limit might be inaccurate'
      warning(msg)
      break
    }
    stat1 <- lmt.stat(bet1, alp, inv.alp, pi, inv.pi, level)
    b <- c(bet1, b)
    stat <- c(stat1, stat)
    if(stat0 * stat1 < 0){
      rt1 <- uniroot(lmt.stat, lower = bet1, upper = bet, 
                     extendInt = 'yes', check.conv = TRUE, trace = 10, 
                     alp=alp, inv.alp=inv.alp, pi=pi, inv.pi=inv.pi, level=level)
      break
    }
  }
  
  iter <- 0
  while(TRUE){
    iter <- iter + 1
    bet1 <- bet + iter * se/2
    if(bet1 > bet + 5*fac*se){
      rt2 <- list(root = bet1)
      msg <- 'upper confidence limit might be inaccurate'
      warning(msg)
      break
    }
    stat1 <- lmt.stat(bet1, alp, inv.alp, pi, inv.pi, level)
    stat <- c(stat, stat1)
    b <- c(b, bet1)
    if(stat0 * stat1 < 0){
      rt2 <- uniroot(lmt.stat, lower = bet, upper = bet1, 
                     extendInt = 'yes', check.conv = TRUE, trace = 10, 
                     alp=alp, inv.alp=inv.alp, pi=pi, inv.pi=inv.pi, level=level)
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
  
  ci <- data.frame(LCL = rt1$root, UCL = rt2$root, P = p)
  colnames(ci)[1:2] <- paste((1+c(-1,1)*level)/2*100, '%')
  rownames(ci) <- 'LM'
  
  ci
  
}

