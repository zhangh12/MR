
summary.mr <- function(object, ...){
  
  TAB <- data.frame(Estimate = object$coefficients, 
                    StdErr = object$se, 
                    z.value = object$coefficients / object$se, 
                    p.value = pchisq((object$coefficients / object$se)^2, df = 1, lower.tail = FALSE), 
                    stringsAsFactors = FALSE)
  
  colnames(TAB) <- c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')
  rownames(TAB) <- 'Causal'
  
  res <- list(call = object$call, 
              coefficients = TAB)
  
  class(res) <- 'summary.mr'
  res
  
}

