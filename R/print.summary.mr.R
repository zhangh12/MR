
print.summary.mr <- function(x, ...){
  
  cat('Call:\n')
  print(x$call)
  cat('\n')
  
  printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  
}

