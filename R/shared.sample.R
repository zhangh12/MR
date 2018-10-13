

shared.sample <- function(study, n){
  
  tmp <- strsplit(study, split = '')
  tmp <- sapply(tmp, function(x){paste(ifelse(x=='?',0,1),sep='',collapse = '')})
  tmp <- strsplit(tmp, '')
  tmp <- sapply(tmp, as.integer)
  if('integer' %in% class(tmp)){
    tmp <- matrix(tmp, nrow = 1)
  }
  
  n <- t(tmp) %*% (tmp * n)
  
  n
  
}

