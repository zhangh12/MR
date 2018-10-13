

reformat <- function(exposure, n, outcome, n1, n0, ref){
  
  n <- as.vector(n)
  n1 <- as.vector(n1)
  n0 <- as.vector(n0)
  ref <- as.data.frame(ref)
  
  if(length(n1) != length(n0)){
    msg <- 'n1 and n0 should be identical dimension'
    stop(msg)
  }
  
  colnames(exposure) <- tolower(colnames(exposure))
  colnames(outcome) <- tolower(colnames(outcome))
  e.header <- colnames(exposure)
  o.header <- colnames(outcome)
  
  header1 <- c('iv', 'beta', 'se')
  if(!all(header1 %in% e.header)){
    msg <- paste(paste(setdiff(header1, e.header), collapse = ', '), 
                  'are not in exposure')
    stop(msg)
  }
  
  header2 <- c('iv', 'beta')
  if(!all(header2 %in% o.header)){
    msg <- paste(paste(setdiff(header2, o.header), collapse = ', '), 
                 'are not in outcome')
    stop(msg)
  }
  
  if(!('study' %in% colnames(exposure))){
    
    if(length(n) == 1){
      exposure$study <- '*'
    }else{
      msg <- 'n is not correctly specified'
      stop(msg)
    }
  }
  
  if(!('study' %in% colnames(outcome))){
    
    if(length(n1) == 1 & length(n0) == 1){
      outcome$study <- '*'
    }else{
      msg <- 'n1 or n0 is not correctly specified'
      stop(msg)
    }
  }
  
  ## from now on, both outcome and exposure have a column named study
  exposure <- exposure[, c('iv', 'beta', 'se', 'study')]
  rownames(exposure) <- exposure$iv
  outcome <- outcome[, c('iv', 'beta', 'study')]
  rownames(outcome) <- outcome$iv
  if(any(nchar(exposure$study) != length(n))){
    msg <- 'study is not correctly specified in exposure'
    stop(msg)
  }
  
  if(any(nchar(outcome$study) != nrow(n1))){
    msg <- 'study is not correctly specified in outcome'
    stop(msg)
  }
  
  ## extract common IVs
  iv <- Reduce(intersect, list(exposure$iv, outcome$iv, colnames(ref)))
  if(length(iv) == 0){
    msg <- 'no sufficient information for any of IVs'
    stop(msg)
  }
  
  iv <- sort(unique(iv))
  exposure <- exposure[iv, ]
  outcome <- outcome[iv, ]
  ref <- ref[, iv, drop = FALSE]
  
  rownames(exposure) <- NULL
  rownames(outcome) <- NULL
  
  ## sharing sample size
  n <- shared.sample(exposure$study, n)
  n1 <- shared.sample(outcome$study, n1)
  n0 <- shared.sample(outcome$study, n0)
  
  return(list(exposure = exposure, 
              n = n, 
              outcome = outcome, 
              n1 = n1, 
              n0 = n0, 
              ref = ref))
  
}
