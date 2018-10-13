
## Wheeler Bill developed this code, I did not change it

# X - reference panel
# n1Mat - matrix of sharing sample sizes of cases among gam
# n0Mat - matrix of sharing sample sizes of controls among gam
# gamVec - gamma, marginal logOR, nominal SE of gam is not needed
NewtonRaphson <- function(parms, X, n1Mat, n0Mat, gamVec) {
  
  tmp  <- PSI0(parms, X, n1Mat, n0Mat, gamVec)
  f    <- tmp$PSI
  Jinv <- try(solve(tmp$J), silent=FALSE)
  if (!("try-error" %in% class(Jinv))) {
    ret  <- parms - Jinv %*% matrix(f, ncol=1)
    ret  <- as.vector(ret)
  } else {
    # For now, add random small numbers to parms
    print(Jinv)
    ret <- as.vector(parms) + 0.01*runif(length(parms))
  } 
  
  ret
  
} # END: NewtonRaphson
