

# Function to compute SUM(Psi) over X and the Jacobian
# true model: logit P(Y=1 | X) = a + X_1 pi_1 + ... + X_L pi_l
# marginal model: logit P(Y=1 | X) = mu_l + X_l gam_l


# X - reference panel
# n1Mat - matrix of sharing sample sizes of cases among gam
# n0Mat - matrix of sharing sample sizes of controls among gam
# gamVec - gamma, marginal logOR, nominal SE of gam is not needed
# computeJ - compute Jacobian matrix J for Newtown-Raphson algorithm
# computeSigma.gamma - compute covariance matrix of gamma, marginal logOR 
#                    - This is important because we do not require SE of gam as input
# computeSigma.theta - compute covariance matrix of recovered theta = (a, pi, mu), where pi is joint logOR
PSI0 <- function(parms, X, n1Mat, n0Mat, gamVec, computeJ=TRUE, 
                 computeSigma.gamma=FALSE, computeSigma.theta=FALSE) {
  
  # parms:       Vector of parms in the form (a, pi_1, ..., pi_L, mu_1, ..., mu_L)
  
  # n1Mat      : LxL matrix of the number of overlapping cases form SNPs i and j
  # n0Mat      : LxL matrix of the number of overlapping controls form SNPs i and j
  # X:           NxL matrix, reference
  # gamVec:  Vector of marginal log odd ratios
  
  if (computeSigma.theta) {
    computeJ           <- TRUE
    computeSigma.gamma <- TRUE
  }
  
  N      <- nrow(X) # sample size of reference
  L      <- length(gamVec) # number of variables
  L2     <- 2*L
  M      <- L2 + 1 # length of psi; Epsi = 0
  ids1   <- 2:(L+1)
  ids2   <- (L+2):M
  a      <- parms[1]
  piVec  <- matrix(parms[ids1], nrow=L, ncol=1)
  muVec  <- parms[ids2]
  rhoMat <- n1Mat/n0Mat
  rhoVec <- diag(rhoMat)
  RHOMAT <- matrix(rhoVec, nrow=N, ncol=L, byrow=TRUE)
  
  DELTA  <- exp(a + X %*% piVec) # Nx1
  tmp    <- matrix(muVec, nrow=N, ncol=L, byrow=TRUE) +
    X*matrix(gamVec, nrow=N, ncol=L, byrow=TRUE)
  delta  <- exp(tmp)
  
  # Define PSI in 3 parts, columns 1:L, (L+1):2, and 2L+1
  onePlusRhoDelta          <- 1 + delta*RHOMAT
  deltaOverOnePlusRhoDelta <- delta/onePlusRhoDelta
  DELTAMAT <- matrix(DELTA, nrow=N, ncol=L, byrow=FALSE)
  PSI1     <- (DELTAMAT - delta)/onePlusRhoDelta
  PSI2     <- PSI1*X
  PSI3     <- DELTA - 1
  PSI      <- c(myColSums(PSI1, N, L), myColSums(PSI2, N, L), sum(PSI3))
  if (!computeJ) return(PSI)
  
  #########################
  # Define the Jacobian J
  #########################
  J       <- matrix(0, nrow=M, ncol=M)
  tmp     <- DELTAMAT/onePlusRhoDelta
  tmpX    <- tmp*X
  J[, 1]  <- c(myColSums(tmp, N, L), myColSums(tmpX, N, L), sum(DELTA)) # psi.a
  
  # For pi parms
  tmpX <- DELTAMAT*X
  for (l in 1:L) {
    tmp      <- tmpX[, l]/onePlusRhoDelta
    J[, l+1] <- c(myColSums(tmp, N, L), myColSums(tmp*X, N, L), sum(tmpX[, l]))
  }
  
  # For mu parms
  K    <- L + 1
  for (l in 1:L) {
    rows         <- c(l, L + l)
    col          <- K + l
    vec          <- delta[, l]
    tmp          <- 1 + rhoVec[l]*vec
    tmpX         <- 1 + rhoVec[l]*DELTA
    tmp2         <- -vec*tmpX/(tmp*tmp)
    J[rows, col] <- c(sum(tmp2), sum(tmp2*X[, l])) 
  }
  ##########################
  
  if (!computeSigma.gamma) return(list(PSI=PSI, J=J))
  
  # (mu_1, ..., mu_L, gam_1, ..., gam_L), LxL
  J.gam  <- matrix(0, nrow=L2, ncol=L2) # J.mu,gam used for Sigma.gam, the covariance of marginal estimate gam
  I.gam  <- J.gam
  J.gam2 <- matrix(0, nrow=M, ncol=L)   # Used for cov(theta)
  
  # J.gam
  tmp  <- delta*(1 + RHOMAT*DELTAMAT)/(onePlusRhoDelta*onePlusRhoDelta)
  tmpX <- tmp*X
  tmp2 <- tmpX*X
  vec                       <- 1:L
  vec2                      <- (L+1):L2
  cvec                      <- -diag(n1Mat)/N
  diag(J.gam[vec, vec])   <- myColSums(tmp, N, L)*cvec
  tmpvec                    <- myColSums(tmpX, N, L)
  diag(J.gam[vec, vec2])  <- tmpvec*cvec
  diag(J.gam2[vec, vec])  <- -tmpvec
  diag(J.gam[vec2, vec])  <- diag(J.gam[vec, vec2])
  tmpvec                    <- myColSums(tmp2, N, L)
  diag(J.gam[vec2, vec2]) <- tmpvec*cvec
  diag(J.gam2[vec2, vec]) <- -tmpvec
  J.gam2                  <- J.gam2/N
  
  # I.gam
  for (i in 1:L) {
    tmp   <- matrix(delta[, i], nrow=N, ncol=L, byrow=FALSE)
    tmp2  <- matrix(onePlusRhoDelta[, i], nrow=N, ncol=L, byrow=FALSE)
    num   <- matrix(n1Mat[i, ], nrow=N, ncol=L, byrow=TRUE)*DELTAMAT + 
      matrix(n0Mat[i, ], nrow=N, ncol=L, byrow=TRUE)*rhoVec[i]*matrix(rhoVec, nrow=N, ncol=L, byrow=TRUE)*tmp*delta
    denom <- tmp2*onePlusRhoDelta
    tmp   <- num/denom
    tmpX  <- tmp*X
    tmp2  <- tmpX*matrix(X[, i], nrow=N, ncol=L, byrow=FALSE)
    
    I.gam[i, vec]    <- myColSums(tmp, N, L)
    I.gam[i, vec2]   <- myColSums(tmpX, N, L)
    I.gam[vec2, i]   <- I.gam[i, vec2] 
    I.gam[L+i, vec2] <- myColSums(tmp2, N, L)
  }
  I.gam <- I.gam/N
  
  tmp         <- solve(J.gam)
  tmp2        <- tmp %*% I.gam %*% tmp
  Sigma.gam <- tmp2[vec2, vec2]
  
  if (!computeSigma.theta) {
    return(list(PSI=PSI, J=J, J.gam=J.gam, I.gam=I.gam, Sigma.gam=Sigma.gam))
  }
  
  # COV(theta)
  J.theta     <- J/N
  J.theta.inv <- solve(J.theta)
  I.psi       <- cov(cbind(PSI1, PSI2, PSI3))
  tmp         <- J.gam2 %*% Sigma.gam %*% t(J.gam2)
  Sigma.theta <- J.theta.inv %*% (I.psi/N + tmp) %*% t(J.theta.inv)
  
  list(PSI=PSI, J=J, J.gam=J.gam, I.gam=I.gam, Sigma.gam=Sigma.gam, 
       Sigma.the=Sigma.theta)
  
} # END: PSI0



