
cov.gim <- function(para, map, the0, inv.the, pi, inv.pi, ref){
  
  n <- nrow(ref)
  L <- (length(para) - 1) / 3
  
  bet <- para[map$bet]
  alp <- para[map$alp]
  the <- para[map$the]
  lam <- para[map$lam]
  
  U <- ind.score(para, map, the0, inv.the, pi, inv.pi, ref)
  v1 <- (n - 1) * cov(U)
  v2 <- matrix(0, 3 * L + 1, 3 * L + 1)
  v2[map$bet, map$bet] <- t(alp) %*% inv.pi %*% alp
  v2[map$bet, map$alp] <- bet * t(alp) %*% inv.pi
  v2[map$alp, map$bet] <- t(v2[map$bet, map$alp])
  v2[map$alp, map$alp] <- bet^2 * inv.pi
  v2[map$the, map$the] <- inv.the
  
  C <- v1 + v2
  
  R11 <- rcpp_R11(bet, alp, the, lam, the0, inv.the, pi, inv.pi, ref)
  R12 <- rcpp_R12(bet, alp, the, lam, the0, inv.the, pi, inv.pi, ref)
  R22 <- rcpp_R22(bet, alp, the, lam, the0, inv.the, pi, inv.pi, ref)
  
  H <- matrix(0, 3 * L + 1, 3 * L + 1)
  H[1:(2 * L + 1), 1:(2 * L + 1)] <- R11
  H[1:(2 * L + 1), (2 * L + 2):(3 * L + 1)] <- R12
  H[(2 * L + 2):(3 * L + 1), 1:(2 * L + 1)] <- t(R12)
  H[(2 * L + 2):(3 * L + 1), (2 * L + 2):(3 * L + 1)] <- R22
  
  V <- solve(H) %*% C %*% solve(H)
  
  V
  
}
