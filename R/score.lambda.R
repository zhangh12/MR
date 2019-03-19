
score.lambda <- function(lam, bet, a, alp, the, mu, gam, ref, rho){
  
  rcpp_score_lambda(bet, a, alp, the, mu, gam, lam, ref, rho)
  
}

