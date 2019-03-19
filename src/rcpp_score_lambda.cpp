
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

/*
 * this function returns components in rcpp_score that corresponds to lambda
 * will be used to resolve dual problem subject to lambda
 */

NumericVector rcpp_score_lambda(double bet, 
                                double a, 
                                NumericVector alp, 
                                NumericVector the, 
                                NumericVector mu, 
                                NumericVector gam, 
                                NumericVector lam, 
                                NumericMatrix x, 
                                NumericVector rho) {
  
  int L = x.ncol();
  int n = x.nrow();
  
  NumericVector sc(3 * L + 1); // to be returned, score of lambda
  
  NumericVector lin(n); // X^T alp
  NumericVector Delta(n);
  NumericMatrix delta(n, L); 
  NumericMatrix one_plus_rho_delta(n, L);
  
  for(int i = 0; i < n; ++i){
    for(int l = 0; l < L; ++l){
      lin[i] += x(i, l) * alp[l];
      delta(i, l) = exp(mu[l] + x(i, l) * gam[l]);
      one_plus_rho_delta(i, l) = 1 + rho[l] * delta(i, l);
    }
    Delta[i] = exp(a + bet * lin[i]); 
  }
  
  for(int i = 0; i < n; ++i){
    NumericVector g(3 * L + 1);
    for(int l = 0; l < L; ++l){
      g[l] = (lin[i] - x(i, l) * the[l]) * x(i, l); 
      g[l + L] = (Delta[i] - delta(i, l)) / one_plus_rho_delta(i, l);
      g[l + 2 * L] = g[l + L] * x(i, l);
    }
    g[3 * L] = Delta[i] - 1; 
    
    double one_plus_g_lam = 1.0;
    for(int j = 0; j < 3 * L + 1; ++j){
      one_plus_g_lam += g[j] * lam[j];
    }
    
    for(int j = 0; j < 3 * L + 1; ++j){
      // lam
      sc[j] -= g[j] / one_plus_g_lam;
    }
  }
  
  return sc; 
  
}
