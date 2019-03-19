
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rcpp_hess_lambda(double bet, 
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
  
  NumericVector sc(7 * L + 3);
  
  NumericVector lin(n); // X^T alp
  NumericVector Delta(n);
  NumericMatrix delta(n, L); 
  NumericMatrix one_plus_rho_delta(n, L);
  
  NumericMatrix g(n, 3 * L + 1); 
  NumericMatrix h(3 * L + 1, 3 * L + 1);
  
  for(int i = 0; i < n; ++i){
    for(int l = 0; l < L; ++l){
      lin[i] += x(i, l) * alp[l];
      delta(i, l) = exp(mu[l] + x(i, l) * gam[l]);
      one_plus_rho_delta(i, l) = 1 + rho[l] * delta(i, l);
    }
    Delta[i] = exp(a + bet * lin[i]); 
  }
  
  for(int i = 0; i < n; ++i){
    for(int l = 0; l < L; ++l){
      g(i, l) = (lin[i] - x(i, l) * the[l]) * x(i, l); 
      g(i, l + L) = (Delta[i] - delta(i, l)) / one_plus_rho_delta(i, l);
      g(i, l + 2 * L) = g(i, l + L) * x(i, l);
    }
    g(i, 3 * L) = Delta[i] - 1;
    
    double one_plus_g_lam = 1.0;
    for(int j = 0; j < 3 * L + 1; ++j){
      one_plus_g_lam += g(i, j) * lam[j];
    }
    
    for(int j = 0; j < 3 * L + 1; ++j){
      h(j, j) += g(i, j) * g(i, j) / one_plus_g_lam; 
      for(int k = j + 1; k < 3 * L + 1; ++k){
        double tmp = g(i, j) * g(i, k) / one_plus_g_lam;
        h(j, k) += tmp; 
        h(k, j) += tmp;
      }
    }
  }
  
  return h; 
  
}
