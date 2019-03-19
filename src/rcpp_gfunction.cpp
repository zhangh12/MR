
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rcpp_gfunction(double bet, 
                             double a, 
                             NumericVector alp, 
                             NumericVector the, 
                             NumericVector mu, 
                             NumericVector gam, 
                             NumericVector lam, 
                             NumericMatrix x, 
                             NumericVector rho, 
                             NumericVector pi) {
  
  int L = x.ncol();
  int n = x.nrow();
  
  NumericVector lin(n); // X^T alp
  NumericVector Delta(n);
  NumericMatrix delta(n, L); 
  NumericMatrix one_plus_rho_delta(n, L);
  
  NumericMatrix g(n, 3 * L + 1); 
  
  for(int i = 0; i < n; ++i){
    double tmp = .0;
    for(int l = 0; l < L; ++l){
      tmp += x(i, l) * pi[l]; 
      lin[i] += x(i, l) * alp[l];
      delta(i, l) = exp(mu[l] + x(i, l) * gam[l]);
      one_plus_rho_delta(i, l) = 1 + rho[l] * delta(i, l);
    }
    Delta[i] = exp(a + bet * lin[i]); 
    //Delta[i] = exp(a + tmp);
  }
  
  for(int i = 0; i < n; ++i){
    for(int l = 0; l < L; ++l){
      g(i, l) = (lin[i] - x(i, l) * the[l]) * x(i, l); 
      g(i, l + L) = (Delta[i] - delta(i, l)) / one_plus_rho_delta(i, l);
      g(i, l + 2 * L) = g(i, l + L) * x(i, l);
    }
    g(i, 3 * L) = Delta[i] - 1;
  }
  
  return g; 
  
}
