
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rcpp_R22(double bet, 
                       NumericVector alp, 
                       NumericVector the, 
                       NumericVector lam, 
                       NumericVector the0, 
                       NumericMatrix inv_the, 
                       NumericVector pi, 
                       NumericMatrix inv_pi, 
                       NumericMatrix x) {
  
  int L = x.ncol();
  int n = x.nrow();
  
  NumericMatrix R22(L, L);
  
  for(int i = 0; i < n; ++i){
    double lin = .0;
    for(int l = 0; l < L; ++l){
      lin += x(i, l) * alp[l];
    }
    
    double one_plus_lam_g = 1.0;
    NumericVector g(L);
    for(int l = 0; l < L; ++l){
      g[l] = (lin - x(i, l) * the[l]) * x(i, l);
      one_plus_lam_g += lam[l] * g[l]; 
    }
    double one_plus_lam_g2 = one_plus_lam_g * one_plus_lam_g;
    
    for(int l = 0; l < L; ++l){
      for(int k = 0; k < L; ++k){
        R22(l, k) += g[l] * g[k] / one_plus_lam_g2; 
      }
    }
  }
  
  return R22; 
  
}
