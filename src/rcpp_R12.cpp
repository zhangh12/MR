
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rcpp_R12(double bet, 
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
  
  NumericMatrix R12(2 * L + 1, L);
  
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
    
    NumericVector lam_g_alp(L);
    NumericVector lam_g_the(L);
    for(int l = 0; l < L; ++l){
      for(int k = 0; k < L; ++k){
        lam_g_alp[l] += lam[k] * x(i, k) * x(i, l);
      }
      lam_g_the[l] = - x(i, l) * x(i, l) * lam[l];
    }
    
    for(int l = 0; l < L; ++l){
      for(int k = 0; k < L; ++k){
        R12(1 + l, k) += -x(i, l) * x(i, k) / one_plus_lam_g + lam_g_alp[l] * g[k] / one_plus_lam_g2;
        R12(1 + L + l, k) += lam_g_the[l] * g[k] / one_plus_lam_g2;
        if(l == k){
          R12(1 + L + l, k) -= -x(i, l) * x(i, l) / one_plus_lam_g;
        }
      }
    }
  }
  
  return R12; 
  
}
