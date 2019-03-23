
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rcpp_R11(double bet, 
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
  
  NumericMatrix R11(2 * L + 1, 2 * L + 1);
  
  // initialized R11
  for(int l = 0; l < L; ++l){
    double inv_pi_alp = .0;
    for(int k = 0; k < L; ++k){
      inv_pi_alp += inv_pi(l, k) * alp[k];
      R11(0, 1 + l) += (pi[k] - 2 * bet * alp[k]) * inv_pi(k, l);
      R11(1 + l, 1 + k) = - bet * bet * inv_pi(l, k); 
      R11(1 + L + l, 1 + L + k) = -inv_the(l, k);
    }
    R11(0, 0) -= alp[l] * inv_pi_alp; 
    R11(1 + l, 0) = R11(0, 1 + l);
  }
  
  
  for(int i = 0; i < n; ++i){
    double lin = .0;
    for(int l = 0; l < L; ++l){
      lin += x(i, l) * alp[l];
    }
    
    double one_plus_lam_g = 1.0;
    for(int l = 0; l < L; ++l){
      double g = (lin - x(i, l) * the[l]) * x(i, l);
      one_plus_lam_g += lam[l] * g; 
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
        double tmp = lam_g_alp[l] * lam_g_alp[k] / one_plus_lam_g2;
        R11(1 + l, 1 + k) += tmp;
        
        tmp = lam_g_alp[l] * lam_g_the[k] / one_plus_lam_g2;
        R11(1 + l, 1 + L + k) += tmp;
        
        tmp = lam_g_the[l] * lam_g_alp[k] / one_plus_lam_g2;
        R11(1 + L + l, 1 + k) += tmp;
        
        tmp = lam_g_the[l] * lam_g_the[k] / one_plus_lam_g2;
        R11(1 + L + l, 1 + L + k) += tmp; 
      }
    }
  }
  
  return R11; 
  
}
