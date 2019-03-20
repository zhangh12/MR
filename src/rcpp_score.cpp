
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rcpp_score(double bet, 
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
  
  NumericVector sc(3 * L + 1);
  
  // bet
  for(int l = 0; l < L; ++l){
    double tmp = .0;
    for(int j = 0; j < L; ++j){
      tmp += (bet * alp[j] - pi[j]) * inv_pi(j, l);
      sc[1 + L + l] -= inv_the(l, j) * (the[j] - the0[j]); // for part of the_l
    }
    sc[0] -= tmp * alp[l]; // done for bet
    sc[1 + l] -= tmp * bet; // for part of alp_l
  }
  
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
    
    for(int l = 0; l < L; ++l){
      sc[1 + 2 * L + l] -= g[l] / one_plus_lam_g; // done for lam_l
      double lam_g_alp = .0;
      for(int k = 0; k < L; ++k){
        lam_g_alp += lam[k] * x(i, k) * x(i, l);
      }
      sc[1 + l] -= lam_g_alp / one_plus_lam_g; // done for alp_l
      
      double lam_g_the = - x(i, l) * x(i, l) * lam[l];
      sc[1 + L + l] -= lam_g_the / one_plus_lam_g; // done for the_l
    }
  }
  
  return sc; 
  
}
