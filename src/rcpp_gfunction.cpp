
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rcpp_gfunction(NumericVector alp, 
                             NumericVector the, 
                             NumericMatrix x) {
  
  int L = x.ncol();
  int n = x.nrow();
  
  NumericVector lin(n); // X^T alp
  
  NumericMatrix g(n, L); 
  
  for(int i = 0; i < n; ++i){
    for(int l = 0; l < L; ++l){
      lin[i] += x(i, l) * alp[l];
    }
    
    for(int l = 0; l < L; ++l){
      g(i, l) = (lin[i] - x(i, l) * the[l]) * x(i, l);
    }
  }
  
  return g; 
  
}


