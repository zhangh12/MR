#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double rcpp_obj(double bet, 
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
  
  double obj = .0;
  for(int i = 0; i < n; ++i){
    double lin = .0; // X^T alp
    for(int l = 0; l < L; ++l){
      lin += x(i, l) * alp[l];
    }
    
    double one_plus_g_lam = 1.0; 
    for(int l = 0; l < L; ++l){
      one_plus_g_lam += (lin - x(i, l) * the[l]) * x(i, l) * lam[l]; 
    }
    
    obj -= log(one_plus_g_lam); 
  }
  
  NumericVector d_the(L);
  NumericVector d_pi(L);
  
  for(int l = 0; l < L; ++l){
    for(int j = 0; j < L; ++j){
      d_the[l] += (the0[j] - the[j]) * inv_the(j, l); 
      d_pi[l] += (pi[j] - bet * alp[j]) * inv_pi(j, l); 
    }
    
    obj -= .5 * d_the[l] * (the0[l] - the[l]) + .5 * d_pi[l] * (pi[l] - bet * alp[l]);
  }
  
  return obj; 
  
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
