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
                double a, 
                NumericVector alp, 
                NumericVector the, 
                NumericVector mu, 
                NumericVector gam, 
                NumericVector lam, 
                NumericMatrix x, 
                NumericVector rho, 
                NumericVector the0, 
                NumericMatrix inv_the, 
                NumericVector gam0, 
                NumericMatrix inv_gam) {
  
  
  int L = x.ncol();
  int n = x.nrow();
  
  NumericVector lin(n); // X^T alp
  NumericVector Delta(n);
  NumericMatrix delta(n, L); 
  NumericMatrix one_plus_rho_delta(n, L);
  
  NumericMatrix g(n, 3 * L + 1); 
  
  for(int i = 0; i < n; ++i){
    for(int l = 0; l < L; ++l){
      lin[i] += x(i, l) * alp[l];
      delta(i, l) = exp(mu[l] + x(i, l) * gam[l]);
      one_plus_rho_delta(i, l) = 1 + rho[l] * delta(i, l);
    }
    Delta[i] = exp(a + bet * lin[i]); 
  }
  
  double obj = .0; 
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
    
    obj -= log(one_plus_g_lam); 
  }
  
  NumericVector d_the(L);
  NumericVector d_gam(L);
  
  for(int l = 0; l < L; ++l){
    for(int j = 0; j < L; ++j){
      d_the[l] += (the0[j] - the[j]) * inv_the(j, l); 
      d_gam[l] += (gam0[j] - gam[j]) * inv_gam(j, l); 
    }
    
    obj -= .5 * d_the[l] * (the0[l] - the[l]) + .5 * d_gam[l] * (gam0[l] - gam[l]);
  }
  
  return obj; 
  
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
