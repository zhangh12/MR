
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rcpp_gfunction1(double bet, 
                              double a, 
                              NumericVector alp, 
                              NumericVector the, 
                              NumericVector mu, 
                              NumericVector gam, 
                              NumericVector rho, 
                              NumericMatrix x, 
                              NumericVector pi) {
  
  int L = x.ncol();
  int n = x.nrow();
  
  NumericMatrix g(n, 3 * L + 1); 
  
  for(int i = 0; i < n; ++i){
    double x_alp = .0;
    double tmp = .0;
    for(int l = 0; l < L; ++l){
      x_alp += x(i, l) * alp[l];
      tmp += x(i, l) * pi[l];
    }
    
    double Delta = exp(a + bet * x_alp); 
    //double Delta = exp(a + tmp);
    
    for(int l = 0; l < L; ++l){
      g(i, l) = (x_alp - x(i, l) * the[l]) * x(i, l);
      
      double delta = exp(mu[l] + x(i, l) * gam[l]);
      g(i, L + l) = (Delta - delta) / (1 + rho[l] * delta); 
      g(i, 2 * L + l) = g(i, L + l) * x(i, l);
    }
    
    g(i, 3 * L) = Delta - 1.0;
    
  }
  
  return g; 
  
}


