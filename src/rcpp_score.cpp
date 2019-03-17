
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rcpp_score(double bet, 
                  double a, 
                  NumericVector alp, 
                  NumericVector the, 
                  NumericVector mu, 
                  NumericVector gam, 
                  NumericVector lam, 
                  NumericMatrix x, 
                  NumericVector rho, 
                  NumericVector d_the, 
                  NumericVector d_gam) {
  
  int L = x.ncol();
  int n = x.nrow();
  
  NumericVector sc(7 * L + 3);
  
  NumericVector lin(n); // X^T alp
  NumericVector Delta(n);
  NumericMatrix delta(n, L); 
  NumericMatrix one_plus_rho_delta(n, L);
  
  NumericMatrix g(n, 3 * L + 1); 
  // NumericMatrix g_bet(n, 3 * L + 1);
  NumericMatrix g_a(n, 3 * L + 1);
  
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
      
      g_a(i, l + L) = Delta[i] / one_plus_rho_delta(i, l);
      g_a(i, l + 2 * L) = g_a(i, l + L) * x(i, l);
      
      /*
      g_bet(i, l + L) = g_a(i, l + L) * lin[i];
      g_bet(i, l + 2 * L) = g_a(i, l + 2 * L) * lin[i]; 
      */
      
    }
    g(i, 3 * L) = Delta[i] - 1;
    g_a(i, 3 * L) = Delta[i]; 
    // g_bet(i, 3 * L) = Delta[i] * lin[i]; 
    
    double one_plus_g_lam = 1.0;
    double g_bet_lam = .0;
    double g_a_lam = .0;
    for(int j = 0; j < 3 * L + 1; ++j){
      one_plus_g_lam += g(i, j) * lam[j];
      
      if(j >= L) {
        double tmp = g_a(i, j) * lam[j];
        g_a_lam += tmp;
        g_bet_lam += tmp * lin[i]; 
      }
    }
    
    for(int j = 0; j < 3 * L + 1; ++j){
      // lam
      sc[4 * L + 2 + j] += g(i, j) / one_plus_g_lam;
    }
    
    // bet
    sc[0] += g_bet_lam / one_plus_g_lam; 
    
    // a
    sc[1] += g_a_lam / one_plus_g_lam; 
    
    for(int k = 0; k < L; ++k){
      double g_alp_lam = .0; 
      double tmp = bet * x(i, k); 
      for(int l = 0; l < L; ++l){
        g_alp_lam += x(i, k) * x(i, l) * lam[l];
        g_alp_lam += g_a(i, l + L) * tmp * lam[l + L]; 
        g_alp_lam += g_a(i, l + 2 * L) * tmp * lam[l + 2 * L];
      }
      g_alp_lam += Delta[i] * tmp * lam[3 * L]; 
      
      // alp_k
      sc[2 + k] += g_alp_lam / one_plus_g_lam;
      
      // the_k
      double g_the_lam = -x(i, k) * x(i, k) * lam[k]; 
      sc[2 + L + k] +=  g_the_lam / one_plus_g_lam; 
      
      // mu_k
      double tmp1 = - (1 + rho[k] * Delta[i]) * delta(i, k) / (one_plus_rho_delta(i, k) * one_plus_rho_delta(i, k)); 
      double tmp2 = tmp1 * x(i, k);
      double g_mu_lam = tmp1 * lam[k + L] + tmp2 * lam[k + 2 * L];
      sc[2 + 2 * L + k] += g_mu_lam / one_plus_g_lam; 
      
      // gam_k
      double g_gam_lam = g_mu_lam * x(i, k);
      sc[2 + 3 * L + k] += g_gam_lam / one_plus_g_lam; 
    }
  }
  
  for(int l = 0; l < L; ++l){
    sc[2 + L + l] += d_the[l];
    sc[2 + 3 * L + l] += d_gam[l];
  }
  
  return sc; 
  
}
