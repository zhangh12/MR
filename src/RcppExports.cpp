// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _MR_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpp_score
NumericVector rcpp_score(double bet, double a, NumericVector alp, NumericVector the, NumericVector mu, NumericVector gam, NumericVector lam, NumericMatrix x, NumericVector rho, NumericVector d_the, NumericVector d_gam);
RcppExport SEXP _MR_rcpp_score(SEXP betSEXP, SEXP aSEXP, SEXP alpSEXP, SEXP theSEXP, SEXP muSEXP, SEXP gamSEXP, SEXP lamSEXP, SEXP xSEXP, SEXP rhoSEXP, SEXP d_theSEXP, SEXP d_gamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type bet(betSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alp(alpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type the(theSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gam(gamSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d_the(d_theSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d_gam(d_gamSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_score(bet, a, alp, the, mu, gam, lam, x, rho, d_the, d_gam));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MR_rcpp_hello_world", (DL_FUNC) &_MR_rcpp_hello_world, 0},
    {"_MR_rcpp_score", (DL_FUNC) &_MR_rcpp_score, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_MR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
