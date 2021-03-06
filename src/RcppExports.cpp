// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// fit_mu_lasso
arma::vec fit_mu_lasso(arma::vec omega, double gamma, arma::mat W, arma::vec z, arma::vec betaInit, bool activeSet);
RcppExport SEXP _hdme_fit_mu_lasso(SEXP omegaSEXP, SEXP gammaSEXP, SEXP WSEXP, SEXP zSEXP, SEXP betaInitSEXP, SEXP activeSetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type betaInit(betaInitSEXP);
    Rcpp::traits::input_parameter< bool >::type activeSet(activeSetSEXP);
    rcpp_result_gen = Rcpp::wrap(fit_mu_lasso(omega, gamma, W, z, betaInit, activeSet));
    return rcpp_result_gen;
END_RCPP
}
