// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// softthreshold
double softthreshold(double x, double lambda);
RcppExport SEXP _LogRatioReg_softthreshold(SEXP xSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(softthreshold(x, lambda));
    return rcpp_result_gen;
END_RCPP
}
// gd_naive
arma::vec gd_naive(arma::mat x, arma::vec y, double l, arma::vec beta);
RcppExport SEXP _LogRatioReg_gd_naive(SEXP xSEXP, SEXP ySEXP, SEXP lSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type l(lSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(gd_naive(x, y, l, beta));
    return rcpp_result_gen;
END_RCPP
}
// gd_cov
arma::vec gd_cov(arma::mat xx, arma::vec xy, int n, double l, arma::vec beta);
RcppExport SEXP _LogRatioReg_gd_cov(SEXP xxSEXP, SEXP xySEXP, SEXP nSEXP, SEXP lSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xy(xySEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type l(lSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(gd_cov(xx, xy, n, l, beta));
    return rcpp_result_gen;
END_RCPP
}
// gd_cov_al
arma::vec gd_cov_al(arma::mat xx, arma::vec xy, int n, double l, arma::vec beta, double mu, double alpha, bool intercept);
RcppExport SEXP _LogRatioReg_gd_cov_al(SEXP xxSEXP, SEXP xySEXP, SEXP nSEXP, SEXP lSEXP, SEXP betaSEXP, SEXP muSEXP, SEXP alphaSEXP, SEXP interceptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xy(xySEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type l(lSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type intercept(interceptSEXP);
    rcpp_result_gen = Rcpp::wrap(gd_cov_al(xx, xy, n, l, beta, mu, alpha, intercept));
    return rcpp_result_gen;
END_RCPP
}
// linear_lasso
Rcpp::List linear_lasso(arma::mat x, arma::vec y, int len);
RcppExport SEXP _LogRatioReg_linear_lasso(SEXP xSEXP, SEXP ySEXP, SEXP lenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    rcpp_result_gen = Rcpp::wrap(linear_lasso(x, y, len));
    return rcpp_result_gen;
END_RCPP
}
// linear_lasso_al
Rcpp::List linear_lasso_al(arma::mat x, arma::vec y, int len, double mu, int ub, arma::vec lambda, bool intercept);
RcppExport SEXP _LogRatioReg_linear_lasso_al(SEXP xSEXP, SEXP ySEXP, SEXP lenSEXP, SEXP muSEXP, SEXP ubSEXP, SEXP lambdaSEXP, SEXP interceptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< bool >::type intercept(interceptSEXP);
    rcpp_result_gen = Rcpp::wrap(linear_lasso_al(x, y, len, mu, ub, lambda, intercept));
    return rcpp_result_gen;
END_RCPP
}
// logistic_lasso_al
Rcpp::List logistic_lasso_al(arma::mat x, arma::vec y, int len, double mu, int ub, arma::vec lambda);
RcppExport SEXP _LogRatioReg_logistic_lasso_al(SEXP xSEXP, SEXP ySEXP, SEXP lenSEXP, SEXP muSEXP, SEXP ubSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(logistic_lasso_al(x, y, len, mu, ub, lambda));
    return rcpp_result_gen;
END_RCPP
}
// cox_lasso_al
Rcpp::List cox_lasso_al(arma::mat x, arma::vec t, arma::vec d, arma::vec tj, int len, double mu, int ub, arma::vec lambda, double devnull);
RcppExport SEXP _LogRatioReg_cox_lasso_al(SEXP xSEXP, SEXP tSEXP, SEXP dSEXP, SEXP tjSEXP, SEXP lenSEXP, SEXP muSEXP, SEXP ubSEXP, SEXP lambdaSEXP, SEXP devnullSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tj(tjSEXP);
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type devnull(devnullSEXP);
    rcpp_result_gen = Rcpp::wrap(cox_lasso_al(x, t, d, tj, len, mu, ub, lambda, devnull));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LogRatioReg_softthreshold", (DL_FUNC) &_LogRatioReg_softthreshold, 2},
    {"_LogRatioReg_gd_naive", (DL_FUNC) &_LogRatioReg_gd_naive, 4},
    {"_LogRatioReg_gd_cov", (DL_FUNC) &_LogRatioReg_gd_cov, 5},
    {"_LogRatioReg_gd_cov_al", (DL_FUNC) &_LogRatioReg_gd_cov_al, 8},
    {"_LogRatioReg_linear_lasso", (DL_FUNC) &_LogRatioReg_linear_lasso, 3},
    {"_LogRatioReg_linear_lasso_al", (DL_FUNC) &_LogRatioReg_linear_lasso_al, 7},
    {"_LogRatioReg_logistic_lasso_al", (DL_FUNC) &_LogRatioReg_logistic_lasso_al, 6},
    {"_LogRatioReg_cox_lasso_al", (DL_FUNC) &_LogRatioReg_cox_lasso_al, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_LogRatioReg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
