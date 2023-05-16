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
RcppExport SEXP _FLORAL_softthreshold(SEXP xSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(softthreshold(x, lambda));
    return rcpp_result_gen;
END_RCPP
}
// gd_cov
arma::vec gd_cov(arma::mat xx, arma::vec xy, int n, double l, arma::vec beta);
RcppExport SEXP _FLORAL_gd_cov(SEXP xxSEXP, SEXP xySEXP, SEXP nSEXP, SEXP lSEXP, SEXP betaSEXP) {
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
arma::vec gd_cov_al(arma::mat xx, arma::vec xy, int n, double l, double a, arma::vec beta, double mu, double alpha, bool adjust, unsigned int ncov);
RcppExport SEXP _FLORAL_gd_cov_al(SEXP xxSEXP, SEXP xySEXP, SEXP nSEXP, SEXP lSEXP, SEXP aSEXP, SEXP betaSEXP, SEXP muSEXP, SEXP alphaSEXP, SEXP adjustSEXP, SEXP ncovSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xy(xySEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type l(lSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type adjust(adjustSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ncov(ncovSEXP);
    rcpp_result_gen = Rcpp::wrap(gd_cov_al(xx, xy, n, l, a, beta, mu, alpha, adjust, ncov));
    return rcpp_result_gen;
END_RCPP
}
// linear_enet_al
Rcpp::List linear_enet_al(arma::mat x, arma::vec y, int len, double mu, int ub, arma::vec lambda, double a, bool adjust, unsigned int ncov, bool display_progress);
RcppExport SEXP _FLORAL_linear_enet_al(SEXP xSEXP, SEXP ySEXP, SEXP lenSEXP, SEXP muSEXP, SEXP ubSEXP, SEXP lambdaSEXP, SEXP aSEXP, SEXP adjustSEXP, SEXP ncovSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< bool >::type adjust(adjustSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ncov(ncovSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(linear_enet_al(x, y, len, mu, ub, lambda, a, adjust, ncov, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// logistic_enet_al
Rcpp::List logistic_enet_al(arma::mat x, arma::vec y, int len, double mu, int ub, arma::vec lambda, double a, bool adjust, unsigned int ncov, bool display_progress, bool loop1, bool loop2);
RcppExport SEXP _FLORAL_logistic_enet_al(SEXP xSEXP, SEXP ySEXP, SEXP lenSEXP, SEXP muSEXP, SEXP ubSEXP, SEXP lambdaSEXP, SEXP aSEXP, SEXP adjustSEXP, SEXP ncovSEXP, SEXP display_progressSEXP, SEXP loop1SEXP, SEXP loop2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< bool >::type adjust(adjustSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ncov(ncovSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    Rcpp::traits::input_parameter< bool >::type loop1(loop1SEXP);
    Rcpp::traits::input_parameter< bool >::type loop2(loop2SEXP);
    rcpp_result_gen = Rcpp::wrap(logistic_enet_al(x, y, len, mu, ub, lambda, a, adjust, ncov, display_progress, loop1, loop2));
    return rcpp_result_gen;
END_RCPP
}
// cox_enet_al
Rcpp::List cox_enet_al(arma::mat x, arma::vec t, arma::vec d, arma::vec tj, int len, double mu, int ub, arma::vec lambda, double a, bool adjust, unsigned int ncov, double devnull, bool display_progress, bool loop1, bool loop2, bool notcv);
RcppExport SEXP _FLORAL_cox_enet_al(SEXP xSEXP, SEXP tSEXP, SEXP dSEXP, SEXP tjSEXP, SEXP lenSEXP, SEXP muSEXP, SEXP ubSEXP, SEXP lambdaSEXP, SEXP aSEXP, SEXP adjustSEXP, SEXP ncovSEXP, SEXP devnullSEXP, SEXP display_progressSEXP, SEXP loop1SEXP, SEXP loop2SEXP, SEXP notcvSEXP) {
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
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< bool >::type adjust(adjustSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ncov(ncovSEXP);
    Rcpp::traits::input_parameter< double >::type devnull(devnullSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    Rcpp::traits::input_parameter< bool >::type loop1(loop1SEXP);
    Rcpp::traits::input_parameter< bool >::type loop2(loop2SEXP);
    Rcpp::traits::input_parameter< bool >::type notcv(notcvSEXP);
    rcpp_result_gen = Rcpp::wrap(cox_enet_al(x, t, d, tj, len, mu, ub, lambda, a, adjust, ncov, devnull, display_progress, loop1, loop2, notcv));
    return rcpp_result_gen;
END_RCPP
}
// cox_timedep_enet_al
Rcpp::List cox_timedep_enet_al(arma::mat x, arma::vec t0, arma::vec t1, arma::vec d, arma::vec tj, int len, double mu, int ub, arma::vec lambda, double a, bool adjust, unsigned int ncov, double devnull, bool display_progress);
RcppExport SEXP _FLORAL_cox_timedep_enet_al(SEXP xSEXP, SEXP t0SEXP, SEXP t1SEXP, SEXP dSEXP, SEXP tjSEXP, SEXP lenSEXP, SEXP muSEXP, SEXP ubSEXP, SEXP lambdaSEXP, SEXP aSEXP, SEXP adjustSEXP, SEXP ncovSEXP, SEXP devnullSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tj(tjSEXP);
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< bool >::type adjust(adjustSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ncov(ncovSEXP);
    Rcpp::traits::input_parameter< double >::type devnull(devnullSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(cox_timedep_enet_al(x, t0, t1, d, tj, len, mu, ub, lambda, a, adjust, ncov, devnull, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// fg_enet_al
Rcpp::List fg_enet_al(arma::mat x, arma::vec t0, arma::vec t1, arma::vec d, arma::vec tj, arma::vec w, int len, double mu, int ub, arma::vec lambda, double a, bool adjust, unsigned int ncov, double devnull, bool display_progress);
RcppExport SEXP _FLORAL_fg_enet_al(SEXP xSEXP, SEXP t0SEXP, SEXP t1SEXP, SEXP dSEXP, SEXP tjSEXP, SEXP wSEXP, SEXP lenSEXP, SEXP muSEXP, SEXP ubSEXP, SEXP lambdaSEXP, SEXP aSEXP, SEXP adjustSEXP, SEXP ncovSEXP, SEXP devnullSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tj(tjSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type ub(ubSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< bool >::type adjust(adjustSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ncov(ncovSEXP);
    Rcpp::traits::input_parameter< double >::type devnull(devnullSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(fg_enet_al(x, t0, t1, d, tj, w, len, mu, ub, lambda, a, adjust, ncov, devnull, display_progress));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FLORAL_softthreshold", (DL_FUNC) &_FLORAL_softthreshold, 2},
    {"_FLORAL_gd_cov", (DL_FUNC) &_FLORAL_gd_cov, 5},
    {"_FLORAL_gd_cov_al", (DL_FUNC) &_FLORAL_gd_cov_al, 10},
    {"_FLORAL_linear_enet_al", (DL_FUNC) &_FLORAL_linear_enet_al, 10},
    {"_FLORAL_logistic_enet_al", (DL_FUNC) &_FLORAL_logistic_enet_al, 12},
    {"_FLORAL_cox_enet_al", (DL_FUNC) &_FLORAL_cox_enet_al, 16},
    {"_FLORAL_cox_timedep_enet_al", (DL_FUNC) &_FLORAL_cox_timedep_enet_al, 14},
    {"_FLORAL_fg_enet_al", (DL_FUNC) &_FLORAL_fg_enet_al, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_FLORAL(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
