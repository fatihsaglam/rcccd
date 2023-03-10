// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// f_cover_rcpp
List f_cover_rcpp(IntegerVector cover, double thresh, NumericMatrix M, NumericMatrix dist_main2main, NumericVector dist_main2other);
RcppExport SEXP _rcccd_f_cover_rcpp(SEXP coverSEXP, SEXP threshSEXP, SEXP MSEXP, SEXP dist_main2mainSEXP, SEXP dist_main2otherSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type cover(coverSEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type M(MSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dist_main2main(dist_main2mainSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dist_main2other(dist_main2otherSEXP);
    rcpp_result_gen = Rcpp::wrap(f_cover_rcpp(cover, thresh, M, dist_main2main, dist_main2other));
    return rcpp_result_gen;
END_RCPP
}
// radii_rw
arma::vec radii_rw(const arma::mat& dist_main2all, const arma::mat& dist_main2main, const arma::vec& w, const double m, const int p);
RcppExport SEXP _rcccd_radii_rw(SEXP dist_main2allSEXP, SEXP dist_main2mainSEXP, SEXP wSEXP, SEXP mSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type dist_main2all(dist_main2allSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type dist_main2main(dist_main2mainSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double >::type m(mSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(radii_rw(dist_main2all, dist_main2main, w, m, p));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rcccd_f_cover_rcpp", (DL_FUNC) &_rcccd_f_cover_rcpp, 5},
    {"_rcccd_radii_rw", (DL_FUNC) &_rcccd_radii_rw, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_rcccd(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
