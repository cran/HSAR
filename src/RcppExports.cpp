// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// hsar_cpp_arma
List hsar_cpp_arma(arma::mat X, arma::vec y, arma::sp_mat W, arma::sp_mat M, arma::sp_mat Z, arma::mat detval, arma::mat detvalM, arma::vec Unum, int burnin, int Nsim);
RcppExport SEXP HSAR_hsar_cpp_arma(SEXP XSEXP, SEXP ySEXP, SEXP WSEXP, SEXP MSEXP, SEXP ZSEXP, SEXP detvalSEXP, SEXP detvalMSEXP, SEXP UnumSEXP, SEXP burninSEXP, SEXP NsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type detval(detvalSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type detvalM(detvalMSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Unum(UnumSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type Nsim(NsimSEXP);
    __result = Rcpp::wrap(hsar_cpp_arma(X, y, W, M, Z, detval, detvalM, Unum, burnin, Nsim));
    return __result;
END_RCPP
}
// hsar_cpp_arma_lambda_0
List hsar_cpp_arma_lambda_0(arma::mat X, arma::vec y, arma::sp_mat W, arma::sp_mat Z, arma::mat detval, arma::vec Unum);
RcppExport SEXP HSAR_hsar_cpp_arma_lambda_0(SEXP XSEXP, SEXP ySEXP, SEXP WSEXP, SEXP ZSEXP, SEXP detvalSEXP, SEXP UnumSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type detval(detvalSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Unum(UnumSEXP);
    __result = Rcpp::wrap(hsar_cpp_arma_lambda_0(X, y, W, Z, detval, Unum));
    return __result;
END_RCPP
}
// hsar_cpp_arma_rho_0
List hsar_cpp_arma_rho_0(arma::mat X, arma::vec y, arma::sp_mat M, arma::sp_mat Z, arma::mat detvalM, arma::vec Unum);
RcppExport SEXP HSAR_hsar_cpp_arma_rho_0(SEXP XSEXP, SEXP ySEXP, SEXP MSEXP, SEXP ZSEXP, SEXP detvalMSEXP, SEXP UnumSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type detvalM(detvalMSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Unum(UnumSEXP);
    __result = Rcpp::wrap(hsar_cpp_arma_rho_0(X, y, M, Z, detvalM, Unum));
    return __result;
END_RCPP
}
// sar_cpp_arma
List sar_cpp_arma(arma::mat X, arma::vec y, arma::sp_mat W, arma::mat detval, int burnin, int Nsim);
RcppExport SEXP HSAR_sar_cpp_arma(SEXP XSEXP, SEXP ySEXP, SEXP WSEXP, SEXP detvalSEXP, SEXP burninSEXP, SEXP NsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type detval(detvalSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type Nsim(NsimSEXP);
    __result = Rcpp::wrap(sar_cpp_arma(X, y, W, detval, burnin, Nsim));
    return __result;
END_RCPP
}