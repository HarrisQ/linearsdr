// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// center_cpp
Rcpp::NumericVector center_cpp(Rcpp::NumericVector x, Rcpp::Nullable<Rcpp::NumericVector> center);
RcppExport SEXP _linearsdr_center_cpp(SEXP xSEXP, SEXP centerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type center(centerSEXP);
    rcpp_result_gen = Rcpp::wrap(center_cpp(x, center));
    return rcpp_result_gen;
END_RCPP
}
// stand_vec_cpp
Rcpp::NumericVector stand_vec_cpp(Rcpp::NumericVector x);
RcppExport SEXP _linearsdr_stand_vec_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(stand_vec_cpp(x));
    return rcpp_result_gen;
END_RCPP
}
// normalize_cpp
arma::vec normalize_cpp(arma::colvec x);
RcppExport SEXP _linearsdr_normalize_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(normalize_cpp(x));
    return rcpp_result_gen;
END_RCPP
}
// euc_norm_cpp
double euc_norm_cpp(arma::vec x);
RcppExport SEXP _linearsdr_euc_norm_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(euc_norm_cpp(x));
    return rcpp_result_gen;
END_RCPP
}
// list_sum
arma::mat list_sum(Rcpp::List listA, Rcpp::List listB);
RcppExport SEXP _linearsdr_list_sum(SEXP listASEXP, SEXP listBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type listA(listASEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type listB(listBSEXP);
    rcpp_result_gen = Rcpp::wrap(list_sum(listA, listB));
    return rcpp_result_gen;
END_RCPP
}
// matpower_cpp
arma::mat matpower_cpp(arma::mat A, double alpha);
RcppExport SEXP _linearsdr_matpower_cpp(SEXP ASEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(matpower_cpp(A, alpha));
    return rcpp_result_gen;
END_RCPP
}
// matcenter_cpp
arma::mat matcenter_cpp(arma::mat x_matrix, Rcpp::Nullable<unsigned int> index, Rcpp::Nullable<Rcpp::NumericVector> x0);
RcppExport SEXP _linearsdr_matcenter_cpp(SEXP x_matrixSEXP, SEXP indexSEXP, SEXP x0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x_matrix(x_matrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<unsigned int> >::type index(indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type x0(x0SEXP);
    rcpp_result_gen = Rcpp::wrap(matcenter_cpp(x_matrix, index, x0));
    return rcpp_result_gen;
END_RCPP
}
// eigen_cpp
Rcpp::List eigen_cpp(arma::mat A);
RcppExport SEXP _linearsdr_eigen_cpp(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(eigen_cpp(A));
    return rcpp_result_gen;
END_RCPP
}
// gev_cpp
Rcpp::List gev_cpp(arma::mat A, arma::mat B);
RcppExport SEXP _linearsdr_gev_cpp(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(gev_cpp(A, B));
    return rcpp_result_gen;
END_RCPP
}
// inv_sympd_cpp
arma::mat inv_sympd_cpp(arma::mat A);
RcppExport SEXP _linearsdr_inv_sympd_cpp(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(inv_sympd_cpp(A));
    return rcpp_result_gen;
END_RCPP
}
// sqrtmat_cpp
arma::mat sqrtmat_cpp(arma::mat A);
RcppExport SEXP _linearsdr_sqrtmat_cpp(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(sqrtmat_cpp(A));
    return rcpp_result_gen;
END_RCPP
}
// chol_cpp
arma::mat chol_cpp(arma::mat A);
RcppExport SEXP _linearsdr_chol_cpp(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(chol_cpp(A));
    return rcpp_result_gen;
END_RCPP
}
// solve_cpp
arma::mat solve_cpp(arma::mat A, arma::mat B);
RcppExport SEXP _linearsdr_solve_cpp(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_cpp(A, B));
    return rcpp_result_gen;
END_RCPP
}
// gauss_kern_cpp
arma::vec gauss_kern_cpp(arma::mat centered_data, double bw);
RcppExport SEXP _linearsdr_gauss_kern_cpp(SEXP centered_dataSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type centered_data(centered_dataSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    rcpp_result_gen = Rcpp::wrap(gauss_kern_cpp(centered_data, bw));
    return rcpp_result_gen;
END_RCPP
}
// wls_cpp
arma::mat wls_cpp(arma::mat x_matrix, arma::mat y_matrix, arma::vec weights, double reg);
RcppExport SEXP _linearsdr_wls_cpp(SEXP x_matrixSEXP, SEXP y_matrixSEXP, SEXP weightsSEXP, SEXP regSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x_matrix(x_matrixSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_matrix(y_matrixSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< double >::type reg(regSEXP);
    rcpp_result_gen = Rcpp::wrap(wls_cpp(x_matrix, y_matrix, weights, reg));
    return rcpp_result_gen;
END_RCPP
}
// mgauss_loss_j_made
arma::mat mgauss_loss_j_made(arma::vec c, arma::mat xj, arma::mat y_matrix, arma::vec wj, arma::vec ahat, arma::mat Dhat);
RcppExport SEXP _linearsdr_mgauss_loss_j_made(SEXP cSEXP, SEXP xjSEXP, SEXP y_matrixSEXP, SEXP wjSEXP, SEXP ahatSEXP, SEXP DhatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xj(xjSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_matrix(y_matrixSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type wj(wjSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ahat(ahatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Dhat(DhatSEXP);
    rcpp_result_gen = Rcpp::wrap(mgauss_loss_j_made(c, xj, y_matrix, wj, ahat, Dhat));
    return rcpp_result_gen;
END_RCPP
}
// mgauss_loss_made
arma::mat mgauss_loss_made(arma::vec c, arma::mat x_matrix, arma::mat y_matrix, double bw, Rcpp::List ahat_list, Rcpp::List Dhat_list, arma::mat r_mat);
RcppExport SEXP _linearsdr_mgauss_loss_made(SEXP cSEXP, SEXP x_matrixSEXP, SEXP y_matrixSEXP, SEXP bwSEXP, SEXP ahat_listSEXP, SEXP Dhat_listSEXP, SEXP r_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_matrix(x_matrixSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_matrix(y_matrixSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type ahat_list(ahat_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Dhat_list(Dhat_listSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type r_mat(r_matSEXP);
    rcpp_result_gen = Rcpp::wrap(mgauss_loss_made(c, x_matrix, y_matrix, bw, ahat_list, Dhat_list, r_mat));
    return rcpp_result_gen;
END_RCPP
}
// mgauss_score_j_made
arma::mat mgauss_score_j_made(arma::vec c, arma::mat xj, arma::mat y_matrix, arma::vec wj, arma::vec ahat, arma::mat Dhat);
RcppExport SEXP _linearsdr_mgauss_score_j_made(SEXP cSEXP, SEXP xjSEXP, SEXP y_matrixSEXP, SEXP wjSEXP, SEXP ahatSEXP, SEXP DhatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xj(xjSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_matrix(y_matrixSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type wj(wjSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ahat(ahatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Dhat(DhatSEXP);
    rcpp_result_gen = Rcpp::wrap(mgauss_score_j_made(c, xj, y_matrix, wj, ahat, Dhat));
    return rcpp_result_gen;
END_RCPP
}
// mgauss_info_j_made
arma::mat mgauss_info_j_made(arma::vec c, arma::mat xj, arma::mat y_matrix, arma::vec wj, arma::vec ahat, arma::mat Dhat);
RcppExport SEXP _linearsdr_mgauss_info_j_made(SEXP cSEXP, SEXP xjSEXP, SEXP y_matrixSEXP, SEXP wjSEXP, SEXP ahatSEXP, SEXP DhatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xj(xjSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_matrix(y_matrixSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type wj(wjSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ahat(ahatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Dhat(DhatSEXP);
    rcpp_result_gen = Rcpp::wrap(mgauss_info_j_made(c, xj, y_matrix, wj, ahat, Dhat));
    return rcpp_result_gen;
END_RCPP
}
// mnY_to_mvY
arma::mat mnY_to_mvY(arma::mat mn_y, arma::vec m_classes, Rcpp::String ytype);
RcppExport SEXP _linearsdr_mnY_to_mvY(SEXP mn_ySEXP, SEXP m_classesSEXP, SEXP ytypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mn_y(mn_ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type m_classes(m_classesSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type ytype(ytypeSEXP);
    rcpp_result_gen = Rcpp::wrap(mnY_to_mvY(mn_y, m_classes, ytype));
    return rcpp_result_gen;
END_RCPP
}
// emp_logit
arma::mat emp_logit(arma::mat y_matrix, arma::vec k_vec, double tune);
RcppExport SEXP _linearsdr_emp_logit(SEXP y_matrixSEXP, SEXP k_vecSEXP, SEXP tuneSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y_matrix(y_matrixSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k_vec(k_vecSEXP);
    Rcpp::traits::input_parameter< double >::type tune(tuneSEXP);
    rcpp_result_gen = Rcpp::wrap(emp_logit(y_matrix, k_vec, tune));
    return rcpp_result_gen;
END_RCPP
}
// emp_adcat
arma::mat emp_adcat(arma::mat y_matrix, double tune);
RcppExport SEXP _linearsdr_emp_adcat(SEXP y_matrixSEXP, SEXP tuneSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y_matrix(y_matrixSEXP);
    Rcpp::traits::input_parameter< double >::type tune(tuneSEXP);
    rcpp_result_gen = Rcpp::wrap(emp_adcat(y_matrix, tune));
    return rcpp_result_gen;
END_RCPP
}
// dot_b_multinom
arma::vec dot_b_multinom(arma::vec lin_can_par, int k_i, String link);
RcppExport SEXP _linearsdr_dot_b_multinom(SEXP lin_can_parSEXP, SEXP k_iSEXP, SEXP linkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type lin_can_par(lin_can_parSEXP);
    Rcpp::traits::input_parameter< int >::type k_i(k_iSEXP);
    Rcpp::traits::input_parameter< String >::type link(linkSEXP);
    rcpp_result_gen = Rcpp::wrap(dot_b_multinom(lin_can_par, k_i, link));
    return rcpp_result_gen;
END_RCPP
}
// mn_loss_j
arma::mat mn_loss_j(arma::vec c, arma::mat vj, arma::mat y_datta, arma::vec wj, double lambda, Rcpp::String link, arma::vec k);
RcppExport SEXP _linearsdr_mn_loss_j(SEXP cSEXP, SEXP vjSEXP, SEXP y_dattaSEXP, SEXP wjSEXP, SEXP lambdaSEXP, SEXP linkSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type vj(vjSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_datta(y_dattaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type wj(wjSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type link(linkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(mn_loss_j(c, vj, y_datta, wj, lambda, link, k));
    return rcpp_result_gen;
END_RCPP
}
// mn_score_j
arma::mat mn_score_j(arma::vec c, arma::mat vj, arma::mat y_datta, arma::vec wj, double lambda, Rcpp::String link, arma::vec k);
RcppExport SEXP _linearsdr_mn_score_j(SEXP cSEXP, SEXP vjSEXP, SEXP y_dattaSEXP, SEXP wjSEXP, SEXP lambdaSEXP, SEXP linkSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type vj(vjSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_datta(y_dattaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type wj(wjSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type link(linkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(mn_score_j(c, vj, y_datta, wj, lambda, link, k));
    return rcpp_result_gen;
END_RCPP
}
// mn_info_j
arma::mat mn_info_j(arma::vec c, arma::mat vj, arma::mat y_datta, arma::vec wj, double lambda, Rcpp::String link, arma::vec k);
RcppExport SEXP _linearsdr_mn_info_j(SEXP cSEXP, SEXP vjSEXP, SEXP y_dattaSEXP, SEXP wjSEXP, SEXP lambdaSEXP, SEXP linkSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type vj(vjSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_datta(y_dattaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type wj(wjSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type link(linkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(mn_info_j(c, vj, y_datta, wj, lambda, link, k));
    return rcpp_result_gen;
END_RCPP
}
// mn_loss_j_made
arma::mat mn_loss_j_made(arma::vec c, arma::mat xj, arma::mat y_matrix, arma::vec wj, arma::vec ahat, arma::mat Dhat, Rcpp::String link, arma::vec k);
RcppExport SEXP _linearsdr_mn_loss_j_made(SEXP cSEXP, SEXP xjSEXP, SEXP y_matrixSEXP, SEXP wjSEXP, SEXP ahatSEXP, SEXP DhatSEXP, SEXP linkSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xj(xjSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_matrix(y_matrixSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type wj(wjSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ahat(ahatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Dhat(DhatSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type link(linkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(mn_loss_j_made(c, xj, y_matrix, wj, ahat, Dhat, link, k));
    return rcpp_result_gen;
END_RCPP
}
// mn_loss_made
arma::mat mn_loss_made(arma::vec c, arma::mat x_matrix, arma::mat y_matrix, double bw, Rcpp::List ahat_list, Rcpp::List Dhat_list, Rcpp::String link, arma::vec k, arma::mat r_mat);
RcppExport SEXP _linearsdr_mn_loss_made(SEXP cSEXP, SEXP x_matrixSEXP, SEXP y_matrixSEXP, SEXP bwSEXP, SEXP ahat_listSEXP, SEXP Dhat_listSEXP, SEXP linkSEXP, SEXP kSEXP, SEXP r_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_matrix(x_matrixSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_matrix(y_matrixSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type ahat_list(ahat_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Dhat_list(Dhat_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type link(linkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type r_mat(r_matSEXP);
    rcpp_result_gen = Rcpp::wrap(mn_loss_made(c, x_matrix, y_matrix, bw, ahat_list, Dhat_list, link, k, r_mat));
    return rcpp_result_gen;
END_RCPP
}
// mn_score_j_made
arma::mat mn_score_j_made(arma::vec c, arma::mat xj, arma::mat y_matrix, arma::vec wj, arma::vec ahat, arma::mat Dhat, Rcpp::String link, arma::vec k);
RcppExport SEXP _linearsdr_mn_score_j_made(SEXP cSEXP, SEXP xjSEXP, SEXP y_matrixSEXP, SEXP wjSEXP, SEXP ahatSEXP, SEXP DhatSEXP, SEXP linkSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xj(xjSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_matrix(y_matrixSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type wj(wjSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ahat(ahatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Dhat(DhatSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type link(linkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(mn_score_j_made(c, xj, y_matrix, wj, ahat, Dhat, link, k));
    return rcpp_result_gen;
END_RCPP
}
// mn_score_made
arma::mat mn_score_made(arma::vec c, arma::mat x_matrix, arma::mat y_matrix, double bw, Rcpp::List ahat_list, Rcpp::List Dhat_list, Rcpp::String link, arma::vec k, arma::mat r_mat);
RcppExport SEXP _linearsdr_mn_score_made(SEXP cSEXP, SEXP x_matrixSEXP, SEXP y_matrixSEXP, SEXP bwSEXP, SEXP ahat_listSEXP, SEXP Dhat_listSEXP, SEXP linkSEXP, SEXP kSEXP, SEXP r_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_matrix(x_matrixSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_matrix(y_matrixSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type ahat_list(ahat_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Dhat_list(Dhat_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type link(linkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type r_mat(r_matSEXP);
    rcpp_result_gen = Rcpp::wrap(mn_score_made(c, x_matrix, y_matrix, bw, ahat_list, Dhat_list, link, k, r_mat));
    return rcpp_result_gen;
END_RCPP
}
// mn_info_j_made
arma::mat mn_info_j_made(arma::vec c, arma::mat xj, arma::mat y_matrix, arma::vec wj, arma::vec ahat, arma::mat Dhat, Rcpp::String link, arma::vec k);
RcppExport SEXP _linearsdr_mn_info_j_made(SEXP cSEXP, SEXP xjSEXP, SEXP y_matrixSEXP, SEXP wjSEXP, SEXP ahatSEXP, SEXP DhatSEXP, SEXP linkSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xj(xjSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_matrix(y_matrixSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type wj(wjSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ahat(ahatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Dhat(DhatSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type link(linkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(mn_info_j_made(c, xj, y_matrix, wj, ahat, Dhat, link, k));
    return rcpp_result_gen;
END_RCPP
}
// list_mean
arma::mat list_mean(Rcpp::List est_list);
RcppExport SEXP _linearsdr_list_mean(SEXP est_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type est_list(est_listSEXP);
    rcpp_result_gen = Rcpp::wrap(list_mean(est_list));
    return rcpp_result_gen;
END_RCPP
}
// aD_j_newton
arma::vec aD_j_newton(arma::vec init, arma::mat vj, arma::mat y_datta, arma::vec wj, double lambda, Rcpp::String link, arma::vec k, double tol, arma::uword max_iter, bool test);
RcppExport SEXP _linearsdr_aD_j_newton(SEXP initSEXP, SEXP vjSEXP, SEXP y_dattaSEXP, SEXP wjSEXP, SEXP lambdaSEXP, SEXP linkSEXP, SEXP kSEXP, SEXP tolSEXP, SEXP max_iterSEXP, SEXP testSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type init(initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type vj(vjSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_datta(y_dattaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type wj(wjSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type link(linkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type test(testSEXP);
    rcpp_result_gen = Rcpp::wrap(aD_j_newton(init, vj, y_datta, wj, lambda, link, k, tol, max_iter, test));
    return rcpp_result_gen;
END_RCPP
}
// vecB_hat
arma::mat vecB_hat(arma::vec c0, Rcpp::List score_list, Rcpp::List info_list);
RcppExport SEXP _linearsdr_vecB_hat(SEXP c0SEXP, SEXP score_listSEXP, SEXP info_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type score_list(score_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type info_list(info_listSEXP);
    rcpp_result_gen = Rcpp::wrap(vecB_hat(c0, score_list, info_list));
    return rcpp_result_gen;
END_RCPP
}
// aD_j_cg
arma::vec aD_j_cg(arma::vec init, arma::mat vj, arma::mat y_datta, arma::vec wj, double lambda, Rcpp::String link, arma::vec k, Rcpp::List control_list, bool test);
RcppExport SEXP _linearsdr_aD_j_cg(SEXP initSEXP, SEXP vjSEXP, SEXP y_dattaSEXP, SEXP wjSEXP, SEXP lambdaSEXP, SEXP linkSEXP, SEXP kSEXP, SEXP control_listSEXP, SEXP testSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type init(initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type vj(vjSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_datta(y_dattaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type wj(wjSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type link(linkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control_list(control_listSEXP);
    Rcpp::traits::input_parameter< bool >::type test(testSEXP);
    rcpp_result_gen = Rcpp::wrap(aD_j_cg(init, vj, y_datta, wj, lambda, link, k, control_list, test));
    return rcpp_result_gen;
END_RCPP
}
// vecB_cg
arma::vec vecB_cg(arma::vec init, arma::mat x_datta, arma::mat y_datta, double bw, Rcpp::List ahat_list, Rcpp::List Dhat_list, Rcpp::String link, arma::vec k, arma::mat r_mat, Rcpp::List control_list, bool test);
RcppExport SEXP _linearsdr_vecB_cg(SEXP initSEXP, SEXP x_dattaSEXP, SEXP y_dattaSEXP, SEXP bwSEXP, SEXP ahat_listSEXP, SEXP Dhat_listSEXP, SEXP linkSEXP, SEXP kSEXP, SEXP r_matSEXP, SEXP control_listSEXP, SEXP testSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type init(initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_datta(x_dattaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_datta(y_dattaSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type ahat_list(ahat_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Dhat_list(Dhat_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type link(linkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type r_mat(r_matSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control_list(control_listSEXP);
    Rcpp::traits::input_parameter< bool >::type test(testSEXP);
    rcpp_result_gen = Rcpp::wrap(vecB_cg(init, x_datta, y_datta, bw, ahat_list, Dhat_list, link, k, r_mat, control_list, test));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _linearsdr_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _linearsdr_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _linearsdr_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _linearsdr_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_linearsdr_center_cpp", (DL_FUNC) &_linearsdr_center_cpp, 2},
    {"_linearsdr_stand_vec_cpp", (DL_FUNC) &_linearsdr_stand_vec_cpp, 1},
    {"_linearsdr_normalize_cpp", (DL_FUNC) &_linearsdr_normalize_cpp, 1},
    {"_linearsdr_euc_norm_cpp", (DL_FUNC) &_linearsdr_euc_norm_cpp, 1},
    {"_linearsdr_list_sum", (DL_FUNC) &_linearsdr_list_sum, 2},
    {"_linearsdr_matpower_cpp", (DL_FUNC) &_linearsdr_matpower_cpp, 2},
    {"_linearsdr_matcenter_cpp", (DL_FUNC) &_linearsdr_matcenter_cpp, 3},
    {"_linearsdr_eigen_cpp", (DL_FUNC) &_linearsdr_eigen_cpp, 1},
    {"_linearsdr_gev_cpp", (DL_FUNC) &_linearsdr_gev_cpp, 2},
    {"_linearsdr_inv_sympd_cpp", (DL_FUNC) &_linearsdr_inv_sympd_cpp, 1},
    {"_linearsdr_sqrtmat_cpp", (DL_FUNC) &_linearsdr_sqrtmat_cpp, 1},
    {"_linearsdr_chol_cpp", (DL_FUNC) &_linearsdr_chol_cpp, 1},
    {"_linearsdr_solve_cpp", (DL_FUNC) &_linearsdr_solve_cpp, 2},
    {"_linearsdr_gauss_kern_cpp", (DL_FUNC) &_linearsdr_gauss_kern_cpp, 2},
    {"_linearsdr_wls_cpp", (DL_FUNC) &_linearsdr_wls_cpp, 4},
    {"_linearsdr_mgauss_loss_j_made", (DL_FUNC) &_linearsdr_mgauss_loss_j_made, 6},
    {"_linearsdr_mgauss_loss_made", (DL_FUNC) &_linearsdr_mgauss_loss_made, 7},
    {"_linearsdr_mgauss_score_j_made", (DL_FUNC) &_linearsdr_mgauss_score_j_made, 6},
    {"_linearsdr_mgauss_info_j_made", (DL_FUNC) &_linearsdr_mgauss_info_j_made, 6},
    {"_linearsdr_mnY_to_mvY", (DL_FUNC) &_linearsdr_mnY_to_mvY, 3},
    {"_linearsdr_emp_logit", (DL_FUNC) &_linearsdr_emp_logit, 3},
    {"_linearsdr_emp_adcat", (DL_FUNC) &_linearsdr_emp_adcat, 2},
    {"_linearsdr_dot_b_multinom", (DL_FUNC) &_linearsdr_dot_b_multinom, 3},
    {"_linearsdr_mn_loss_j", (DL_FUNC) &_linearsdr_mn_loss_j, 7},
    {"_linearsdr_mn_score_j", (DL_FUNC) &_linearsdr_mn_score_j, 7},
    {"_linearsdr_mn_info_j", (DL_FUNC) &_linearsdr_mn_info_j, 7},
    {"_linearsdr_mn_loss_j_made", (DL_FUNC) &_linearsdr_mn_loss_j_made, 8},
    {"_linearsdr_mn_loss_made", (DL_FUNC) &_linearsdr_mn_loss_made, 9},
    {"_linearsdr_mn_score_j_made", (DL_FUNC) &_linearsdr_mn_score_j_made, 8},
    {"_linearsdr_mn_score_made", (DL_FUNC) &_linearsdr_mn_score_made, 9},
    {"_linearsdr_mn_info_j_made", (DL_FUNC) &_linearsdr_mn_info_j_made, 8},
    {"_linearsdr_list_mean", (DL_FUNC) &_linearsdr_list_mean, 1},
    {"_linearsdr_aD_j_newton", (DL_FUNC) &_linearsdr_aD_j_newton, 10},
    {"_linearsdr_vecB_hat", (DL_FUNC) &_linearsdr_vecB_hat, 3},
    {"_linearsdr_aD_j_cg", (DL_FUNC) &_linearsdr_aD_j_cg, 9},
    {"_linearsdr_vecB_cg", (DL_FUNC) &_linearsdr_vecB_cg, 11},
    {"_linearsdr_rcpparma_hello_world", (DL_FUNC) &_linearsdr_rcpparma_hello_world, 0},
    {"_linearsdr_rcpparma_outerproduct", (DL_FUNC) &_linearsdr_rcpparma_outerproduct, 1},
    {"_linearsdr_rcpparma_innerproduct", (DL_FUNC) &_linearsdr_rcpparma_innerproduct, 1},
    {"_linearsdr_rcpparma_bothproducts", (DL_FUNC) &_linearsdr_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_linearsdr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}