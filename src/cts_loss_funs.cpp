#include <RcppArmadillo.h>
#include "auxfun.h" 
#include "cts_loss_funs.h" 
#include "mn_loss_funs.h"  //need for now to compile

using namespace Rcpp;
using namespace arma;

// ##################################################################
// ############## Loss Functions for OPCG and MADE ##################
// ##################################################################

// #######################################################################
// #                MADE for Continuous Response (i.e. MAVE)         
// #######################################################################



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mgauss_loss_j_made(arma::vec c, 
                         arma::mat xj, 
                         arma::mat y_matrix, 
                         arma::vec wj, 
                         arma::vec ahat,
                         arma::mat Dhat) { 
  
  
  arma::uword n=y_matrix.n_cols;
  arma::uword m=y_matrix.n_rows + 1;
  arma::mat I(m,m); I.eye();
  
  arma::mat mean_nll_j(1,1); mean_nll_j.zeros();    
  
  // Writing the For loop instead of sapply.
  arma::uword i;
  
  for (i = 0; i < n; i++ ) {
    
    arma::mat tXij_tD=kron( (xj.col(i)).t(),Dhat.t());
    arma::vec lcp=tXij_tD*c;
    
    // gaussian-loss
    // the ahat*y term is free of the parameter, so we get drop it in the estimation
    mean_nll_j += -wj(i)*( lcp.t()*y_matrix.col(i) + (ahat + lcp).t()*(ahat + lcp)/2 ) /n;
    
    // square-loss
    // mean_nll_j += -wj(i)*( (ahat + lcp).t()* y_matrix.col(i) -
    //    (ahat + lcp) ).t()*( (ahat + lcp).t()* y_matrix.col(i) -
    //    (ahat + lcp) )/n;
    
  }
  
  return mean_nll_j;
  
} ;


// A Loss function for whole MADE 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mgauss_loss_made(arma::vec c, 
                           arma::mat x_matrix, 
                           arma::mat y_matrix, 
                           double bw, 
                           Rcpp::List ahat_list,
                           Rcpp::List Dhat_list,
                           arma::mat r_mat) { 
  
  arma::uword n=y_matrix.n_cols;  
  arma::mat mean_nll(1,1); mean_nll.zeros();  
  
  // Writing the For loop instead of sapply.
  arma::uword j;
  
  for (j = 0; j < n; j++ ) {
    arma::vec ahat = ahat_list[j]; 
    arma::mat Dhat = Dhat_list[j];
    
    arma::mat xj = x_matrix;
    xj.each_col() -= xj.col(j); 
    arma::mat Bxj = r_mat.t()*xj;
    
    arma::vec wj = gauss_kern_cpp(Bxj, bw); 
    
    mean_nll += mgauss_loss_j_made(c, xj, y_matrix, wj, ahat, Dhat)/n; 
  }  
  
  return mean_nll;
  
  
} ;



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mgauss_score_j_made(arma::vec c, 
                          arma::mat xj, 
                          arma::mat y_matrix, 
                          arma::vec wj, 
                          arma::vec ahat,
                          arma::mat Dhat) {
  // # c <- c.j.ls
  arma::uword n=y_matrix.n_cols;
  arma::uword m=y_matrix.n_rows;
  arma::uword pm=c.n_elem;
  arma::mat I(m,m); I.eye();
  
  arma::mat mean_score_j(pm,1); mean_score_j.zeros();
  // arma::mat test;
  
  // Link only matters for the mean, the form of the score is
  // always the same;
  // Writing the For loop instead of sapply.
  arma::uword i;
  for (i = 0; i < n; i++ ) {
    
    arma::mat tXij_tD=kron( (xj.col(i)).t(),Dhat.t());
    arma::vec lcp=tXij_tD*c;
    arma::vec mu_ij = ahat + lcp ;
    
    mean_score_j += -wj(i)*tXij_tD.t()*( y_matrix.col(i) - mu_ij)/n;
    
    // test = -wj(i)*tVij_I.t()*( y_datta.col(i) - mu_ij)/n;
  }
  
  return mean_score_j;
  
};

// A score function for whole MADE 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mgauss_score_made(arma::vec c, 
                            arma::mat x_matrix, 
                            arma::mat y_matrix, 
                            double bw,
                            Rcpp::List ahat_list,
                            Rcpp::List Dhat_list,
                            arma::mat r_mat) {
  
  
  arma::uword n=y_matrix.n_cols; //arma::uword m=y_matrix.n_rows;
  arma::uword mp=c.n_elem;
  
  arma::mat mean_score(mp,1); mean_score.zeros();   
  
  // Writing the For loop instead of sapply.
  arma::uword j;
  
  for (j = 0; j < n; j++ ) {
    arma::vec ahat = ahat_list[j]; 
    arma::mat Dhat = Dhat_list[j];
    
    arma::mat xj = x_matrix;
    xj.each_col() -= xj.col(j); 
    arma::mat Bxj = r_mat.t()*xj;
    
    arma::vec wj = gauss_kern_cpp(Bxj, bw); 
    
    mean_score += mgauss_score_j_made(c, xj, y_matrix, wj, ahat, Dhat); 
    
  } 
  
  
  return mean_score; //pow(n,2)*
  
} ;
// 
// arma::mat mgauss_score_made(arma::vec c,
//                             arma::mat x_matrix,
//                             arma::mat y_matrix, 
//                             double bw,
//                             Rcpp::List ahat_list,
//                             Rcpp::List Dhat_list,
//                             arma::mat r_mat) { 
//   
//   arma::uword n=y_matrix.n_cols;
//   arma::uword m=y_matrix.n_rows;
//   arma::uword pm=c.n_elem;
//   arma::mat I(m,m); I.eye();
//   
//   arma::mat mean_score(pm,1); mean_score.zeros();
//   
//   // Writing the For loop instead of sapply.
//   arma::uword j;
//   
//   for (j = 0; j < n; j++ ) {
//     arma::vec ahat = ahat_list[j]; 
//     arma::mat Dhat = Dhat_list[j];
//     
//     arma::mat xj = x_matrix;
//     xj.each_col() -= xj.col(j); 
//     arma::mat Bxj = r_mat.t()*xj;
//     
//     arma::vec wj = gauss_kern_cpp(Bxj, bw); 
//     
//     mean_score += mgauss_score_j_made(c, xj, y_matrix, wj, ahat, Dhat)/n; 
//   }  
//   
//   return mean_score;
//   
//   
// } ;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mgauss_info_j_made(arma::vec c, 
                         arma::mat xj, 
                         arma::mat y_matrix, 
                         arma::vec wj, 
                         arma::vec ahat,
                         arma::mat Dhat) {
  // # c <- c.j.ls
  arma::uword n=y_matrix.n_cols;
  arma::uword m=y_matrix.n_rows;
  arma::uword pm=c.n_elem;
  arma::mat I(m,m); I.eye();
  
  arma::mat mean_info_j(pm,pm); mean_info_j.zeros();
  // arma::mat test;
  
  // Writing the For loop instead of sapply.
  arma::uword i;
  for (i = 0; i < n; i++ ) {
    
    arma::mat tXij_tD=kron( (xj.col(i)).t(),Dhat.t());
    // arma::vec lcp=tXij_tD*c;
    // arma::vec mu_ij = dot_b_multinom( ( ahat + lcp), k(i), link);
    
    mean_info_j += wj(i)*tXij_tD.t()*tXij_tD/n;
    //( diagmat(mu_ij) - mu_ij*mu_ij.t())
    // test = -wj(i)*tVij_I.t()*( y_datta.col(i) - mu_ij)/n;
  }
  
  arma::mat I_j; I_j = (mean_info_j + mean_info_j.t())/2;
  return I_j;
  
};

// A info function for whole MADE 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mgauss_info_made(arma::vec c, 
                           arma::mat x_matrix, 
                           arma::mat y_matrix, 
                           double bw, 
                           Rcpp::List ahat_list,
                           Rcpp::List Dhat_list,
                           arma::mat r_mat) { 
  
  arma::uword n=y_matrix.n_cols;
  // arma::uword m=y_matrix.n_rows;
  arma::uword pm=c.n_elem;
  // arma::mat I(m,m); I.eye();
  
  arma::mat mean_info(pm,pm); mean_info.zeros();
  
  // Writing the For loop instead of sapply.
  arma::uword j;
  
  for (j = 0; j < n; j++ ) {
    arma::vec ahat = ahat_list[j]; 
    arma::mat Dhat = Dhat_list[j];
    
    arma::mat xj = x_matrix;
    xj.each_col() -= xj.col(j); 
    arma::mat Bxj = r_mat.t()*xj;
    
    arma::vec wj = gauss_kern_cpp(Bxj, bw); 
    
    mean_info += mgauss_info_j_made(c, xj, y_matrix, wj, ahat, Dhat); 
  }  
  
  return mean_info;
  
  
} ;