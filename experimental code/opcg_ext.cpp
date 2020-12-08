#include <RcppArmadillo.h>
#include "auxfun.h"
#include "lossfuns.h"

using namespace Rcpp;
using namespace arma;

/***
 * An OPCG-oriented alternative to MAVE/MADE
 */

// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
Rcpp::List rade(Rcpp::List grad_list,  
                arma::mat init_est,
                arma::mat sigma_mat,
                double l2_pen,
                Rcpp::List control_list) {
  
  /***
   * 
   */
  
  
  // Set out controls for the algorithm
  
  double tol = control_list["tol_val"];
  arma::uword max_iter = control_list["max_iter"];
  bool print_iter = control_list["print_iter"];
  
  // beta and sigma in backtracking rule (Bertsekas)
  // inital step size s_k that will be adjusted by backtracking   
  
  
  // Declaring 
  arma::uword n = grad_list.length();
  arma::uword p = init_est.n_rows;
  arma::uword d = init_est.n_cols;
  arma::mat test=grad_list[0]; arma::uword m = test.n_cols;
  arma::uword j;
  
  // nu update, given initial Gamma  
  // arma::mat test2 = grad_list[0];
  // arma::mat test3 = solve( init_est.t()*init_est , init_est.t()*test2 );
  Rcpp::List nu_est0(n);
  for (j = 0; j < n; j++) {
    arma::mat temp_j = grad_list[j];
    nu_est0[j] = solve( init_est.t()*sigma_mat*sigma_mat*init_est, 
                        init_est.t()*sigma_mat*temp_j ) ;
  }
  
 
  // update gamma, given nu estimates 
  arma::mat eye_d(d,d); eye_d.eye();
  arma::mat l2eye_p=l2_pen*eye(p,p);
  // arma::mat l2eye_p(p,p); l2eye_p.diag(l2_pen);
  // arma::mat bigeye_dp(d*p,d*p); bigeye_dp.diag(n*l2_pen);
  arma::mat cons_kron = kron( eye_d, sigma_mat + n*l2eye_p );

  
  arma::mat vvt0 = list_sum(nu_est0,nu_est0);
  // arma::mat inv_term01 = kron( vvt0 - eye_d, sigma_mat + l2eye_p );
  // arma::mat inv_term02 = kron( eye_d - vvt0 , l2eye_p );
  // arma::mat inv_term0 = inv_term01 + inv_term02 + cons_kron;
  arma::mat inv_term0 = kron( vvt0 - eye_d, sigma_mat + l2eye_p ) + 
    kron( eye_d - vvt0 , l2eye_p ) + cons_kron;
  // arma::mat inv_term0 = kron( vvt0 , sigma_mat ) + bigeye_dp;
  arma::mat ninv_term0 = reshape(sigma_mat*list_sum(grad_list, nu_est0 ), p*d, 1);
 
  // updating G is involves inverse a large, sparse kronecker matrix
  // where we need to define/declare A as sp_mat;
  // spsolve(x, A, b, "lapack" );  // use LAPACK  solver
  // spsolve(x, A, b, "superlu");  // use SuperLU solver

  // // opts only for superlu option, which is good for large, sparse, symmetric
  // superlu_opts opts;
  // opts.allow_ugly  = true;
  // opts.equilibrate = false; // normalize rows and columns to 1
  // opts.symmetric = true;

  
  arma::mat gamma_hat0 = solve(inv_term0, ninv_term0);
  arma::mat Gamma_hat0 = reshape(gamma_hat0, p, d);

  arma::mat loss0(1,1); loss0.zeros();
  for (j=0; j < n; j++) {
    arma::mat temp_j0 = grad_list[j];
    arma::mat temp_j1 = nu_est0[j];
    // double pen_norm = as_scalar( gamma_hat0.t()*gamma_hat0 );
    arma::mat diff_term0 = temp_j0 - sigma_mat*Gamma_hat0*temp_j1;

    loss0 =+ (reshape(diff_term0, 1, p*m)*reshape(diff_term0, p*m,1)/n +
      l2_pen*( gamma_hat0.t()*gamma_hat0 ) );
  }

  
  arma::uword iter; arma::uword k; Rcpp::List nu_est1(n);
  for (iter = 0; iter < max_iter; iter++ ) {

    // update nu_est
    // nu update, given initial Gamma 
    for (k = 0; k < n; k++) {
      arma::mat temp_k0= grad_list[k];
      nu_est1[k] = solve( Gamma_hat0.t()*sigma_mat*sigma_mat*Gamma_hat0 ,
                          Gamma_hat0.t()*sigma_mat*temp_k0 );
    }

    // update gamma, given nu estimates
    arma::mat vvt1 =  list_sum( nu_est1, nu_est1 );
    arma::mat inv_term1 = kron( vvt1 - eye_d, sigma_mat + l2eye_p ) + cons_kron +
      kron( eye_d - vvt1 , l2eye_p );
    // arma::mat inv_term1 = kron( vvt1 , sigma_mat ) + bigeye_dp;
    arma::mat ninv_term1 = reshape( sigma_mat*list_sum(grad_list, nu_est1 ), p*d,1);

    arma::mat gamma_hat1 = solve(inv_term1, ninv_term1);
    arma::mat Gamma_hat1 = reshape(gamma_hat1, p, d);

    arma::mat loss1(1,1); loss1.zeros();
    for (k=0; k < n; k++) {
      arma::mat temp_k1=grad_list[k];
      arma::mat temp_k2=nu_est1[k];
      arma::mat diff_term1 = temp_k1 - sigma_mat*Gamma_hat1*temp_k2;

      loss1 =+ reshape(diff_term1, 1, p*m)*reshape(diff_term1, p*m,1)/n +
        l2_pen* gamma_hat1.t()*gamma_hat1 ;
    }

    double loss_dist = as_scalar( loss0 - loss1 );

    if (print_iter) {
      Rcout << "Printing loss_dist, iter: " << loss_dist<< ", "  << iter<< "\n";
    }
    if( loss_dist < tol) {
      if (loss_dist > 0 ) Gamma_hat0 = Gamma_hat1; 
      // if loss is neg, next step is bad
      break;
    } else {
      loss0=loss1; Gamma_hat0 = Gamma_hat1;
    }

  }
  
  
  
  // return Rcpp::List::create(Rcpp::Named("l2eye_p")=l2eye_p,
  //                           Rcpp::Named("eye_d")=eye_d,
  //                           Rcpp::Named("cons_kron")=cons_kron,
  //                           Rcpp::Named("vvt0")=vvt0,
  //                           Rcpp::Named("inv_term01")=inv_term0,
  //                           Rcpp::Named("ninv_term0")=inv_term02);
  
  // return Rcpp::List::create(Rcpp::Named("nu_est1")=nu_est1);
  // Gamma_hat0
  return Rcpp::List::create(Rcpp::Named("Gamma_hat")=Gamma_hat0,
                            Rcpp::Named("nu_hat")=nu_est1);

}
