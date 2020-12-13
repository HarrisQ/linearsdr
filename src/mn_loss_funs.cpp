#include <RcppArmadillo.h>
#include "auxfun.h" 
#include "mn_loss_funs.h" 

using namespace Rcpp;
using namespace arma;

// ##################################################################
// #                    Loss Functions for OPCG and MADE 
// ##################################################################

// ################################################################## 
// #                 Empirical Logit and Culmit Transforms
// ################################################################## 

// # Multinomial Y to Multivariate Y ####
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mnY_to_mvY(arma::mat mn_y,
                     arma::vec m_classes,
                     Rcpp::String ytype) {
  /***
   * Takes in a 1 x n vector with multinomial entries
   * Returns a m x n matrix, with m being the categories
   */
  
  arma::uword n = mn_y.n_cols;
  arma::uword m = m_classes.n_elem;
  
  arma::mat mv_Y(m,n); //mv_Y.zeros(m,n);  //sapply( seq_len(n), mv_Y_i);
  
  if (ytype=="multinomial") {
    // Writing the For loop instead of sapply.
    arma::uword i;
    for (i = 0; i < n; i++ ) {
      
      arma::mat Y_i(m,1) ; Y_i.zeros();
      arma::uvec ids = find(m_classes==mn_y(i));
      Y_i.rows(ids).fill(1);
      
      mv_Y.col(i) = Y_i;
    }
    
  } if (ytype=="ordinal") {
    arma::uword i;
    for (i = 0; i < n; i++ ) {
      
      arma::vec Y_i ; Y_i.zeros(m);
      arma::uvec ids = find(m_classes<=mn_y(i)); 
      Y_i.elem(ids).fill(1);
      arma::vec cum_Y_i = cumsum(Y_i);
      
      mv_Y.col(i) = cum_Y_i/max(cum_Y_i);
    }
  }
  return mv_Y;
};

// # Empirical Logit Transform ####
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat emp_logit(arma::mat y_matrix,
                    arma::vec k_vec,
                    double tune) {
  
  /***
   * y_matrix is m x n
   */
  
  arma::uword n = y_matrix.n_cols;
  arma::uword m = y_matrix.n_rows;
  
  arma::vec k = k_vec;
  
  // set tuning 
  
  arma::mat tildeY0 = y_matrix;
  tildeY0.elem( find( tildeY0 == 0) ).fill(tune);  
  
  // # So that the numerator in log(p/(1-p)) is not zero
  
  arma::mat emp_logit(m,n); //mv_Y.zeros(m,n);  //sapply( seq_len(n), mv_Y_i);
  
  // Writing the For loop instead of sapply.
  arma::uword i;
  for (i = 0; i < n; i++ ) {
    
    arma::vec logit_y_i = tildeY0.col(i)/k(i);
    double d1 = sum(logit_y_i);
    double d;
    if ( d1 >= 1 ){
      d = d1 + tune; 
    } else {
      d = d1; 
    } 
    
    emp_logit.col(i) = log( abs( logit_y_i/(1-d)) );
    // emp_logit.col(i) = d*arma::ones(m);
    
  } 
  
  // recall that cpp starts indexing at 0;
  return emp_logit; 
} ;

// # Empirical Cumulative Transform ####
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat emp_culmit(arma::mat y_matrix,
                     arma::vec k_vec,
                     double tune) {
  
  /***
   * y_matrix is m x n
   */
  
  arma::uword n = y_matrix.n_cols;
  arma::rowvec ones_vec; ones_vec.ones(n);
  arma::vec k = k_vec;
  
  
  // re-fill matrix with last row
  arma::mat tildeY0 ; tildeY0 = join_cols(y_matrix,ones_vec);
  arma::uword m = tildeY0.n_rows;
  arma::mat I; I.eye(m, m);
  
  
  // create matrices with ones; they are difference matrices;
  arma::mat C(m,m - 1); C.eye(); C.diag(1).fill(-1);
  arma::mat D(m,m - 1); D.eye(); D.diag().fill(-1); D.diag(-1).fill(1); 
  
  // applying difference matrices to Y and tuning;
  arma::mat C_Y = C.t()*tildeY0; arma::mat D_Y = D.t()*tildeY0;
  C_Y.elem( find(C_Y == 0) ).fill(tune); 
  D_Y.elem( find(D_Y == 0) ).fill(tune);
  
  // # So that the numerator in log(p/(1-p)) is not zero
  arma::mat emp_culmit(m-1,n);  
  
  // Writing the For loop instead of sapply.
  arma::uword i;
  for (i = 0; i < n; i++ ) {
    emp_culmit.col(i) = k(i)*log( diagmat(1/D_Y.col(i))*C_Y.col(i) );
    
  }
  
  return emp_culmit;   
  
  
};

// #######################################################################
// #                  Multinomial GLM dot.b function 
// #######################################################################

// Derivative of the conditional Mean for the MultiNomial
// Called within par_obj_i, which is within par_next_cpp

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec dot_b_multinom(arma::vec lin_can_par, int k_i, String link ){
  
  /***
   *  Matrix return with columns as the b_i for observation/individual i;
   *  \th_i is the canonical parameter for obs i;
   *  \th_i is p dimensional
   *  th.i = as.vector(tZij.I%*%c.now)
   */ 
  
  arma::uword m=lin_can_par.n_rows;
  arma::vec dot_b;
  
  if ( link == "expit") {
    
    arma::vec e_lcp=exp(lin_can_par);
    // arma::vec ones_vec(m); ones_vec.ones();
    double dem; dem = 1 + sum(e_lcp) ;
    
    arma::vec pi_i=e_lcp/dem;
    arma::uvec ids = find(pi_i < 0); // Find indices
    pi_i.elem(ids).fill(0);   // Assign value '0' to ids identified by condition
    
    dot_b = k_i*pi_i;
    
  } else if (link == "culmit") {
    
    // creating upper/lower triangular matrices;
    arma::mat A; A.ones(m,m);
    arma::mat U=trimatu(A); arma::mat L=trimatl(A);
    
    arma::vec tau_i = L*exp( U*lin_can_par/k_i );
    //cpp index starts at 0, so m-1 is the mth entry
    double dem; dem = 1 + tau_i(m-1); 
    
    dot_b = tau_i/dem;
  }
  return dot_b;
};

// These need to come after dot_b because they use it.  
// # b.expit
double b_expit(arma::vec lin_can_par, int k_i) {
  return k_i*log( 1 + sum( exp( lin_can_par ) ) );
};

// # b.culmit
double b_culmit(arma::vec lin_can_par, int k_i) {
  arma::uword m=lin_can_par.n_rows;
  arma::vec mu_i=dot_b_multinom( lin_can_par, k_i, "culmit"); 
  return -k_i*log(  (1 - mu_i(m-1))  ) ;
};  

// #######################################################################
// #           Multinomial GLM loss function (i.e. Deviance)
// #######################################################################
// 
// ### The loss, score and info functions for OPCG ----

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mn_loss_j(arma::vec c, 
                    arma::mat vj, 
                    arma::mat y_datta, 
                    arma::vec wj, 
                    Rcpp::String link, 
                    arma::vec k) {
  // # c <- c.j.ls
  arma::uword n=y_datta.n_cols;
  arma::uword m=y_datta.n_rows;
  arma::mat I(m,m); I.eye();
  
  arma::mat mean_nll_j(1,1); mean_nll_j.zeros();   
  // arma::mat test;
  
  if (link=="expit"){
    
    // Writing the For loop instead of sapply.
    arma::uword i;
    for (i = 0; i < n; i++ ) {
      
      arma::mat tVij_I=kron( (vj.col(i)).t(),I);
      arma::vec lcp=tVij_I*c;
      
      mean_nll_j += -wj(i)*( lcp.t()*( y_datta.col(i)) - b_expit(lcp, k(i) ) )/n;
      
      // test = -wj(i)*( lcp.t()*( y_datta.col(i)) - b_expit(lcp, k(i) ) );
    } 
    
  } else if (link=="culmit") {
    
    // Writing the For loop instead of sapply.
    arma::uword i;
    for (i = 0; i < n; i++ ) {
      
      arma::mat tVij_I=kron( (vj.col(i)).t(),I);
      arma::vec lcp=tVij_I*c;
      
      mean_nll_j += -wj(i)*( lcp.t()*y_datta.col(i) - b_expit(lcp, k(i) ) )/n;
    } 
  }
  
  return mean_nll_j;
  
} ;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mn_score_j(arma::vec c,
                     arma::mat vj,
                     arma::mat y_datta,
                     arma::vec wj,
                     Rcpp::String link,
                     arma::vec k) {
  // # c <- c.j.ls
  arma::uword n=y_datta.n_cols;
  arma::uword m=y_datta.n_rows;
  arma::uword pm=c.n_elem;
  arma::mat I(m,m); I.eye();
  
  arma::mat mean_score_j(pm,1); mean_score_j.zeros();
  // arma::mat test;
  
  // Link only matters for the mean, the form of the score is
  // always the same;
  // Writing the For loop instead of sapply.
  arma::uword i;
  for (i = 0; i < n; i++ ) {
    
    arma::mat tVij_I=kron( (vj.col(i)).t(),I);
    arma::vec lcp=tVij_I*c;
    arma::vec mu_ij = dot_b_multinom(lcp, k(i), link);
    
    mean_score_j += -wj(i)*tVij_I.t()*( y_datta.col(i) - mu_ij)/n;
    
    // test = -wj(i)*tVij_I.t()*( y_datta.col(i) - mu_ij)/n;
  }
  
  return mean_score_j;
  
};

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mn_info_j(arma::vec c,
                    arma::mat vj,
                    arma::mat y_datta,
                    arma::vec wj,
                    Rcpp::String link,
                    arma::vec k) {
  // # c <- c.j.ls
  arma::uword n=y_datta.n_cols;
  arma::uword m=y_datta.n_rows;
  arma::uword pm=c.n_elem;
  arma::mat I(m,m); I.eye();
  
  arma::mat mean_info_j(pm,pm); mean_info_j.zeros();
  // arma::mat test;
  
  if (link=="expit"){
    
    // Writing the For loop instead of sapply.
    arma::uword i;
    for (i = 0; i < n; i++ ) {
      
      arma::mat tVij_I=kron( (vj.col(i)).t(),I);
      arma::vec lcp=tVij_I*c;
      arma::vec mu_ij = dot_b_multinom(lcp, k(i), link);
      
      mean_info_j += wj(i)*k(i)*tVij_I.t()*( diagmat(mu_ij) - mu_ij*mu_ij.t())*tVij_I/n;
      
      // test = -wj(i)*tVij_I.t()*( y_datta.col(i) - mu_ij)/n;
    }
    
  } else if (link=="culmit") {
    // Writing the For loop instead of sapply.
    arma::uword i;
    for (i = 0; i < n; i++ ) {
      
      arma::mat tVij_I=kron( (vj.col(i)).t(),I);
      arma::vec lcp=tVij_I*c;
      arma::vec mu_ij = dot_b_multinom(lcp, k(i), link);
      
      arma::vec v_m(m); v_m.ones(); //vec of ones
      arma::mat E; E = v_m*mu_ij.t();
      arma::mat E_syml = symmatl(E);
      
      mean_info_j += wj(i)*tVij_I.t()*(E_syml - mu_ij*mu_ij.t())*tVij_I/(k(i)*n);
      
      // test = -wj(i)*tVij_I.t()*( y_datta.col(i) - mu_ij)/n; 
    }
    
    
  }
  
  arma::mat I_j; I_j = (mean_info_j + mean_info_j.t())/2;
  return I_j;
  
};

// ############# The loss, score and info functions for MADE #################

// A local Loss function for MADE 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mn_loss_j_made(arma::vec c, 
                         arma::mat xj, 
                         arma::mat y_matrix, 
                         arma::vec wj, 
                         arma::vec ahat,
                         arma::mat Dhat,
                         Rcpp::String link, 
                         arma::vec k) {
  
  
  arma::uword n=y_matrix.n_cols;
  arma::uword m=y_matrix.n_rows + 1;
  arma::mat I(m,m); I.eye();
  
  arma::mat mean_nll_j(1,1); mean_nll_j.zeros();   
  // arma::mat test;
  
  if (link=="expit"){
    
    // Writing the For loop instead of sapply.
    arma::uword i;
    
    for (i = 0; i < n; i++ ) {
      
      arma::mat tXij_tD=kron( (xj.col(i)).t(),Dhat.t());
      arma::vec lcp=tXij_tD*c;
      
      mean_nll_j += -wj(i)*( (ahat + lcp).t()* y_matrix.col(i) - 
        b_expit((ahat + lcp), k(i) ) )/n;
      // test = -wj(i)*( lcp.t()*( y_datta.col(i)) - b_expit(lcp, k(i) ) );
    } 
    
  } else if (link=="culmit") {
    
    // Writing the For loop instead of sapply.
    arma::uword i;
    for (i = 0; i < n; i++ ) {
      
      arma::mat tXij_tD=kron( (xj.col(i)).t(),Dhat.t());
      arma::vec lcp=tXij_tD*c;
      
      mean_nll_j += -wj(i)*( (ahat + lcp).t()* y_matrix.col(i) - 
        b_culmit((ahat + lcp), k(i) ) )/n;
      
    } 
  }
  
  return mean_nll_j;
  
} ; 


// A Loss function for whole MADE 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mn_loss_made(arma::vec c, 
                       arma::mat x_matrix, 
                       arma::mat y_matrix, 
                       double bw,
                       Rcpp::List ahat_list,
                       Rcpp::List Dhat_list,
                       Rcpp::String link, 
                       arma::vec k,
                       arma::mat r_mat) {
  
  
  arma::uword n=y_matrix.n_cols; //arma::uword m=y_matrix.n_rows;  
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
    
    // not really a mean
    mean_nll += mn_loss_j_made(c, xj, y_matrix, wj, ahat, Dhat, link, k)/n;
    // test = -wj(i)*( lcp.t()*( y_datta.col(i)) - b_expit(lcp, k(i) ) );
  } 
  
  
  return mean_nll;
  
} ;

// A Local Score function for MADE 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mn_score_j_made(arma::vec c, 
                          arma::mat xj, 
                          arma::mat y_matrix, 
                          arma::vec wj, 
                          arma::vec ahat,
                          arma::mat Dhat,
                          Rcpp::String link, 
                          arma::vec k) {
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
    arma::vec mu_ij = dot_b_multinom( ( ahat + lcp), k(i), link);
    
   
    mean_score_j += -wj(i)*tXij_tD.t()*( y_matrix.col(i) - mu_ij)/n;
    
    // test = -wj(i)*tVij_I.t()*( y_datta.col(i) - mu_ij)/n;
  }
  
  return mean_score_j;
  
};


// A Score function for whole MADE 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mn_score_made(arma::vec c, 
                        arma::mat x_matrix, 
                        arma::mat y_matrix, 
                        double bw,
                        Rcpp::List ahat_list,
                        Rcpp::List Dhat_list,
                        Rcpp::String link, 
                        arma::vec k,
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
    
    // not really a mean, just a sum
    mean_score += mn_score_j_made(c, xj, y_matrix, wj, ahat, Dhat, link, k)/n;
    // test = -wj(i)*( lcp.t()*( y_datta.col(i)) - b_expit(lcp, k(i) ) );
    
  } 
  
  
  return mean_score;
  
} ;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mn_info_j_made(arma::vec c, 
                         arma::mat xj, 
                         arma::mat y_matrix, 
                         arma::vec wj, 
                         arma::vec ahat,
                         arma::mat Dhat,
                         Rcpp::String link, 
                         arma::vec k) {
  // # c <- c.j.ls
  arma::uword n=y_matrix.n_cols;
  arma::uword m=y_matrix.n_rows;
  arma::uword pm=c.n_elem;
  arma::mat I(m,m); I.eye();
  
  arma::mat mean_info_j(pm,pm); mean_info_j.zeros();
  // arma::mat test;
  
  if (link=="expit"){
    
    // Writing the For loop instead of sapply.
    arma::uword i;
    for (i = 0; i < n; i++ ) {
      
      arma::mat tXij_tD=kron( (xj.col(i)).t(),Dhat.t());
      arma::vec lcp=tXij_tD*c;
      arma::vec mu_ij = dot_b_multinom( ( ahat + lcp), k(i), link);
      
      mean_info_j += wj(i)*k(i)*tXij_tD.t()*( diagmat(mu_ij) - mu_ij*mu_ij.t())*tXij_tD/n;
      
      // test = -wj(i)*tVij_I.t()*( y_datta.col(i) - mu_ij)/n;
    }
    
  } else if (link=="culmit") {
    // Writing the For loop instead of sapply.
    arma::uword i;
    for (i = 0; i < n; i++ ) {
      
      arma::mat tXij_tD=kron( (xj.col(i)).t(),Dhat.t());
      arma::vec lcp=tXij_tD*c;
      arma::vec mu_ij = dot_b_multinom( ( ahat + lcp), k(i), link);
      
      arma::vec v_m(m); v_m.ones(); //vec of ones
      arma::mat E; E = v_m*mu_ij.t();
      arma::mat E_syml = symmatl(E);
      
      mean_info_j += wj(i)*tXij_tD.t()*(E_syml - mu_ij*mu_ij.t())*tXij_tD/(k(i)*n);
      
      // test = -wj(i)*tVij_I.t()*( y_datta.col(i) - mu_ij)/n;  
    }
    
    
  }
  
  arma::mat I_j; I_j = (mean_info_j + mean_info_j.t())/2;
  return I_j;
  
};