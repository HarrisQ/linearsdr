#include <RcppArmadillo.h>
#include "auxfun.h"

using namespace Rcpp;
using namespace arma;

// #####################################################
// ############## Auxiliary Functions ##################
// #####################################################

// ##########  Centering function for a vector
//' @noRd
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector center_cpp(Rcpp::NumericVector x, 
                         Rcpp::Nullable<Rcpp::NumericVector> center = R_NilValue) {
  
  Rcpp::NumericVector x_cen;
  if (center.isNotNull() ) {
    Rcpp::NumericVector m; m=center;
    x_cen = x - m;
  } else {
    double m = mean(x); 
    x_cen = x - m;
  }
  return x_cen;
}

// ##########  Standardizing a Vector
//' @noRd
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector stand_vec_cpp(Rcpp::NumericVector x) {

  // Note that luckily, R has a mean function that will be called globally from R

  // Declare types
  // double m, s;

  // Compute Mean
  double m = mean(x); double s = sd(x);

  return (x - m)/s;
}


// ##########  Normalizing a Vector 
//' @noRd
//' @export
// [[Rcpp::export]]
arma::vec normalize_cpp(arma::colvec x) {
  //  double norm = arma::as_scalar( pow( x.t()*x , 0.5) );
  double sum_sq = as_scalar( x.t()*x );
  arma::vec w;
  if (sum_sq==0) {
     w = x ;
  } else {
    w = x/sum(x) ;
  }
  
  return w;
}

// ##########  Euclidean Norm
//' @noRd
//' @export
// [[Rcpp::export]]
double euc_norm_cpp(arma::vec x) {
  double norm = arma::as_scalar( pow( x.t()*x , 0.5) );
  return norm;
}

// For taking in a list and returning the candidate matrix for OPCG

//' @noRd
//' @export
// [[Rcpp::export]]
arma::mat list_sum(Rcpp::List listA, Rcpp::List listB){
  /***
   * Used for creating the sum of product from list A and list B, 
   * where the elements in list B are transposed
   */
  
  arma::uword n=listA.length();
  
  // elements are all matrices and we want D*D.t()
  arma::mat testA = listA[0]; arma::mat testB = listB[0];
  arma::uword p=testA.n_rows; //arma::uword d=testA.n_cols;
  arma::uword m=testB.n_rows; 
  
  arma::mat out(p, m); out.zeros();
  arma::uword j;
  for (j = 0; j < n; j++ ) {
    
    arma::mat Aj=listA[j];
    arma::mat Bj=listB[j];
    arma::mat ABj = Aj*Bj.t() ;
    out += ABj; 
  }
  return out;  
}


//////////////////////////////////////////////////////////////////
// Functions related to Eigen Decomp and roots and powers of matrices
/////////////////////////////////////////////////////////////////

// ##########  Matrix Power  

//' General Matrix Power Function
//' 
//' This is an internal function used for computing matrix powers
//' 
//' @param a a square matrix 
//' @param alpha power of matrix to be taken 
//' @param lead number of leading eigenvalues to consider; default is all
//' @param ignore decimal place to ignore; default is 10^(-15)
//'
//' @return a matrix to the power alpha
//'
//' @keywords internal
//' 
//' @noRd
//' @export
// [[Rcpp::export]]
arma::mat matpower_cpp(arma::mat a, double alpha,
                       Rcpp::Nullable<double> lead = R_NilValue,
                       Rcpp::Nullable<double> ignore = R_NilValue){
  
  // Declare Types (in Armadillo)
  arma::vec eigval; arma::mat eigvec; 
  
  // Run EigenDecomp on matrix a, symmetrized
  eig_sym(eigval, eigvec, 0.5*( a + a.t() ) ) ;
  eigvec=fliplr(eigvec); eigval=sort(eigval, "descend"); 
  
  
  arma::uword p; p = eigval.n_elem; 
  
  // eigval.elem(find(eigval<1)).fill(0);
  //arma::uvec ids(p);  
  
  
  if ( ignore.isNotNull() ) { 
    
    // Need to assign the ignore value to a double;
    double ignore_val = Rcpp::as<double>(ignore);
    
    eigval.elem(find(eigval<ignore_val)).fill(0);
    
  } 
  
  if ( lead.isNotNull() ) { 
    
    // Need to assign the ignore value to a double;
    uword trunc_index = Rcpp::as<arma::uword>(lead);
    
    eigval.subvec(trunc_index,p-1).fill(0);
  }
  
  
  // Return Result
  return eigvec*diagmat(pow(eigval, alpha))*eigvec.t();
}

// ##########  Matrix Centering  
// centers the matrix
//' @noRd
//' @export
// [[Rcpp::export]]
arma::mat matcenter_cpp(arma::mat x_matrix, 
                        Rcpp::Nullable<unsigned int> index = R_NilValue,
                        Rcpp::Nullable<Rcpp::NumericVector> x0 = R_NilValue) {
  
  // recall cpp centering starts at 0; 
  if ( index.isNotNull() ) { 
    arma::uword j; j = Rcpp::as<arma::uword>(index) - 1;
    x_matrix.each_col() -= x_matrix.col(j); 
  }
  if ( x0.isNotNull() ) { 
    arma::vec x00; x00 = Rcpp::as<arma::vec>(x0);
    x_matrix.each_col() -= x00;
  }
  return x_matrix;
}

// For doing eigen decomp of sym pos-def mat
//' @noRd
//' @export
// [[Rcpp::export]]
Rcpp::List eigen_cpp(arma::mat A){
  /***
   * Used for conducting eigen decomp
   */
  
  // Declare Types (in Armadillo)
  arma::vec eigval; arma::mat eigvec;
  
  // Run EigenDecomp on matrix A; we symmetrize A to ensure 
  // positive-definiteness
  eig_sym(eigval, eigvec, A = 0.5*( A + A.t() ) ) ;
  
  // Return Result
  return Rcpp::List::create(Rcpp::Named("values")=sort(eigval, "descend"),
                            Rcpp::Named("vectors")=fliplr(eigvec) ); 
}

//' @noRd
//' @export
// [[Rcpp::export]]
Rcpp::List gev_cpp(arma::mat A, 
                   arma::mat B){ 
  /***
   * B is the metric, A is the matrix we want the eigenvalues of
   */
  
  // Declare Types (in Armadillo)
  arma::cx_vec eigval; 
  arma::cx_mat  eigvec; // this is the left eigenvector as well
  arma::cx_mat  leigvec;
  arma::cx_mat  reigvec;
  
  // Solving Generalized Eigen problem 
  // leigvec, reigvec, - in place of eigvec if we want left and right eigvec 
  // separately
  // eig_pair( eigval, eigvec, M, G );
  eig_pair(eigval, leigvec, reigvec, A, B);
  
  // Return Result
  return Rcpp::List::create(Rcpp::Named("values")=sort(eigval, "descend"),
                            Rcpp::Named("vectors")=fliplr(leigvec) ); 
}


// For doing inverse of sym pos-def mat
//' @noRd
//' @export
// [[Rcpp::export]]
arma::mat inv_sympd_cpp(arma::mat A){ 
  
  arma::mat C = inv_sympd(A);
  
  // Return Result
  return C; 
}

// For doing sqrt of sym pos-def mat
//' @noRd
//' @export
// [[Rcpp::export]]
arma::mat sqrtmat_cpp(arma::mat A){ 
  
  arma::mat C = sqrtmat_sympd(A);
  
  // Return Result
  return C; 
}


// For doing chol decom of mat
//' @noRd
//' @export
// [[Rcpp::export]]
arma::mat chol_cpp(arma::mat A ){ 
  
  // chol layout is default to upper triangular, so C.t()*C = A
  arma::mat C = chol(A);
  
  // Return Result
  return C; 
}

// For doing solving linear systems
// Standardizing a random matrix or variable is an example of this
//' @noRd
//' @export
// [[Rcpp::export]]
arma::mat solve_cpp(arma::mat A,
                    arma::mat B){ 
  
  // C = inv(A)*B
  arma::mat C = solve(A, B);
  
  // Return Result
  return C; 
}

//////////////////////////////////////////
// ########## Gaussian kernel Weight Function 
///////////////////////////////////////////
//' @noRd
//' @export
// [[Rcpp::export]]
arma::vec gauss_kern_cpp(arma::mat centered_data, double bw ) {
  // do not need the exact pdf since the bandwidth h is proporitional anywas
  // can do what Bing does and just compute exp(-0.5*(X-X_0).t()*(X-X_0) ) 
  // But computationally inefficient compared to just using pdf norm
  
  // arma::mat Wj; Wj = log(normpdf(centered_data, 0, bw));
  // arma::rowvec WW; WW = sum(Wj); //sum adds over rows
  // //since WW is a rowvec, need .t() to make it colvec, but still not vec
  // Wj = normalise(exp(WW), 1).t();
  
  arma::mat XtX = centered_data.t()*centered_data;
  double cj = -.5*pow(bw,-2);
  arma::vec XtX_diag = cj* diagvec(XtX);
  arma::vec Wj;
  Wj = normalise(exp(XtX_diag), 1) ;
  
  return Wj; 
}


// #####################################################
// Weighted Least Squares
// #####################################################
//' @noRd
//' @export
// [[Rcpp::export]]
arma::mat wls_cpp(arma::mat x_matrix, 
                  arma::mat y_matrix,
                  arma::vec weights,
                  double reg) {
  /***
   * matrices are entered as p x n and m x n
   * If centering is necessary, the X matrix is already centered at the jth 
   * observation; if intercept is necessary, it is already done
   */
  
  // create identity matrix
  arma::uword p; p=x_matrix.n_rows;
  arma::mat reg_mat(p,p); reg_mat.eye();
  
  arma::mat W; W = diagmat( weights );
  // Solve is equivalent to X*Beta = Y, where Beta is unknown
  
  // Enforce symmetry of the weighted design matrix
  arma::mat xwx; xwx = x_matrix*W*x_matrix.t();
  arma::mat sX; sX = 0.5*(xwx + xwx.t() ) + reg*reg_mat;
  
  return arma::solve( sX, x_matrix*W*y_matrix.t() );
  
}
 