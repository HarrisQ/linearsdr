#ifndef auxfun_H    
#define auxfun_H

using namespace Rcpp;
using namespace arma;

// #####################################################
// ############## Auxiliary Functions ##################
// #####################################################

// ##########  Centering function for a vector 

Rcpp::NumericVector center_cpp(Rcpp::NumericVector x, 
                         Rcpp::Nullable<Rcpp::NumericVector> center); 

// ##########  Standardizing a Vector
 
Rcpp::NumericVector stand_vec_cpp(Rcpp::NumericVector x);

// ##########  Normalizing a Vector
 
arma::vec normalize_cpp(arma::vec x);

// ##########  Euclidean Norm 
 
double euc_norm_cpp(arma::vec x);

// For taking in a list and returning the candidate matrix for OPCG

arma::mat list_sum(Rcpp::List listA, Rcpp::List listB);


//////////////////////////////////////////////////////////////////
// Functions related to Eigen Decomp and roots and powers of matrices
/////////////////////////////////////////////////////////////////

    
// ##########  Matrix Power  
// the matpower function  
 
arma::mat matpower_cpp(arma::mat A, double alpha); 

// ##########  Matrix Centering  
// centers the matrix 
 
arma::mat matcenter_cpp(arma::mat x_matrix, 
                        Rcpp::Nullable<unsigned int> index,
                        Rcpp::Nullable<Rcpp::NumericVector> x0);  



// Generalized Eigen Value Problem
Rcpp::List gev_cpp(arma::mat A, arma::mat B);

// For doing inverse of sym pos-def mat 
arma::mat inv_sympd_cpp(arma::mat A);

// For doing sqrt of sym pos-def mat
arma::mat sqrtmat_cpp(arma::mat A);


// For doing chol decom of mat 
arma::mat chol_cpp(arma::mat A );

// For doing solving linear systems
// Standardizing a random matrix or variable is an example of this
arma::mat solve_cpp(arma::mat A, arma::mat B);

// ############################################################
// ########## Gaussian kernel Weight Function 
 
arma::vec gauss_kern_cpp(arma::mat centered_data, double bw ); 

// #####################################################
// Weighted Least Squares
// #####################################################
 
arma::mat wls_cpp(arma::mat x_matrix, arma::mat y_matrix,
                  arma::vec weights); 


#endif