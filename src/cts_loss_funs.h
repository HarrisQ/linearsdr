#ifndef cts_loss_funs_H    
#define cts_loss_funs_H

using namespace Rcpp;
using namespace arma;

// ##################################################################
// ############## Loss Functions for OPCG and MADE ##################
// ##################################################################

// #######################################################################
// #                MADE for Continuous Response (i.e. MAVE)    
// #######################################################################
 
arma::mat nm_loss_made(arma::vec c, 
                       arma::mat x_matrix, 
                       arma::mat y_matrix, 
                       double bw, 
                       Rcpp::List ahat_list,
                       Rcpp::List Dhat_list,
                       arma::mat r_mat);  
  
arma::mat nm_loss_j_made(arma::vec c, 
                         arma::mat xj, 
                         arma::mat y_matrix, 
                         arma::vec wj, 
                         arma::vec ahat,
                         arma::mat Dhat);

arma::mat nm_score_j_made(arma::vec c, 
                          arma::mat xj, 
                          arma::mat y_matrix, 
                          arma::vec wj, 
                          arma::vec ahat,
                          arma::mat Dhat);

arma::mat nm_info_j_made(arma::vec c, 
                         arma::mat xj, 
                         arma::mat y_matrix, 
                         arma::vec wj, 
                         arma::vec ahat,
                         arma::mat Dhat);

#endif