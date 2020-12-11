#ifndef mn_loss_funs_H    
#define mn_loss_funs_H

using namespace Rcpp;
using namespace arma;

// ##################################################################
// ############## Loss Functions for OPCG and MADE ##################
// ##################################################################
 
// #######################################################################
// #                  Empirical Logit and Culmit Transforms
// #######################################################################

// # Multinomial Y to Multivariate Y ####
 
arma::mat mnY_to_mvY(arma::mat mn_y,
                     arma::vec m_classes,
                     Rcpp::String ytype);  

// # Empirical Logit Transform ####
 
arma::mat emp_logit(arma::mat y_matrix,
                    arma::vec k_vec,
                    double tune);  

// # Empirical Cumulative Transform ####
 
arma::mat emp_culmit(arma::mat y_matrix,
                     arma::vec k_vec,
                     double tune);  

// #######################################################################
// #                  Multinomial GLM dot.b function 
// #######################################################################

// Derivative of the conditional Mean for the MultiNomial
// Called within par_obj_i, which is within par_next_cpp

 
arma::vec dot_b_multinom(arma::vec lin_can_par, int k_i, String link ); 

// These need to come after dot_b because they use it.  
// # b.expit
double b_expit(arma::vec lin_can_par, int k_i);  

// # b.culmit
double b_culmit(arma::vec lin_can_par, int k_i); 

// #######################################################################
// #           Multinomial GLM loss function (i.e. Deviance)
// #######################################################################
// 
// ### The loss, score and info functions for OPCG ----
 
arma::mat mn_loss_j(arma::vec c, 
                    arma::mat vj, 
                    arma::mat y_datta, 
                    arma::vec wj, 
                    Rcpp::String link, 
                    arma::vec k);  

  
arma::mat mn_score_j(arma::vec c,
                     arma::mat vj,
                     arma::mat y_datta,
                     arma::vec wj,
                     Rcpp::String link,
                     arma::vec k);  

 
arma::mat mn_info_j(arma::vec c,
                    arma::mat vj,
                    arma::mat y_datta,
                    arma::vec wj,
                    Rcpp::String link,
                    arma::vec k);  


// ################## Loss functions for MADE ################# 

arma::mat mn_loss_made(arma::vec c, 
                       arma::mat x_matrix, 
                       arma::mat y_matrix, 
                       double bw,
                       Rcpp::List ahat_list,
                       Rcpp::List Dhat_list,
                       Rcpp::String link, 
                       arma::vec k,
                       arma::mat r_mat); 
  
arma::mat mn_loss_j_made(arma::vec c, 
                         arma::mat xj, 
                         arma::mat y_matrix, 
                         arma::vec wj, 
                         arma::vec ahat,
                         arma::mat Dhat,
                         Rcpp::String link, 
                         arma::vec k);

 
arma::mat mn_score_made(arma::vec c, 
                        arma::mat x_matrix, 
                        arma::mat y_matrix, 
                        double bw,
                        Rcpp::List ahat_list,
                        Rcpp::List Dhat_list,
                        Rcpp::String link, 
                        arma::vec k,
                        arma::mat r_mat);
  
arma::mat mn_score_j_made(arma::vec c, 
                          arma::mat xj, 
                          arma::mat y_matrix, 
                          arma::vec wj, 
                          arma::vec ahat,
                          arma::mat Dhat,
                          Rcpp::String link, 
                          arma::vec k);

arma::mat mn_info_j_made(arma::vec c, 
                         arma::mat xj, 
                         arma::mat y_matrix, 
                         arma::vec wj, 
                         arma::vec ahat,
                         arma::mat Dhat,
                         Rcpp::String link, 
                         arma::vec k);


#endif