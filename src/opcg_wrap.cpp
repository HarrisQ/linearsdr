#include <RcppArmadillo.h>
#include "auxfun.h"
#include "cts_loss_funs.h"
#include "mn_loss_funs.h"

using namespace Rcpp;
using namespace arma;

  
//////////////////////////////////////////////////////////////////
///////////// CPP functions for OPCG-MADE Wrappers ///////////////
//////////////////////////////////////////////////////////////////

/***
 * My goal of implementing Rcpp is not to replace R code entirely,
 * but to augment the high level operations in R with faster computations
 * for the heavy operations in the function
 */ 

// For taking in a list and returning the candidate matrix for OPCG

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat list_mean(Rcpp::List est_list){
  /***
   * Used for creating the average the OUTER PRODUCT from a list
   */
  
  arma::uword n=est_list.length();
  
  // elements are all matrices and we want D*D.t()
  arma::mat test = est_list[0];
  arma::uword p=test.n_rows; //arma::uword d=test.n_cols;
  
  // calculating the mean of the gradients.
  // arma::mat mean_grad(p, d); mean_grad.zeros();
  // arma::uword l; // just need to declare this type once.
  // for (l = 0; l < n; l++ ) { 
  //   arma::mat grad_l = est_list[l];
  //   mean_grad += grad_l/n; 
  // }
  
  // Calculating the Outer Product
  // calculating the centered outer product
  // centering is not appropriate here since the gradient, as a function
  // is also random, and E g_hat (X) is not necessarily E g(X).
  
  arma::mat mean_outer_prod(p, p); mean_outer_prod.zeros();
  arma::uword j;
  for (j = 0; j < n; j++ ) {
    
    arma::mat Dj=est_list[j]; 
    // arma::mat cen_Dj= Dj - mean_grad;
    // arma::mat DD = cen_Dj*cen_Dj.t()/n;
    arma::mat DD = Dj*Dj.t()/n;
    mean_outer_prod += DD;
    
    // test = -wj(i)*tVij_I.t()*( y_datta.col(i) - mu_ij)/n;
  }
  return mean_outer_prod; //mean_grad; //
}


//////////////////////////////////////////////////////////////////////////
//////////// The Newton-Raphson Iteration for OPCG (and MADE) ////////////
//////////////////////////////////////////////////////////////////////////

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec aD_j_newton(arma::vec init,  
                      arma::mat vj, 
                      arma::mat y_datta, 
                      arma::vec wj, 
                      Rcpp::String link, 
                      arma::vec k,
                      double tol,
                      arma::uword max_iter,
                      bool test
                      ) {
  
  /***
   * The Newton Raphson Iteration for OPCG and MADE
   * 
   * Can either call the loss, score and info functions from
   * R for each iteration (calling from R is costly computation)
   * or load the arguments into the function and call within cpp, 
   * but this leads to messier code
   * 
   * arma::vec c, 
   arma::mat vj, 
   arma::mat y_datta, 
   arma::vec wj, 
   Rcpp::String link, 
   arma::vec k
   * 
   * 
   * 
   */
  
  arma::vec c_now = init; 
  
  arma::uword iter;
  arma::vec c_next;
  for (iter = 0; iter < max_iter; iter++ ) { 
    
    arma::mat nll_now(1,1); nll_now = mn_loss_j(c_now,vj,y_datta,wj,link,k); 
    
    arma::vec S_j = mn_score_j(c_now,vj,y_datta,wj,link,k);
    arma::mat I_j = mn_info_j(c_now,vj,y_datta,wj,link,k);
    
    c_next = c_now - pinv(I_j)*S_j;
    
    arma::mat nll_next(1,1); nll_next = mn_loss_j(c_next,vj,y_datta,wj,link,k); 
    
    double ll_dist; ll_dist = as_scalar( nll_now - nll_next);
    double eu_dist; eu_dist = euc_norm_cpp(c_now - c_next)/euc_norm_cpp(c_now);
    
    if (test) {
      // Rprintf("Printing: iter %iter, ll Dist %ll_dist, eu Dist %eu_dist ", 
      //         iter, ll_dist, eu_dist);
      Rcout << "Printing eu_Dist, ll_dist, iter: " << eu_dist  << ", " << 
        ll_dist<< ", "  << iter<< "\n";
    }
    if( ll_dist < tol) { 
      break;
    } else {
      c_now = c_next;
    }   

  }
  return c_next;
  
}

////////////////////////////////////////////////////////////////////////
/////////////////////// The Iteration for MADE /////////////////////////
////////////////////////////////////////////////////////////////////////


// Construction of the estimate

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat vecB_hat(arma::vec c0,
                   Rcpp::List score_list,
                   Rcpp::List info_list){
  
  arma::uword n=score_list.length();
  arma::uword p=c0.n_elem;
  
  // elements are all matrices and we want D*D.t()
  arma::vec score(p); score.zeros();
  arma::mat info(p,p); info.zeros();
  
  arma::uword j;
  for (j = 0; j < n; j++ ) {

    arma::vec score_j = score_list[j];
    arma::mat info_j = info_list[j];
    
    score += score_j;
    info += info_j;
  }
  
  // arma::mat reg_mat(p,p); reg_mat.ones();
  // + reg*reg_mat; + reg*c0
  arma::vec c1 = c0 - pinv(info )*(score );
  
  return c1;
}

//////////////////////////////////////////////////////////////////////////
//////////// The Conjugate Gradients for OPCG (and MADE) ////////////
//////////////////////////////////////////////////////////////////////////

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec aD_j_cg(arma::vec init,
                       arma::mat vj,
                       arma::mat y_datta,
                       arma::vec wj,
                       Rcpp::String link,
                       arma::vec k, 
                       Rcpp::List control_list,
                       bool test
) {
  
  /***
   * A Conjugate Gradient Algorithm for OPCG and MADE
   * 
   * Can either call the loss, score and info functions from
   * R for each iteration (calling from R is costly computation)
   * or load the arguments into the function and call within cpp, 
   * but this leads to messier code
   * 
   * arma::vec c, 
   arma::mat vj, 
   arma::mat y_datta, 
   arma::vec wj, 
   Rcpp::String link, 
   arma::vec k
   * 
   * 
   * 
   */
  
  // Set out controls for the algorithm
  
  double tol = control_list["tol_val"];
  arma::uword max_iter = control_list["max_iter"];
  
  // beta and sigma in backtracking rule (Bertsekas)
  // inital step size s_k that will be adjusted by backtracking 
  arma::vec s=control_list["init_stepsize"]; 
  double beta_bt=control_list["beta_bt"]; //backtracking beta
  double c_ag1=control_list["c_ag1"]; //parameter for Armijo condition, delta in Hybrid-DY
  double c_ag2=control_list["c_ag2"]; //parameter for second Armijo bound in nesterov 
  double c_wolfe=control_list["c_wolfe"]; //parameter for curvature condition, sigma in H-DY
  int max_iter_line=control_list["max_iter_line"];
  //int l2_pen=control_list["l2_pen"]; //penalty parameter for ridge regression
  
  // Step 0: Set the initial value, and compute the loss and score
  // Also set the initial p_0 to the gradient
  
  arma::uword m; m= y_datta.n_rows;
  arma::uword n; n= y_datta.n_cols;
  arma::uword p; p= vj.n_rows;
  arma::vec c_now = init; 
  
  arma::mat nll_now(1,1);  
  arma::vec grad_now; 
  nll_now = mn_loss_j(c_now,vj,y_datta,wj,link,k);
  grad_now = mn_score_j(c_now,vj,y_datta,wj,link,k);
  arma::vec p_now = -grad_now; 
  
  arma::vec c_next;
  arma::uword iter;
  // arma::vec things(7);
  
  for (iter = 0; iter < max_iter; iter++ ) { 
    double s_now = s(iter); 
    // Step 1: Line search 
    int m_ag=0;
    int m_cg;
    for (m_cg = 1; m_cg < max_iter_line; m_cg++ ) {
      
      // the Armijo rule
      // the computational effort is not just dependent on the m_cg search
      // it also depends on the step-size; large s_now means the m_cg will have
      // to search farther.
      // So we would like to have the smallest s_now to allow for fastest m_cg find. 
      double armijo_bound = as_scalar(c_ag1*pow(beta_bt,m_cg)*s_now*
                                      (p_now.t()*grad_now));
      
      // things(2)=armijo_bound; 
      
      // evaluation for armijo condition
      arma::vec c_search; c_search = c_now + pow(beta_bt,m_cg)*s_now*p_now;
      
      double suff_dec_ag;
      suff_dec_ag = as_scalar( mn_loss_j(c_search,vj,y_datta,wj,link,k) - 
        mn_loss_j(c_now,vj,y_datta,wj,link,k) );
      int armijo_cond = (suff_dec_ag <= armijo_bound);
      
      // things(3)=suff_dec_ag;
      
      int armijo_cond2=0;
      if (c_ag2 > 0) {
        // the second bound in armijo-goldstein in nesterovs intro to conv opt text
        double armijo_bound2 = as_scalar(c_ag2*pow(beta_bt,m_cg)*s_now*
                                         (p_now.t()*grad_now));
        // second sufficient descent bound uses the same suff_dec_ag
        armijo_cond2 =+ (suff_dec_ag >= armijo_bound2);

        // things(4)=armijo_bound2;
      } else if (c_ag2==0) {
        armijo_cond2=1;
      } 
      
      int wolfe_cond=0;
      if (c_wolfe > 0) {
        // the weak Wolfe condition
        double wolfe_bound = as_scalar(c_wolfe*(p_now.t()*grad_now));
        
        // evaluation for curvature in weak wolfe
        double curv_wolfe;
        curv_wolfe = as_scalar( p_now.t()*mn_score_j(c_search,vj,y_datta,wj,link,k) );
        
        wolfe_cond =+ (curv_wolfe <= wolfe_bound);
        
        // things(5)=wolfe_bound;things(6)= curv_wolfe;
      } else if (c_wolfe==0) {
        wolfe_cond=1;
      } 
      
      // things(1)=as_scalar(m_cg);
      
      
      if ( armijo_cond + armijo_cond2 + wolfe_cond == 3 ) { //  + armijo_cond2== 2
        m_ag = m_cg;
        break;
      }
    }  
    // things(0)=as_scalar(m_ag);
    
    double h_now = as_scalar( pow(beta_bt,m_ag)*s_now );
    c_next = c_now + h_now*p_now;
    
    // #Step 2a: Compute Loss;
    arma::mat nll_next(1,1); 
    nll_next=mn_loss_j(c_next,vj,y_datta,wj,link,k); 
    double nll_dist; nll_dist = as_scalar( nll_now - nll_next);
    
    if (test) {
      // Rprintf("Printing: iter %iter, ll Dist %ll_dist, eu Dist %eu_dist ", 
      //         iter, ll_dist, eu_dist);
      Rcout << "Printing nll_dist: " << nll_dist<< ", m_ag:"  << m_ag << ", iter:"  << iter << "\n ";
    }
    
    if ( nll_dist < tol ) {
      break;
    } else if( nll_dist.is_finite() ) {
      c_next=c_now;
      break;
      
    } else { 
      
      // #Step 2b: Compute gradient;
      arma::vec grad_next = mn_score_j(c_next,vj,y_datta,wj,link,k);
      
      // Step 3: Compute the coeffiecient
      // Fletcher-Reeves
      double beta_cg_fr = as_scalar( ( grad_next.t()*grad_next )/
                                      ( grad_now.t()*grad_now ) );  
      
      // Dai-Yuan
      double beta_cg_dy = as_scalar( ( grad_next.t()*grad_next )/
                                      ( p_now.t()*(grad_next-grad_now) ) );
      
      // Hestenes-Stiefel
      double beta_cg_hs= as_scalar( ( grad_next.t()*(grad_next-grad_now) )/
                                     ( p_now.t()*(grad_next-grad_now) ) );
      
      // Hybrid
      arma::vec beta2(2); beta2(0) = beta_cg_dy; beta2(1)=beta_cg_hs;
      arma::vec beta3(2); beta3(0) = 0; beta3(1)=min(beta2);
      double beta_cg_hybrid=max(beta3);
      
      // Step 4: Update p
      arma::vec p_next = -grad_next + beta_cg_dy*p_now;
      
      // Update all inputs for next iteration
      c_now=c_next;
      nll_now=nll_next;
      grad_now=grad_next;
      p_now=p_next;
    }
    
  } // end of cg iter  
  
  return c_next;//c_now;//things;
  
  

  
}

// A CG for MADE step
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec vecB_cg(arma::vec init,
                  arma::mat x_datta,
                  arma::mat y_datta,
                  double bw,
                  Rcpp::List ahat_list, 
                  Rcpp::List Dhat_list,
                  Rcpp::String link,
                  arma::vec k,
                  arma::mat r_mat, 
                  Rcpp::List control_list,
                  bool test
) {
  
  /***
   * A Conjugate Gradient Algorithm for OPCG and MADE
   * 
   * Can either call the loss, score and info functions from
   * R for each iteration (calling from R is costly computation)
   * or load the arguments into the function and call within cpp, 
   * but this leads to messier code
   * 
   * arma::vec c, 
   arma::mat vj, 
   arma::mat y_datta, 
   arma::vec wj, 
   Rcpp::String link, 
   arma::vec k
   * 
   * 
   * 
   */
  
  // Set out controls for the algorithm
  
  double tol = control_list["tol_val"];
  arma::uword max_iter = control_list["max_iter"];

  // beta and sigma in backtracking rule (Bertsekas)
  // inital step size s_k that will be adjusted by backtracking
  arma::vec s=control_list["init_stepsize"];
  double beta_bt=control_list["beta_bt"]; //backtracking beta
  double c_ag1=control_list["c_ag1"]; //parameter for Armijo condition, delta in Hybrid-DY
  double c_ag2=control_list["c_ag2"]; //parameter for Armijo condition, delta in Hybrid-DY
  double c_wolfe=control_list["c_wolfe"]; //parameter for curvature condition, sigma in H-DY
  int max_iter_line=control_list["max_iter_line"];  
  
  // Step 0: Set the initial value, and compute the loss and score
  // Also set the initial p_0 to the gradient
  
  arma::uword m; m= y_datta.n_rows;
  arma::uword n; n= y_datta.n_cols;
  arma::uword p; p= x_datta.n_rows;
  arma::vec c_now = init; 
  
  arma::mat nll_now(1,1);  
  arma::vec grad_now;
  nll_now = mn_loss_made(c_now,x_datta,y_datta,bw,ahat_list, Dhat_list,link,k,r_mat);
  grad_now = mn_score_made(c_now,x_datta,y_datta,bw,ahat_list, Dhat_list,link,k,r_mat);
  arma::vec p_now = -grad_now; 
  
  arma::vec c_next;
  arma::uword iter;
  for (iter = 0; iter < max_iter; iter++ ) { 
      
    // Need to pick initial stepsize for this iteration
    double s_now = s(iter);
    
    // Step 1: Line search
    int m_ag=0;
    int m_cg;
    for (m_cg = 1; m_cg < max_iter_line; m_cg++ ) {

      // the Armijo rule
      // the computational effort is not just dependent on the m_cg search
      // it also depends on the step-size; large s_now means the m_cg will have
      // to search farther.
      // So we would like to have the smallest s_now to allow for fastest m_cg find. 
      double armijo_bound = as_scalar(c_ag1*pow(beta_bt,m_cg)*s_now*
                                      (p_now.t()*grad_now));
      
      // evaluation for armijo condition
      arma::vec c_search; c_search = c_now + pow(beta_bt,m_cg)*s_now*p_now;
      
      double suff_dec_ag;
      suff_dec_ag = as_scalar( mn_loss_made(c_search,x_datta,y_datta,bw,ahat_list, Dhat_list,link,k,r_mat) 
                                 - mn_loss_made(c_now,x_datta,y_datta,bw,ahat_list, Dhat_list,link,k,r_mat) );
      int armijo_cond = (suff_dec_ag <= armijo_bound);
      
      // things(3)=suff_dec_ag;
      
      int armijo_cond2=0;
      if (c_ag2 > 0) {
        // the second bound in armijo-goldstein in nesterovs intro to conv opt text
        double armijo_bound2 = as_scalar(c_ag2*pow(beta_bt,m_cg)*s_now*
                                         (p_now.t()*grad_now));
        // second sufficient descent bound uses the same suff_dec_ag
        armijo_cond2 =+ (suff_dec_ag >= armijo_bound2);
        
        // things(4)=armijo_bound2;
      } else if (c_ag2==0) {
        armijo_cond2=1;
      } 
      
      int wolfe_cond=0;
      if (c_wolfe > 0) {
        // the weak Wolfe condition
        double wolfe_bound = as_scalar(c_wolfe*(p_now.t()*grad_now));
        
        // evaluation for curvature in weak wolfe
        double curv_wolfe;
        curv_wolfe = as_scalar( p_now.t()*mn_score_made(c_search,x_datta,y_datta,bw,ahat_list, Dhat_list,link,k,r_mat) );
        
        wolfe_cond =+ (curv_wolfe <= wolfe_bound);
        
        // things(5)=wolfe_bound;things(6)= curv_wolfe;
      } else if (c_wolfe==0) {
        wolfe_cond=1;
      } 
      
      // things(1)=as_scalar(m_cg);
      
      
      if ( armijo_cond + armijo_cond2 + wolfe_cond == 3 ) { //  + armijo_cond2== 2
        m_ag = m_cg;
        break;
      }
      
    } // End of line search loop
    
    
    double h_now = as_scalar( pow(beta_bt,m_ag)*s_now );
    c_next = c_now + h_now*p_now;
    
    // #Step 2a: Compute Loss;
    arma::mat nll_next(1,1);
    nll_next=mn_loss_made(c_next,x_datta,y_datta,bw,ahat_list, Dhat_list,link,k,r_mat) ;

    double nll_dist; nll_dist = as_scalar( nll_now - nll_next);
    
    if (test) {
      // Rprintf("Printing: iter %iter, ll Dist %ll_dist, eu Dist %eu_dist ", 
      //         iter, ll_dist, eu_dist);
      Rcout << "Printing nll_dist: " << nll_dist<< ", m_ag:"  << m_ag << ", iter:"  << iter << "\n ";
    }

    if( nll_dist < tol) {
      break;
    } else {

      // #Step 2b: Compute gradient;
      arma::vec grad_next = mn_score_made(c_next,x_datta,y_datta,
                                            bw,ahat_list, Dhat_list,link,k,r_mat) ;

      // Step 3: Compute the coeffiecient
      // Fletcher-Reeves
      double beta_cg_fr = as_scalar( ( grad_next.t()*grad_next )/
                                     ( grad_now.t()*grad_now ) );  
      
      // Dai-Yuan
      double beta_cg_dy = as_scalar( ( grad_next.t()*grad_next )/
                                     ( p_now.t()*(grad_next-grad_now) ) );
      
      // Hestenes-Stiefel
      double beta_cg_hs= as_scalar( ( grad_next.t()*(grad_next-grad_now) )/
                                    ( p_now.t()*(grad_next-grad_now) ) );
      
      // Hybird
      arma::vec beta2(2); beta2(0) = beta_cg_dy; beta2(1)=beta_cg_hs;
      arma::vec beta3(2); beta3(0) = 0; beta3(1)=min(beta2);
      double beta_cg_hybrid=max(beta3);
      
      // Step 4: Update p
      arma::vec p_next = -grad_next + beta_cg_dy*p_now;
      
      // Update all inputs for next iteration
      c_now=c_next;
      nll_now=nll_next;
      grad_now=grad_next;
      p_now=p_next;
    }
  } // End of CG Iterations
  
  return c_next;
  
}
