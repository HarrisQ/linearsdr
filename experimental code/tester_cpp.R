
#### cg ----

sourceCpp(code='
  #include <RcppArmadillo.h>
  using namespace Rcpp;
  using namespace arma;



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec dot_b_multinom(arma::vec lin_can_par, int k_i, String link ){
  
  arma::uword m=lin_can_par.n_rows;
  arma::vec dot_b;
  
  if ( link == "expit") {
    
    arma::vec e_lcp=exp(lin_can_par);
    // arma::vec ones_vec(m); ones_vec.ones();
    double dem; dem = 1 + sum(e_lcp) ;
    
    arma::vec pi_i=e_lcp/dem;
    arma::uvec ids = find(pi_i < 0);  
    pi_i.elem(ids).fill(0);    
    
    dot_b = k_i*pi_i;
    
    return dot_b;
    
  } else if (link == "culmit") {
    
    // creating upper/lower triangular matrices;
    arma::mat A; A.ones(m,m);
    //arma::mat U=trimatu(A); 
    arma::mat L=trimatl(A);
    
    // creating Permutation and Differencing Matrices, P, Q
    arma::mat I(m,m); I.eye();
    arma::mat P(m,m); P.zeros(); P.cols(0,m-2) = I.cols(1,m-1); P.row(0) = I.row(m-1);
    arma::mat Q(m,m); Q = -I; Q.col(0).fill(1); 
    
    arma::vec phi_i = L*exp( L*lin_can_par  );  
    double dem; dem = 1 + phi_i(m-1);  
    arma::vec num = Q*P*phi_i;
    
    dot_b = num/dem;
    
    return dot_b;
  }
  
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
  return -k_i*log(  (1 - mu_i(0))  ) ;
};   
          
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mn_loss_j(arma::vec c, 
                    arma::mat vj, 
                    arma::mat y_datta, 
                    arma::vec wj, 
                    Rcpp::String link, 
                    arma::vec k) {
  
  arma::uword pm=c.n_elem; 
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
      
      // Without constant gradient per class
      arma::mat tVij_I=kron( (vj.col(i)).t(),I);
      
      // Imposing constant gradient per class
      // matrix(Vj[,1], nrow=m-1, ncol=p+(m-1) ) 
      // arma::mat tVij_I=reshape(vj.col(i),m,pm );
      
      // Creating lin_can_parameter
      arma::vec lcp=tVij_I*c;
      
      mean_nll_j += -wj(i)*( lcp.t()*y_datta.col(i) - b_culmit(lcp, k(i) ) )/n;
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
  
  arma::uword pm=c.n_elem;
  arma::uword n=y_datta.n_cols;
  arma::uword m=y_datta.n_rows;
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
  }

  return mean_score_j;
  
};
  
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(name = "aD_j_cg_test")]]
arma::vec aD_j_cg_test(arma::vec init,
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
  double c_ag=control_list["c_ag"]; //parameter for Armijo condition, delta in Hybrid-DY
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
      double armijo_bound = as_scalar(c_ag*pow(beta_bt,m_cg)*s_now*
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
        
        wolfe_cond =+ (curv_wolfe >= wolfe_bound);
        
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
      Rcout << "Printing nll_dist: " << nll_dist<< ", m_ag:"  << m_ag << ", iter:"  << iter << " \ ";
    }
    
    if( nll_dist < tol) {
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
    
  } // end of cg iter  
  
  return c_now;//things;
  
};
')




#### cg quadratic test----
sourceCpp(code='
  #include <RcppArmadillo.h>
  using namespace Rcpp;
  using namespace arma;
  
  // [[Rcpp::depends(RcppArmadillo)]]
  // [[Rcpp::export(name="sq_loss")]]
  arma::vec sq_loss(arma::vec init,
                  arma::mat vj,
                  arma::mat y_datta
  ) {

    arma::uword m; m= y_datta.n_rows;
    arma::uword n; n= y_datta.n_cols;
    arma::uword p; p= vj.n_rows;
    arma::mat I(m,m); I.eye(); 
    
    arma::vec c_now = init;
    arma::vec y_vec = reshape(y_datta, n*m,1); 
    
    arma::mat nll_now(1,1);  
    nll_now = (y_vec-kron(vj.t(),I)*c_now).t()*(y_vec-kron(vj.t(),I)*c_now)/n;
    return nll_now;
  };
  
  // [[Rcpp::depends(RcppArmadillo)]]
  // [[Rcpp::export(name = "sq_grad")]]
  arma::mat sq_grad(arma::vec init,
                  arma::mat vj,
                  arma::mat y_datta
  ) {
  
    arma::uword m; m= y_datta.n_rows;
    arma::uword n; n= y_datta.n_cols;
    arma::uword p; p=vj.n_rows;
    arma::mat I(m,m); I.eye(); 
    
    // arma::mat c_mat; c_mat=reshape(init, p, m); 
    arma::vec c_now = init;
    arma::vec y_vec = reshape(y_datta, n*m,1); 
    
    arma::vec grad;
    grad = -2*kron( vj.t(),I).t()*(y_vec-kron( vj.t(),I)*c_now )/n;//
    return grad;//c_vec;//c_now; 
  };

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(name = "aD_j_cg_test")]]
arma::vec aD_j_cg_test(arma::vec init,
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
  double c_ag=control_list["c_ag"]; //parameter for Armijo condition, delta in Hybrid-DY
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
  nll_now = sq_loss(c_now, vj, y_datta); 
  grad_now = sq_grad(c_now,vj,y_datta);
  arma::vec p_now = grad_now; 
  
  arma::vec c_next;
  arma::uword iter;
  arma::vec things(7);

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
      double armijo_bound = as_scalar(c_ag*pow(beta_bt,m_cg)*s_now*
                                      (p_now.t()*grad_now));
      
      things(2)=armijo_bound; 
      
      // evaluation for armijo condition
      arma::vec c_search; c_search = c_now - pow(beta_bt,m_cg)*s_now*p_now;
      
      double suff_dec_ag;
      suff_dec_ag = as_scalar( sq_loss(c_now,vj,y_datta) - 
        sq_loss(c_search,vj,y_datta) );
      int armijo_cond = (suff_dec_ag >= armijo_bound);
      
      things(3)=suff_dec_ag;
      
      int armijo_cond2=0;
      if (c_ag2 > 0) {
        // the second bound in armijo-goldstein in nesterovs intro to conv opt text
        double armijo_bound2 = as_scalar(c_ag2*pow(beta_bt,m_cg)*s_now*
                                          (p_now.t()*grad_now));
        // second sufficient descent bound uses the same suff_dec_ag   
        armijo_cond2 =+ (suff_dec_ag <= armijo_bound2);
        
        things(4)=armijo_bound2;
      }  
      
      int wolfe_cond=0;
      if (c_wolfe > 0) {
        // the weak Wolfe condition
        double wolfe_bound = as_scalar(c_wolfe*(p_now.t()*grad_now));
        
        // evaluation for curvature in weak wolfe
        double curv_wolfe;
        curv_wolfe = as_scalar( p_now.t()*sq_grad(c_search,vj,y_datta) );
        
        wolfe_cond =+ (curv_wolfe >= wolfe_bound);
        
        things(5)=wolfe_bound;things(6)= curv_wolfe;
      }  
      
       things(1)=as_scalar(m_cg);
      
      
      if ( armijo_cond + armijo_cond2 == 2 ) { //== 2 + wolfe_cond 
        m_ag = m_cg;
        break;
      }
    }  
    things(0)=as_scalar(m_ag);
    
    double h_now = as_scalar( pow(beta_bt,m_ag)*s_now );
    c_next = c_now - h_now*p_now;
    
    // #Step 2a: Compute Loss;
    arma::mat nll_next(1,1); 
    nll_next=sq_loss(c_next,vj,y_datta); 
    double nll_dist; nll_dist = as_scalar( nll_now - nll_next);
    
    if (test) {
        // Rprintf("Printing: iter %iter, ll Dist %ll_dist, eu Dist %eu_dist ", 
        //         iter, ll_dist, eu_dist);
        Rcout << "Printing nll_dist, iter: " << nll_dist<< ", "  << iter << "\";
    }
    
    if( nll_dist < tol) {
      break;
    } else {
      
      // #Step 2b: Compute gradient;
      arma::vec grad_next = sq_grad(c_next,vj,y_datta);
      
      // Step 3: Compute the coeffiecient
      // Fletcher-Reeves
      double beta_cg_fr = -as_scalar( ( grad_next.t()*grad_next )/
                                      ( grad_now.t()*grad_now ) );  
      
      // Dai-Yuan
      double beta_cg_dy = -as_scalar( ( grad_next.t()*grad_next )/
                                     ( p_now.t()*(grad_next-grad_now) ) );
      
      // Hestenes-Stiefel
      double beta_cg_hs= -as_scalar( ( grad_next.t()*(grad_next-grad_now) )/
                                    ( p_now.t()*(grad_next-grad_now) ) );

      // Hybird
      arma::vec beta2(2); beta2(0) = beta_cg_dy; beta2(1)=beta_cg_hs;
      arma::vec beta3(2); beta3(0) = 0; beta3(1)=min(beta2);
      double beta_cg_hybrid=max(beta3);
      
      // Step 4: Update p
      arma::vec p_next = grad_next - beta_cg_fr*p_now;
      
      // Update all inputs for next iteration
      c_now=c_next;
      nll_now=nll_next;
      grad_now=grad_next;
      p_now=p_next;
    }
  
  } // end of cg iter  
 
  return c_now;
  
}')



#### dot_b ----
sourceCpp(code='// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export(name = "dot_b_multinom")]]
arma::vec dot_b_multinom(arma::vec lin_can_par, int k_i, String link ){
  arma::uword m=lin_can_par.n_rows;
  arma::vec dot_b;
  if ( link == "expit") {
    arma::vec e_lcp=exp(lin_can_par);
    // arma::vec ones_vec(m); ones_vec.ones();
    double dem; dem = 1 + sum(e_lcp) ;
    arma::vec pi_i=e_lcp/dem;
    arma::uvec ids = find(pi_i < 0);
    pi_i.elem(ids).fill(0);
    dot_b = k_i*pi_i;
    return dot_b;
  } else if (link == "culmit") {
    // creating upper/lower triangular matrices;
    arma::mat A; A.ones(m,m);
    //arma::mat U=trimatu(A);
    arma::mat L=trimatl(A);
    // creating Permutation and Differencing Matrices, P, Q
    arma::mat I(m,m); I.eye();
    arma::mat P(m,m); P.zeros(); P.cols(0,m-2) = I.cols(1,m-1); P.row(0) = I.row(m-1);
    arma::mat Q(m,m); Q = -I; Q.col(0).fill(1);
    arma::vec phi_i = L*exp( L*lin_can_par  );
    double dem; dem = 1 + phi_i(m-1);
    arma::vec num = Q*P*phi_i;
    dot_b = num/dem;
    return dot_b;
  }
};
// [[Rcpp::export(name = "b_culmit")]]
double b_culmit(arma::vec lin_can_par, int k_i) {
  arma::uword m=lin_can_par.n_rows;
  arma::vec mu_i=dot_b_multinom( lin_can_par, k_i, "culmit");
  return -k_i*log(  (1 - mu_i(0))  ) ;
};
// [[Rcpp::export(name = "mn_loss_j")]]
arma::mat mn_loss_j(arma::vec c, 
            arma::mat vj, 
            arma::mat y_datta, 
            arma::vec wj, 
            Rcpp::String link, 
            arma::vec k) {

  arma::uword pm=c.n_elem;
  arma::uword n=y_datta.n_cols;
  arma::uword m=y_datta.n_rows;
  arma::mat I(m,m); I.eye();
  
  arma::mat mean_nll_j(1,1); mean_nll_j.zeros();   
  // arma::mat test;
  
  if (link=="culmit") {
    
    // Writing the For loop instead of sapply.
    arma::uword i;
    for (i = 0; i < n; i++ ) {
      
      // Without constant gradient per class
      // arma::mat tVij_I=kron( (vj.col(i)).t(),I);
      
      // Imposing constant gradient per class
      // matrix(Vj[,1], nrow=m-1, ncol=p+(m-1) ) 
      arma::mat tVij_I=reshape(vj.col(i),m,pm );
      
      // Creating lin_can_parameter
      arma::vec lcp=tVij_I*c;
      
      mean_nll_j += -wj(i)*( lcp.t()*y_datta.col(i) - b_culmit(lcp, k(i) ) )/n;
    } 
  }
  
  return mean_nll_j;
      
};
//[[Rcpp::export(name = "mn_score_j")]]
arma::mat mn_score_j(arma::vec c,
                     arma::mat vj,
                     arma::mat y_datta,
                     arma::vec wj,
                     Rcpp::String link,
                     arma::vec k) { 
  arma::uword n=y_datta.n_cols;
  arma::uword m=y_datta.n_rows;
  arma::uword pm=c.n_elem;
  arma::mat I(m,m); I.eye();
  
  arma::mat mean_score_j(pm,1); mean_score_j.zeros();
  // arma::mat test;
  
  // Link only matters for the mean, the form of the score is
  // always the same;
  // Writing the For loop instead of sapply.
  
  
  if (link=="culmit"){
    
    arma::uword i;
    for (i = 0; i < n; i++ ) {
      
      // Without constant gradient per class
      // arma::mat tVij_I=kron( (vj.col(i)).t(),I);
      
      // Imposing constant gradient per class
      // matrix(Vj[,1], nrow=m-1, ncol=p+(m-1) ) 
      arma::mat tVij_I=reshape(vj.col(i),m,pm );
      
      // Creating lin_can_parameter
      arma::vec lcp=tVij_I*c;
      
      arma::vec mu_ij = dot_b_multinom(lcp, k(i), link);
      
      mean_score_j += -wj(i)*tVij_I.t()*( y_datta.col(i) - mu_ij)/n; 
    }
    
  } else {
    
    arma::uword i;
    for (i = 0; i < n; i++ ) {
      
      arma::mat tVij_I=kron( (vj.col(i)).t(),I);
      arma::vec lcp=tVij_I*c;
      arma::vec mu_ij = dot_b_multinom(lcp, k(i), link);
      
      mean_score_j += -wj(i)*tVij_I.t()*( y_datta.col(i) - mu_ij)/n; 
    }
    
    
  }
  
  return mean_score_j;
  
};
// [[Rcpp::export(name = "aD_j_cg")]]
arma::vec aD_j_cg(arma::vec init,
                  arma::mat vj,
                  arma::mat y_datta,
                  arma::vec wj,
                  Rcpp::String link,
                  arma::vec k, 
                  Rcpp::List control_list,
                  bool test
) {
  // Set out controls for the algorithm
  
  double tol = control_list["tol_val"];
  arma::uword max_iter = control_list["max_iter"];
  
  // beta and sigma in backtracking rule (Bertsekas)
  // inital step size s_k that will be adjusted by backtracking 
  arma::vec s=control_list["init_stepsize"]; 
  double beta_bt=control_list["beta_bt"]; //backtracking beta
  double c_ag=control_list["c_ag"]; //parameter for Armijo condition, delta in Hybrid-DY
  double c_ag2=control_list["c_ag2"]; //parameter for second Armijo bound in nesterov 
  double c_wolfe=control_list["c_wolfe"]; //parameter for curvature condition, sigma in H-DY
  int max_iter_line=control_list["max_iter_line"];
  //int l2_pen=control_list["l2_pen"]; //penalty parameter for ridge regression
  
  // Step 0: Set the initial value, and compute the loss and score
  // Also set the initial p_0 to the gradient
  
  arma::vec c_now = init; 
  
  arma::mat nll_now(1,1);  
  arma::vec grad_now; 
  nll_now = mn_loss_j(c_now,vj,y_datta,wj,link,k);
  grad_now = mn_score_j(c_now,vj,y_datta,wj,link,k);
  arma::vec p_now = grad_now; 
  
  // //   grad_now=score_fn(c_now);
  // //   p_now = grad_now; 
  
  arma::vec c_next;
  arma::uword iter;
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
      double armijo_bound = as_scalar(c_ag*pow(beta_bt,m_cg)*s_now*
                                        (p_now.t()*grad_now));
      
      
      
      // evaluation for armijo condition
      arma::vec c_search; c_search = c_now - pow(beta_bt,m_cg)*s_now*p_now;
      
      double suff_dec_ag;
      suff_dec_ag = as_scalar( mn_loss_j(c_now,vj,y_datta,wj,link,k) - 
                                 mn_loss_j(c_search,vj,y_datta,wj,link,k) );
      int armijo_cond = (suff_dec_ag >= armijo_bound);
      
      
      int armijo_cond2=0;
      
      double armijo_bound2 = as_scalar(c_ag2*pow(beta_bt,m_cg)*s_now*
                                           (p_now.t()*grad_now));
      // second sufficient descent bound uses the same suff_dec_ag   
      armijo_cond2 =+ (suff_dec_ag <= armijo_bound2);
      
      //if (c_ag2 > 0) {
        // the second bound in armijo-goldstein in nesterovs intro to conv opt text
      //}  
      
      int wolfe_cond=0;
      if (c_wolfe > 0) {
        // the weak Wolfe condition
        double wolfe_bound = as_scalar(c_wolfe*(p_now.t()*grad_now));
        
        // evaluation for curvature in weak wolfe
        double curv_wolfe;
        curv_wolfe = as_scalar( p_now.t()*mn_score_j(c_search,vj,y_datta,wj,link,k) );
        
        wolfe_cond =+ (curv_wolfe >= wolfe_bound);
      }  
      
      if ( armijo_cond == 1 ) { //+ wolfe_cond
        m_ag = m_cg;
        break;
      }
      
    }
    
    double h_now = as_scalar(pow(beta_bt,m_ag)*s_now);
    c_next = c_now - h_now*p_now;
    
    
    // #Step 2a: Compute Loss;
      arma::mat nll_next(1,1);
      nll_next=mn_loss_j(c_next,vj,y_datta,wj,link,k); 
    
    
    double nll_dist; nll_dist = as_scalar( nll_now - nll_next);
    
    if (test) { 
      Rcout << "Printing nll_dist, iter: " << nll_dist<< ", "  << iter << "\" ;
    };
    
    if( nll_dist < tol) {
      break;
    } else {
      
      // #Step 2b: Compute gradient;
        arma::vec grad_next = mn_score_j(c_next,vj,y_datta,wj,link,k);
        
        // Step 3: Compute the coeffiecient
        // Fletcher-Reeves
        // double beta_cg_fr = as_scalar( ( grad_next.t()*grad_next )/
        //                                ( grad_now.t()*grad_now ) );
        
        // Dai-Yuan
        double beta_cg_dy = as_scalar( -( grad_next.t()*grad_next )/
                                         ( p_now.t()*(grad_next-grad_now) ) );
        
        // Hestenes-Stiefel
        double beta_cg_hs= as_scalar( -( grad_next.t()*(grad_next-grad_now) )/
                                        ( p_now.t()*(grad_next-grad_now) ) );
        
        // Hybird
        arma::vec beta2(2); beta2(0) = beta_cg_dy; beta2(1)=beta_cg_hs;
        arma::vec beta3(2); beta3(0) = 0; beta3(1)=min(beta2);
        double beta_cg_hybrid=max(beta3);
        
        // Step 4: Update p
        arma::vec p_next = grad_next - beta_cg_hybrid*p_now;
        
        // Update all inputs for next iteration
        c_now=c_next;
        nll_now=nll_next;
        grad_now=grad_next;
        p_now=p_next;
    }
  }
  
  return c_next;
  
};')



