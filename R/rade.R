###############################################################################
################ RADE and Regularized OPCG - Ridge and LASSO ##################
###############################################################################

#' Regularized Regressing Average Derivative Estimates
#'
#' @param x_matrix a 'pxn' matrix of predictors;
#' @param y_matrix a 'mxn' matrix response 
#' @param d specified the reduced dimension 
#' @param bw the bandwidth parameter for the kernel; the default kernel is gaussian
#' @param ytype specify the response as 'continuous', 'multinomial', or 'ordinal' 
#' @param method Default is 'newton'. Select the method for optimization. Strongly 
#' suggest the conjugate gradient method 'cg' when 'p' of moderate size. 
#' @param parallelize Default is False; to run in parallel, you will need to have
#' foreach and some parallel backend loaded; parallelization is strongly recommended
#' and encouraged.
#' @param l2_pen Default is '0', i.e. no ridge penalty. 
#' @param l1_pen Default is '0', i.e no LASSO Penalty. The LASSO estimation is carried out
#' through the 'glmnet' package. To select the L1 penalty through Cross-Validation in 
#' glmnet, set 'l1_pen = -1'. Otherwise, provide a sequence of penalty values (glmnet 
#' strongly discourages supplying a single value.)
#' @param control_list a list of control parameters for the Newton-Raphson 
#' or Conjugate Gradient methods
#' \itemize{
#'   \item opcg - A 'pxd' matrix that estimates a basis for the central subspace.
#'   \item opcg_wls - A 'pxd' matrix that estimates a basis for the central subspace based 
#'   on the initial value of the optimization problem; useful for examining bad starting 
#'   values.
#'   \item cand_mat - A list that contains both the candidate matrix for OPCG and for
#'   the initial value; this is used in other functions for order determination
#'   \item gradients - The estimated local gradients; used in regularization of OPCG
#'   \item weights - The kernel weights in the local-linear GLM. 
#' }
#' @param bw2 Default is NULL; NULL refers to not using refined weights in the MADE step
#' If bw2 is a scalar, then refined weights is used with bandwidth bw2
#' @param init_mat Default is NULL; refers to the initial matrix in the MADE step. NULL
#' means that the OPCG estimate is used as the initial value; 
#'
#' @return
#' 
#' @export
#'
#' @examples
#' 
rade<-function(x_matrix, y_matrix, d, bw, bw2=NULL, init_mat=NULL, ytype='continuous', 
               method="newton", parallelize=F, r_mat=NULL,
               l2_pen=0, l1_pen=0, print_opcg=F, control_list=list() ){
  
  # x_matrix=X; y_matrix=Y; d; bw; bw2=bw; ytype="multinomial";#"continuous";
  # tol_val= 1e-07; max_iter=25;
  # init_mat = NULL;  method="cg"; parallelize=T; r_mat=NULL;
  # control_list=list(); l2_pen=0; l1_pen=0;
  # ytype="ordinal"
  
  # Control Parameter Defaults
  control_args=control_list; control_names=names(control_args);
  tol_val=if ( "tol_val" %in% control_names ) control_args$tol_val else 1e-7; 
  max_iter=if ( "max_iter" %in% control_names ) control_args$max_iter else 25 ; 
  print_iter=if ( "print_iter" %in% control_names ) control_args$print_iter else F ; 
  
  # Setting parameters
  p = dim(x_matrix)[1]; n = dim(x_matrix)[2];  
  eye_d = diag(1,d,d);
  
  # l2_pen=1; l1_pen=0;
  # l2_pen=0; l1_pen=1;
  
  opcg_obj=opcg_wrap(x_matrix, y_matrix, d, bw, ytype,
                     method, parallelize, r_mat,
                     control_list)
  
  b_hat_opcg=opcg_obj$opcg;
  opcg_grad=opcg_obj$gradients;
  
  if (is.null(init_mat) ) init_mat=b_hat_opcg; # Not needed
  
  if (is.null(bw2)) { #diag(1, p, d)
    Bx_matrix0=t(b_hat_opcg)%*%x_matrix;
    Bx_matrix=matpower_cpp(cov(t(Bx_matrix0)), -1/2)%*%
      t(sapply(1:d, FUN= function(j) center_cpp(Bx_matrix0[j,], NULL) ) );
    
    made_grad=opcg_made(Bx_matrix, y_matrix, bw, B_mat=NULL, ytype, method, 
                        parallelize, r_mat=NULL, control_list)$Dhat;
  } else {
    Bx_matrix0=t(b_hat_opcg)%*%x_matrix;
    Bx_matrix=matpower_cpp(cov(t(Bx_matrix0)), -1/2)%*%
      t(sapply(1:d, FUN= function(j) center_cpp(Bx_matrix0[j,], NULL) ) );
    
    made_grad=opcg_made(Bx_matrix, y_matrix, bw2, B_mat=NULL, ytype, method, 
                        parallelize, r_mat=NULL, control_list)$Dhat;
  }
  
  # if (is.null(init_mat) ) init_mat=b_hat_opcg;
  # 
  # if (is.null(bw2)) { #diag(1, p, d)
  #   made_grad=opcg_made(x_matrix, y_matrix, bw, B_mat=init_mat, ytype, method, 
  #                       parallelize, r_mat=NULL, control_list)$Dhat;
  # } else {
  #   made_grad=opcg_made(x_matrix, y_matrix, bw2, B_mat=init_mat, ytype, method, 
  #                       parallelize, r_mat=b_hat_opcg, control_list)$Dhat;
  # }
  
  
  # made_grad_kp = lapply(1:n, function(j) kronecker( t(made_grad[[j]]), diag(1,p,p) ));
  # made_grad_tkp = lapply(1:n, function(j) t(made_grad_kp[[j]]) );
  # opcg_grad_tvec=lapply(1:n, function(j) t(matrix(c(opcg_grad[[j]]) ) ) ) ;
  
  # gamma_hat= solve_cpp(list_sum(made_grad_tkp,made_grad_kp) + l2_pen*diag(1, p*d, p*d),
  #                      list_sum(made_grad_kp, opcg_grad_tvec) )
  # 
  # Gamma_hat= matrix(gamma_hat, p, d)
  
  made_grad_t = lapply(1:n, function(j) t(made_grad[[j]]) );
  
  
  # LASSO update for Gamma
  
  if(l1_pen[1] != 0){
    # i=1
    vvt0=list_sum(made_grad,made_grad);
    Gamma_hat=t(solve_cpp(vvt0, list_sum(made_grad, opcg_grad) ) )
    
    vvt0_rt= matpower_cpp(vvt0, .5);
    x_star = lapply(1:p, function(i) rbind( vvt0_rt, sqrt(l2_pen)*eye_d ))
    y_star = lapply(1:p, function(i) c( vvt0_rt%*%Gamma_hat[i,], rep(0,d)) ) 
    
    if(l1_pen[1] < 0){
      Gamma_hat_l0=foreach(i=1:p, .combine='rbind', .packages =  "glmnet" ) %dopar% {
        # i=1
        # pen_obj=glmnet(x_star[[i]], y_star[[i]], intercept=F, family="gaussian", alpha=1)
        # print(pen_obj)
        # pen_obj$lambda
        cv_obj=glmnet::cv.glmnet(x_star[[i]], y_star[[i]], intercept=F, 
                                 family="gaussian", alpha=1, standardize = F)
        # coef(pen_obj, s=cv_obj$lambda.min)
        # pen_obj$beta
        gamma_hat_l0_i=as.vector(glmnet::glmnet(x_star[[i]], y_star[[i]], intercept=F, 
                                                family="gaussian", 
                                        alpha=1, lambda=cv_obj$lambda.min,
                                        standardize = F)$beta)
        return(gamma_hat_l0_i);
      }
      Gamma_hat = matrix(Gamma_hat_l0, p, d);
    } else if (l1_pen[1] > 0) {
      Gamma_hat_l0=foreach(i=1:p, .combine = 'rbind', .packages =  "glmnet" ) %dopar% {
        # i=1
        # gamma_hat_l0_i=as.vector(glmnet::glmnet(x_star[[i]], y_star[[i]], intercept=F,
        #                                 family="gaussian", alpha=1, lambda = l1_pen,
        #                                 standardize = F)$beta)
        cv_obj=glmnet::cv.glmnet(x_star[[i]], y_star[[i]], intercept=F, 
                                 family="gaussian", alpha=1, lambda = l1_pen,
                                 standardize = F)
        # coef(pen_obj, s=cv_obj$lambda.min)
        # pen_obj$beta
        gamma_hat_l0_i=as.vector(glmnet::glmnet(x_star[[i]], y_star[[i]], intercept=F, 
                                                family="gaussian", 
                                                alpha=1, lambda=cv_obj$lambda.min,
                                                standardize = F)$beta)
        
        return(gamma_hat_l0_i);
      }
      Gamma_hat = Gamma_hat_l0;
    }
    
  } else if (l1_pen==0 & l2_pen < 0) { # Just L2 penalty, using glmnet
    # made_grad_kp = lapply(1:n, function(j) kronecker( t(made_grad[[j]]), diag(1,p,p) ));
    x_star = lapply(1:p, function(i) kronecker( t(made_grad[[i]]), diag(1,p,p) ) )
    y_star = lapply(1:p, function(i) c( opcg_grad[[i]] ) )  
    
    Gamma_hat_l1=foreach(i=1:p, .combine='rbind', .packages =  "glmnet" ) %dopar% {
      # i=1
      # pen_obj=glmnet(x_star[[i]], y_star[[i]], intercept=F, family="gaussian", alpha=1)
      # print(pen_obj)
      # pen_obj$lambda
      cv_obj=glmnet::cv.glmnet(x_star[[i]], y_star[[i]], intercept=F, 
                               family="gaussian", alpha=0, standardize = F)
      # coef(pen_obj, s=cv_obj$lambda.min)
      # pen_obj$beta
      gamma_hat_l0_i=as.vector(glmnet::glmnet(x_star[[i]], y_star[[i]], intercept=F, 
                                              family="gaussian", 
                                              alpha=0, lambda=cv_obj$lambda.min,
                                              standardize = F)$beta)
      return(gamma_hat_l0_i);
    }
    Gamma_hat = matrix(Gamma_hat_l1, p, d);
    
  } else if (l1_pen==0 & l2_pen >= 0) { # Just L2 penalty, if even
    vvt0=list_sum(made_grad,made_grad) + l2_pen*eye_d;
    Gamma_hat=t(solve_cpp(vvt0, list_sum(made_grad, opcg_grad) ) )
  }
  
  
  Gamma_hat = apply(Gamma_hat, 2, normalize_cpp) 
  if (print_opcg) {
    return(list(opcg=b_hat_opcg, rade=Gamma_hat) );
  } else if (!print_opcg) {
    return(Gamma_hat);
  }
  
} 


# A wrapper for rRADE

rrade=function(x_matrix, y_matrix, d, bw, bw2=NULL, init_mat=NULL, ytype='continuous', 
               method="newton", parallelize=F, 
               l2_pen=0, l1_pen=0, control_list=list() ){
  
  est0=rade(x_matrix, y_matrix, d, bw, bw2, init_mat, ytype, 
            method, parallelize, r_mat = NULL,
            l2_pen, l1_pen, control_list )
  
  refined_est0 = rade(x_matrix, y_matrix, d, bw2, bw2, init_mat, ytype, 
                  method, parallelize, r_mat=est0,
                  l2_pen, l1_pen, control_list )
  
  for(iter in 1:5) {
    refined_est1 = rade(x_matrix, y_matrix, d, bw2, bw2, init_mat, ytype, 
                        method, parallelize, r_mat=refined_est0,
                        l2_pen, l1_pen, control_list )
    
    dist = mat_dist(refined_est1, refined_est0);
    
    print(c("rOPCG: dist dist is", dist, iter));
    if( dist < 1e-7 ) {
      refined_est0=refined_est1;
      break();
    } else{
      # The new B_0 for next iteration
      # B_0 = B_1;
      refined_est0=refined_est1;
    }  
    if(iter==25) print("0 - non-convergence");
    
  } # end of iterations
  
  return(refined_est0)
  
}

