#######################################################################
#           A generic Wrapper for OPCG and MADE
#######################################################################

# This will make use of doParallel and foreach, so we will need to call
# or import these objects/functions for our use
# Coatless has an example here:
# devtools::install_github("r-pkg-examples/rcpp-and-doparallel")
# library("Rcpp2doParallel")
# that we will replicate

# source("opcg_wrap_cpp.R")

############################## MADE wrappers ##################################


# x_matrix=X; y_matrix=Y; bw=1.25; ytype="multinomial"; B_mat = diag(1,p,d);
# method=list(opcg="newton", made="newton"); parallelize = T; r_mat=NULL; control_list=list(c_ag2=.9);
# aD_list=opcg_made(x_matrix, y_matrix, bw, B_mat, ytype,
#                 method=method$opcg, parallelize, r_mat=NULL,
#                 control_list);
# ahat_list = aD_list$ahat;  Dhat_list = aD_list$Dhat;

#### MADE Block update ----
#                 method, parallelize, r_mat=NULL,
#                 control_list);
#' 
#'
#'
#' @keywords internal
#' @noRd
#' 
made_update = function(x_matrix, y_matrix, d, bw, aD_list ,B_mat,  ytype="continuous",
                       method, parallelize=F, r_mat=NULL,
                       control_list=list()) {
  
  # Control Parameter Defaults
  control_args=control_list; control_names=names(control_args); 
  tol_val=if ( "tol_val" %in% control_names ) control_args$tol_val else 1e-7; 
  max_iter=if ( "max_iter" %in% control_names ) control_args$max_iter else 25 ;
  max_iter_made=if ( "max_iter_made" %in% control_names ) {
    control_args$max_iter_made
  } else {
    25
  };
  print_iter_made=if ( "print_iter_made" %in% control_names ) control_args$print_iter_made else F;
  
  # Setting parameters
  p=dim(x_matrix)[1]; n=dim(x_matrix)[2];  
  ahat_list = aD_list$ahat;  Dhat_list = aD_list$Dhat;
  if(is.null(r_mat)) r_mat = diag(1,p,p);

  # Initial c_param value
  c_init = as.vector(t(B_mat));
  
  # Setting up the response ----
  if (ytype=='continuous'){
    # y_matrix should be m x n, m can be >= 1
    # This is just OPG, with the WLS as the exact solution per j
    mv_Y=y_matrix; m=dim(y_matrix)[1];
    
    # # Loss function
    # loss_made = function(c_param){ 
    #   mn_loss_made(c_param, x_matrix, mv_Y, bw, ahat_list, Dhat_list,
    #                link=linktype, k=k_vec, r_mat)
    # } # End of Loss function
    # loss_made(c_init)
    

  } else if (ytype %in% c( "multinomial", "ordinal")) {
    # y.matrix should be 1 x n

    m_classes=as.numeric(levels(as.factor(y_matrix)));
    m=length(m_classes);
    mv_Y = mnY_to_mvY( y_matrix, m_classes, ytype);

    # Don't need Empirical Link Transforms for MADE block step
    if (ytype=="multinomial" ) {
      linktype="expit";
    } else if (ytype=="ordinal" ) {
      linktype="culmit";
    };
    k_vec = rep(1, n);  
    mv_Y=matrix(mv_Y[1:(m-1),], m-1, n)
    
    # Loss function
    loss_made = function(c_param){ 
      mn_loss_made(c_param, x_matrix, mv_Y, bw, ahat_list, Dhat_list,
                   link=linktype, k=k_vec, r_mat)
    } # End of Loss function
    
    # loss_made(c_init)
    
    # # Score function
    # score_fn = function(c_param) {
    #   mn_score_made(c_param, x_matrix, mv_Y, bw, ahat_list, Dhat_list,
    #                 link=linktype, k=k_vec, r_mat)
    # } # End of score; score_fn(c_init)
    
  } # end of setting up response
  
  
  # Running Optimization algorithm
  
  if (method=="newton") {
    # Estimation using Newton-Raphson
    
    # This function returns the next newton iterate, i.e. one Newton Step
    c_next = function(c_param) { 
      
      # Compute the loss, score and info 
      # This will need all loss, score and info over j and i
      
      # For each j
      c_list_newton_j=function(c_param, j){ # j=10; c_param=c_init;
        
        Xj=linearsdr:::matcenter_cpp(x_matrix, index=j,x0=NULL);
        Wj=linearsdr:::gauss_kern_cpp( t(r_mat)%*%Xj,bw)
        
        # Loss functions and related functions
        
        score_j=linearsdr:::mn_score_j_made(c_param, Xj, mv_Y, Wj,
                                            ahat=ahat_list[[j]], Dhat = Dhat_list[[j]],
                                            link=linktype, k=k_vec);
        info_j=linearsdr:::mn_info_j_made(c_param, Xj, mv_Y, Wj,
                                          ahat=(ahat_list[[j]]), Dhat = Dhat_list[[j]],
                                          link=linktype, k=k_vec) ;
        
        return(list(score=score_j, info=info_j))
      } # c_list_newton_j(c_init,231)
      
      # Computing all local loss, score, info
      # Note the loss function is evaluated at c_0, not c_1
      if (parallelize) {
        # Should return three lists, each length n; list of loss, list of score, 
        # list of info
        c_list_newton = foreach::foreach(j = iterators::icount(n),
                                  .packages = "linearsdr",
                                  .export=c("x_matrix",
                                            "mv_Y",
                                            "r_mat", "bw",
                                            "ahat_list" ,"Dhat_list",
                                            "linktype", "k_vec")
                                  ) %dopar% {
                                    # Release results
                                    return(c_list_newton_j(c_param,j))
                                  }
      } else { # No parallel, so for loop 
        c_list_newton = list();
        for(j in 1:n) {
          c_list_newton[[j]] = c_list_newton_j(c_param, j);
        }; 
      } # end of computing loss, score, info
    
      score=lapply( 1:n , FUN=function(j) c_list_newton[[j]][[1]] );
      info=lapply( 1:n , FUN=function(j) c_list_newton[[j]][[2]] );
      # c_1 = vecB_hat(c_param, score, info);
      
      c_1=c_param - ginv(Reduce('+', info))%*%Reduce('+', score)

      return(c_1)
      
    } # End of function returning next newton iterate. 
    
    
    c_0 = c_next(c_init);
    loss_0 = loss_made(c_0);
    
    for(iter in 1:max_iter_made){ #print_B_iter=T
      
      B_0=t(matrix(c_0, nrow=d, ncol=p));
      
      c_1=c_next(c_0); loss_1 = loss_made(c_1);
      
      B_1=t(matrix(c_1, nrow=d, ncol=p));
      
      # A Matrix distance of B_1 to B_0;
      subspace_dist=mat_dist(B_1, B_0);
      
      # Vector Distance
      euc_dist=euc_norm_cpp(c_0 - c_1)/euc_norm_cpp(c_0);
      
      # Loss Distance 
      loss_dist = (loss_0 - loss_1)/loss_0;
      
      if(print_iter_made) print(c("B Update: loss_dist is", loss_dist, 
                               "euc_dist is", euc_dist,
                               "subspace_dist is", subspace_dist,
                               "Iter:", iter));
      if( loss_dist < tol_val ) {
        break();
      } else {
        # The new B_0 for next iteration
        # B_0 = B_1;
        c_0=c_1; loss_0=loss_1; 
      }
      
    }
    
  # End of Newton Algorithm  
  } else if (method=="cg") { 
    # Estimation using Conjugate Gradients
    
    # Control Parameter Defaults
    control_args=control_list; control_names=names(control_args); 
    test=if ( "test" %in% control_names ) control_args$test else F ; 
    max_iter_made=if ( "max_iter_made" %in% control_names ) control_args$max_iter_made else 25 ; 
    init_stepsize_made=if ( "init_stepsize_made" %in% control_names ) control_args$init_stepsize_made else rep(n,max_iter_made); 
    beta_bt_made=if ( "beta_bt_made" %in% control_names ) control_args$beta_bt_made else 0.5;
    c_ag1_made=if ( "c_ag1_made" %in% control_names ) control_args$c_ag1_made else 10e-3;
    c_ag2_made=if ( "c_ag2_made" %in% control_names ) control_args$c_ag2_made else 0.9;
    c_wolfe_made=if ( "c_wolfe_made" %in% control_names ) control_args$c_wolfe_made else 0; # 0.1 wiki-recom 
    max_iter_line_made=if ( "max_iter_line_made" %in% control_names ) control_args$max_iter_line_made else 100;
    
    # Initial c_param value
    c_init = as.vector(t(B_mat));
    
    c_0=vecB_cg(c_init, x_matrix, mv_Y, 
                   bw, ahat_list, Dhat_list,
                   link=linktype, k=k_vec, r_mat,
                   control_list=list(tol_val=tol_val,
                                     max_iter=max_iter_made, 
                                     init_stepsize=init_stepsize_made,
                                     beta_bt=beta_bt_made,
                                     c_ag1=c_ag1_made,
                                     c_ag2=c_ag2_made,
                                     c_wolfe=c_wolfe_made, 
                                     max_iter_line=max_iter_line_made),
                   test) 
  } # End of CG   
  return(c_0);

} # End of made update

# made_update(x_matrix, y_matrix, d, bw, aD_list ,B_mat,  ytype="multinomial",
#             method= method$made, parallelize=T, r_mat=NULL,
#             control_list=list( ))

#### MADE; An alternating Block Co-ordinate Descent Approach ----
# 

#' Minimum Average Deviance Estimation
#'
#' This implements the Outer Product of Canonical Gradients (OPCG) in a forth coming
#' paper Quach and Li (2021).
#' 
#' This version of MADE differs from that of Adragni and is currently available
#' for continuous, multinomial, or ordinal response so far.
#' 
#' For scalar continuous response, the estimation is identical to OPG.
#' 
#' The kernel for the local linear regression is fixed at a gaussian kernel.
#' 
#' For large 'p', we strongly recommend using the Conjugate Gradients implement, 
#' by setting method="cg".
#' For method="cg", the hybrid conjugate gradient of Dai and Yuan is implemented,
#' but only the armijo rule is implemented through backtracking, like in Bertsekas'
#' "Convex Optimization Algorithms".
#' A weak Wolfe condition can also be enforced by adding setting c_wolfe > 0 
#' in the control_list, but since c_wolfe is usually set to 0.1 (Wikipedia)
#' and this drastically slows down the algorithm relative to newton for small to 
#' moderate p, we leave the default as not enforcing a Wolfe condition, since we 
#' assume that our link function gives us a close enough initial point that local
#' convergence is satisfactory. Should the initial values be suspect, then maybe
#' enforcing the Wolfe condition is a reasonable trade-off.  
#' 
#' @param d specified the reduced dimension 
#' @param x_matrix a 'pxn' matrix of predictors;
#' @param y_matrix a 'mxn' matrix response 
#' @param bw the bandwidth parameter for the kernel; the default kernel is gaussian
#' @param method "newton" or "cg" methods; for carrying out the optimization using
#' the standard newton-raphson (i.e. Fisher Scoring) or using Congugate Gradients 
#' @param parallelize Default is False; to run in parallel, you will need to have
#' foreach and some parallel backend loaded; parallelization is strongly recommended
#' and encouraged.
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
#' @param ytype specify the response as 'continuous', 'multinomial', or 'ordinal' 
#' @param B_mat 
#' @param r_mat 
#'
#' @return A list containing both the estimate and candidate matrix.
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
#'  
#' @export
#' 
made <- function(x_matrix, y_matrix, d, bw, B_mat=NULL, ytype="continuous",
                 method=list(opcg="newton", made="newton"), parallelize=F, r_mat=NULL,
                 control_list=list()) {
  
  # x_matrix=X; y_matrix=Y; d; bw; B_mat=B_hat_opcg; ytype="multinomial";
  # x_matrix=X; y_matrix=Y; d; bw; B_mat=NULL; ytype="multinomial";
  # method=list(opcg="newton", made="newton"); parallelize=T; r_mat=NULL;
  # control_list=list(c_ag2=.9, print_iter=T, max_iter_made=25)
  
  # Control Parameter Defaults
  control_args=control_list; control_names=names(control_args);
  tol_val=if ( "tol_val" %in% control_names ) control_args$tol_val else 1e-7;
  max_iter_made=if ( "max_iter_made" %in% control_names ) {
    control_args$max_iter_made
  } else {
    25
  };
  print_iter=if ( "print_iter" %in% control_names ) control_args$print_iter else F;
  
  
  
  # Setting parameters
  p = dim(x_matrix)[1];
  if (is.null(B_mat)) B_mat= diag(1,p,d);
  if (is.null(r_mat)) r_mat= diag(1,p,p);
  
  # Block step for aD parameter
  aDhat=opcg_made(x_matrix, y_matrix, bw, B_mat, ytype,
                  method=method$opcg, parallelize, r_mat=NULL,
                  control_list);
  # Block step for B parameter
  c_0=made_update(x_matrix, y_matrix, d, bw,
                  aD_list = aDhat,
                  B_mat, ytype, method=method$made, parallelize, r_mat,
                  control_list);
  
  if (ytype=="continuous") {
    mv_Y=y_matrix;
  } else if (ytype=="multinomial" ) { 
    linktype="expit";
    mv_Y=linearsdr:::mnY_to_mvY( y_matrix, m_classes, ytype);
    
    k_vec = colSums(mv_Y);
    mv_Y=matrix(mv_Y[1:(m-1),], m-1, n)
    
    
  } else if (ytype=="ordinal" ) {
    linktype="culmit";
    mv_Y=linearsdr:::mnY_to_mvY( y_matrix, m_classes, ytype);
    
    k_vec = rep(1, n) #as.vector(y_matrix);
    mv_Y=matrix(mv_Y[2:(m),], m-1, n) # Drop the first row now cause its all 1

    
  }

  
  loss_0 = mn_loss_made(c_0, x_matrix, mv_Y, bw,
                        ahat_list=aDhat$ahat, Dhat_list=aDhat$Dhat,
                        link=linktype, k=k_vec, r_mat)
  
  for(iter in 1:max_iter_made){
    B_0=t(matrix(c_0, nrow=d, ncol=p));

    # aD-block update
    aDhat1=opcg_made(x_matrix, y_matrix, bw, B_mat=B_0, ytype,
                     method=method$opcg, parallelize, r_mat=NULL,
                    control_list);
    # B-block update
    
    c_1=made_update(x_matrix, y_matrix, d,bw,
                    aD_list = aDhat1, B_0,
                    ytype, method=method$made, parallelize, r_mat=NULL,
                    control_list);
    
    # Loss function
    loss_1 = mn_loss_made(c_1, x_matrix, mv_Y, bw,
                          ahat_list=aDhat1$ahat, Dhat_list=aDhat1$Dhat,
                          link=linktype, k=k_vec, r_mat)
    
    # A Matrix distance of B_1 to B_0;
    # B_1=t(matrix(c_1, nrow=d, ncol=p));
    # subspace_dist=mat_dist(B_1, B_0);
    
    euc_dist=euc_norm_cpp(c_0 - c_1)/euc_norm_cpp(c_0);
    
    # Loss Distance 
    loss_dist = (loss_0 - loss_1)/loss_0;
    # subspace_dist; euc_dist;
    if(print_iter) print(c("MADE: loss_dist dist is", euc_dist, iter));
    if( loss_dist < tol_val ) {
      break();
    } else{
      # The new B_0 for next iteration
      # B_0 = B_1;
      c_0=c_1;
    }  
    if(iter==max_iter_made) print("0 - non-convergence")

  }

  B_hat_made = t(matrix(c_1, nrow=d, ncol=p));
  B_hat_made = apply(B_hat_made, 2, normalize_cpp);
  return(B_hat_made)
  # return(list( estimate=B_hat_made, loss=B_hat_made_loss))
};



