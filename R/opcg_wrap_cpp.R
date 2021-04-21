#######################################################################
#           A generic Wrapper for OPCG and MADE
#######################################################################

# This will make use of doParallel and foreach, so we will need to call
# or import these objects/functions for our use
# Coatless has an example here: 
# devtools::install_github("r-pkg-examples/rcpp-and-doparallel")
# library("Rcpp2doParallel")
# that we will replicate

# OPCG-MADE Wrapper ####
#' This is an internal function called by opcg 
#' 
#'
#'
#' @keywords internal
#' @noRd
#' 
opcg_made <- function(x_matrix, y_matrix, bw, B_mat=NULL, ytype='continuous', 
                      method="newton", parallelize=F, r_mat=NULL, 
                      control_list=list()) {
  # y.matrix should be m x n, with m depending on ytype
  
  # Supported ytypes are 
  # continuous: Defaults to Squared Loss
  # multinomial: Multinomial-Categorical GLM with categorical link
  # ordinal: Multinomial-Ordinal GLM with ordinal link  
  # other: Custom Loss Functions - to be done at a later date
  
   
  # x_matrix=X; y_matrix=Y;#matrix(Y[2,],1,n)#Y;
  # bw;  ytype="clogit"; #"continuous";#"multinomial";
  # tol_val= 1e-07; max_iter=25;
  # B_mat = NULL ; method="cg"; parallelize=T; r_mat=NULL; control_list=list();
  # control_list = list() #control_list = list(); # B_mat=init_mat;
  
  
  # Parameters for the problem/model
  p <- dim(x_matrix)[1]; n <- dim(x_matrix)[2]; 
  B <- if(is.null(B_mat)) diag(1,p,p) else B_mat; 
  d <- dim(as.matrix(B))[2];  
  
  # Control Parameter Defaults
  control_args=control_list; control_names=names(control_args);
  wls_reg=if ( "wls_reg" %in% control_names ) control_args$wls_reg else 0;
  tol_val=if ( "tol_val" %in% control_names ) control_args$tol_val else 1e-7; 
  max_iter=if ( "max_iter" %in% control_names ) control_args$max_iter else 25; 
  max_iter_cg=if ( "max_iter_cg" %in% control_names ) control_args$max_iter_cg else 50; 
  #no more n needed
  init_stepsize_cg=if ( "init_stepsize_cg" %in% control_names ) control_args$init_stepsize_cg else rep(n,max_iter_cg); 
  beta_bt=if ( "beta_bt" %in% control_names ) control_args$beta_bt else 0.5;
  c_ag1=if ( "c_ag1" %in% control_names ) control_args$c_ag1 else 1e-3; #1e-3 too small?
  c_ag2=if ( "c_ag2" %in% control_names ) control_args$c_ag2 else 0.9;
  c_wolfe=if ( "c_wolfe" %in% control_names ) control_args$c_wolfe else 0; 
  max_iter_line=if ( "max_iter_line" %in% control_names ) control_args$max_iter_line else 100; 

  
  
  # Setting up the response and gradient functions ----
  if (ytype=='continuous'){
    # y.matrix should be m x n, m can be >= 1
    # This is just OPG, with the WLS as the exact solution per j
    
    
    mv_Y=y_matrix
    m=dim(y_matrix)[1]
    
    # j=1; test=T
    aD_j = function(j, test=F) {
      # centering data at obs j 
      Xj <- linearsdr:::matcenter_cpp(x_matrix,j,x0=NULL);
      B_Xj <- t(B)%*%Xj
      Vj <- rbind(rep(1,n), B_Xj)
      
      # Kernel Weights at j
      if (is.null(r_mat)) {
        Wj=linearsdr:::gauss_kern_cpp(Xj,bw) 
      } else { 
        rXj = t(r_mat)%*%Xj; 
        Wj=linearsdr:::gauss_kern_cpp(rXj,bw)
      }
      
      # Initial Value
      # WLS gives a (d+1)x(m-1) matrix; 
      # We want its transpose, a (m-1)x(d+1) matrix 
      c_j_ls=as.vector(t(linearsdr:::wls_cpp(Vj,mv_Y,Wj, reg=wls_reg))); 
      
      # Don't need to return the wls starting values for this  
      ## for least squares, we undo the vec to get (m-1)x(d+1), which is t(Aj)
      tA_hatj_ls <- matrix(c_j_ls, nrow = m, ncol = d+1) 
      
      a_hatj_ls <- tA_hatj_ls[,1]  
      D_hatj_ls <- t( tA_hatj_ls[,2:(d+1)] )  
      
      
      if (m > 1) {
        
        return( list( ahat=a_hatj_ls, Dhat=D_hatj_ls,
                      Dhat_ls=D_hatj_ls, 
                      weights=Wj) );
      } else { # i.e. m==1; 
        
        
        
        return( list( ahat=a_hatj_ls, Dhat=t(D_hatj_ls),
                      Dhat_ls=t(D_hatj_ls), 
                      weights=Wj) ) ;
      }
      
       
    }
    # hi=aD_j(1,T)  
  } else if (ytype %in% c("cat", "ord-cat") ) {
    # y.matrix should be 1 x n 
    
    m_classes=as.numeric(levels(as.factor(y_matrix)));
    m=length(m_classes); 
    mv_Y = linearsdr:::mnY_to_mvY( y_matrix, m_classes, ytype);
    if (ytype=="cat" ) { 
      linktype="expit";
      
      k_vec = colSums(mv_Y);
      mv_Y=matrix(mv_Y[1:(m-1),], m-1, n)
      
      # Empirical Logit Transform of the response
      link_mv_y=linearsdr:::emp_logit( mv_Y, k_vec, tune=0.05 ) ;
      
      
    } else if (ytype=="ord-cat" ) {
      linktype="ad-cat";
      
      k_vec = rep(1, n) #as.vector(y_matrix);
      mv_Y=matrix(mv_Y[2:(m),], m-1, n) # Drop the first row now cause its all 1
      # mv_Y[1,] # mv_Y[,1:20]
      # Empirical ad-cat Transform of the reponse
      link_mv_y=linearsdr:::emp_adcat( mv_Y, tune=0.05 ); #
      # emp_adcat( mv_Y[,980:1000], k_vec, tune=0.05 );
    } else if (ytype=="clogit" ) {
      linktype="clogit";
      
      k_vec = rep(1, n) #as.vector(y_matrix);
      mv_Y=matrix(mv_Y[2:(m),], m-1, n) # Drop the first row now cause its all 1
      # mv_Y[1,] # mv_Y[,1:20]
      # Empirical ad-cat Transform of the reponse
      link_mv_y=linearsdr:::emp_adcat( mv_Y, tune=0.05 ); #
      # emp_adcat( mv_Y[,980:1000], k_vec, tune=0.05 );
    }
    
    
    # j=1; test=T
    aD_j = function(j, test=F) {
      
      # centering data at obs j 
      Xj = linearsdr:::matcenter_cpp(x_matrix, index=j,x0=NULL); 
      B_Xj=t(B)%*%Xj;
      
      Vj=rbind(rep(1,n), B_Xj);  
      
      # Kernel Weights at j
      if (is.null(r_mat)) {
        # Wj=  exp(colSums(dnorm(Xj,0, sd=bw, log = T))) 
        # normalize_cpp( exp(colSums(dnorm( rbind(c(1,2), c(1,2),c(1,2) ) , 0, sd=bw, log = T))))
        # gauss_kern_cpp( rbind(c(1,2), c(1,2), c(1,2) ) ,bw)
        Wj=linearsdr:::gauss_kern_cpp(Xj,bw)
      } else { 
        rXj = t(r_mat)%*%Xj; 
        Wj=linearsdr:::gauss_kern_cpp(rXj,bw)
      }
      
      # Initial Value
      # WLS gives a (d+1)x(m-1) matrix; 
      # We want its transpose, a (m-1)x(d+1) matrix wls_reg
      
      c_j_ls=as.vector(t(linearsdr:::wls_cpp(Vj,link_mv_y,Wj, reg= wls_reg)));
      
      
      if (method=="cg") { 
        # Run Conjugate Gradients
        c_j_1=linearsdr:::aD_j_cg(c_j_ls, Vj, mv_Y, Wj, linktype, k_vec,
                      control_list=list(tol_val=tol_val,
                                        max_iter=max_iter_cg,
                                        init_stepsize=init_stepsize_cg,
                                        beta_bt=beta_bt,
                                        c_ag1=c_ag1,
                                        c_ag2=c_ag2,
                                        c_wolfe=c_wolfe,
                                        max_iter_line=max_iter_line ),
                      test);
        
        
        # mn_loss_j(c_j_ls, Vj, mv_Y, Wj, linktype, k_vec) 
        # 
        # mn_score_j(c_j_ls, Vj, mv_Y, Wj, link="clogit", k_vec) 
        # 
        # linearsdr:::dot_b_multinom(c_j_ls, 1, "expit")
          
        # c_j_1=aD_j_cg_test(c_j_ls, Vj, mv_Y, Wj, linktype, k_vec,
        #                    control_list=list(tol_val=tol_val,
        #                                      max_iter=max_iter_cg,
        #                                      init_stepsize=init_stepsize_cg,
        #                                      beta_bt=beta_bt,
        #                                      c_ag1=c_ag1,
        #                                      c_ag2=c_ag2,
        #                                      c_wolfe=c_wolfe,
        #                                      max_iter_line=max_iter_line ),
        #                    test)
        
      } else if (method=="newton") {
        # Run Newton-Raphson
        c_j_1=aD_j_newton(c_j_ls, Vj, mv_Y, Wj, linktype, k_vec, 
                          tol_val,max_iter, test);
      }
      
      
      # Don't need to return the wls starting values for this  
      ## for least squares, we undo the vec to get (m-1)x(p+1), which is t(Aj)
      
      tA_hatj=matrix(c_j_1, nrow = m-1, ncol = d+1);
      a_hatj=tA_hatj[,1]; D_hatj=t( tA_hatj[,2:(d+1)] );
      D_hatj_ls=t( matrix(c_j_ls, nrow = m-1, ncol = d+1)[,2:(d+1)] )
      
      if (m > 2) {
        return( list( ahat=a_hatj, Dhat=D_hatj,
                      Dhat_ls=D_hatj_ls,
                      weights=Wj) );
        # return( list( ahat=a_hatj ) );
      } else {
        return( list( ahat=a_hatj, Dhat=t(D_hatj),
                      Dhat_ls=t(D_hatj_ls),
                      weights=Wj) );
      }
      
    } # hi=aD_j(150,T)
    
  }
  # hi=aD_j(1, test=T)
  
  # Computing the candidate matrix ----
  # Use version with no parallel
  # and one with parallel
  
  if (parallelize) {
    
    # Compute estimates
    # aD_list = foreach::foreach(j = iterators::icount(n), 
    #                            .packages = "linearsdr" ) %dopar%{ 
    #                            # Release results
    #                            return(aD_j(j));
    #                            }
    aD_list = foreach::foreach(j = iterators::icount(n) ) %dopar%{ 
                                 # Release results
                                 return(aD_j(j));
                               }
  } else {
    aD_list = list();
    
    for(j in 1:n) {
      aD_list[[j]] = aD_j(j);
    }
    
  }
  

  
  ahat_list <- lapply( 1:n , FUN=function(j) aD_list[[j]][[1]] ) 
  Dhat_list <- lapply( 1:n , FUN=function(j) aD_list[[j]][[2]] )
  Dhat_ls_list <- lapply( 1:n , FUN=function(j) aD_list[[j]][[3]] )
  W_list <- lapply( 1:n , FUN=function(j) aD_list[[j]][[4]] )
  
  return(list( ahat=ahat_list,  
               Dhat=Dhat_list,
               Dhat_ls=Dhat_ls_list,
               weights=W_list) ) 
  
} 

# opcg_made(x_matrix, y_matrix, bw, B_mat=NULL, ytype="multinomial", #"ordinal",
#           method="cg", parallelize, r_mat, control_list)$Dhat

# opcg_made(x_matrix, y_matrix, bw, B_mat=NULL, ytype="ordinal",
#           method="newton", parallelize, r_mat, control_list)$Dhat

############### OPCG Candidate Matrix #########################
#' This is an internal function called by opcg 
#' 
#'
#'
#' @keywords internal
#' @noRd
#' 
opcg_DD <- function(x_matrix, y_matrix, bw, ytype='continuous', 
                    method="newton", parallelize=F, r_mat=NULL, 
                    control_list=list()) {
  
  
  # x_matrix=X; y_matrix=Y;#matrix(Y[2,],1,n)#Y;
  # bw;  ytype="ordinal"; #"continuous";#"multinomial";
  # tol_val= 1e-07; max_iter=25;
  # B_mat = NULL ; method="cg"; parallelize=T; r_mat=NULL; control_list=list();
  # control_list = list(); # B_mat=init_mat;
  
  # Rcpp::sourceCpp("../forwardsdr/src/opcg_wrap.cpp")
  
  DD_list=opcg_made(x_matrix, y_matrix, bw, B_mat=NULL, ytype, 
                    method, parallelize, r_mat, control_list); 
  # DD_mean=list_mean(list(matrix(0,3,2), matrix(1,3,2), matrix(2,3,2)))
  # DD_mean=list_mean(lapply(1:n, function(j) matrix(rowMeans(DD_list$Dhat[[j]]),ncol=1) ))
  
  DD_mean=list_mean(DD_list$Dhat);
  DD_mean_ls=list_mean(DD_list$Dhat_ls);
  
  # symmetrize the candidate matrices
  DD_mean = .5*(DD_mean + t(DD_mean));
  DD_mean_ls = .5*(DD_mean_ls + t(DD_mean_ls));
  
  return(list( opcg_mat=DD_mean,  
               wls_mat=DD_mean_ls, 
               gradients=DD_list$Dhat,
               weights=DD_list$weights) ) 
  
} 

############### OPCG wrapper #########################

opcg_wrap <- function(x_matrix, y_matrix, d, bw, ytype='continuous',
                 method="newton", parallelize=F, r_mat = NULL,
                 control_list=list()) {
  
  # Rcpp::sourceCpp("../forwardsdr/src/opcg_wrap.cpp")
  cand_mat=opcg_DD(x_matrix, y_matrix, bw, ytype, method, 
                   parallelize, r_mat , control_list)
  
  # x_matrix=X; y_matrix=Y; d=2; bw=4; ytype='multinomial';
  # method="cg"; lambda2a=1;lambda2b=5; parallelize=T; control_list=list();
  
  opcg_sdr=eigen_cpp(cand_mat$opcg_mat)$vec[,1:d]
  wls_sdr=eigen_cpp(cand_mat$wls_mat)$vec[,1:d]
  
  #' @return A list containing both the estimate and candidate matrix.
  #' \itemize{
  #'   \item opcg - 
  #'   \item opcg_wls - A 'pxd' matrix that estimates a basis for the central subspace based 
  #'   on the initial value of the optimization problem; useful for examining bad starting 
  #'   values.
  #'   \item cand_mat - A list that contains both the candidate matrix for OPCG and for
  #'   the initial value; this is used in other functions for order determination
  #'   \item gradients - The estimated local gradients; used in regularization of OPCG
  #'   \item weights - The kernel weights in the local-linear GLM. 
  #' }
  
  
  return( list( opcg=opcg_sdr,  
                # opcg_wls=wls_sdr,
                cand_mat=cand_mat,
                # ,
                gradients=cand_mat$gradients#,
                # weights=cand_mat$weights
  )  )
  
}


#' Outer Product of Canonical Gradients
#'
#' This implements the Outer Product of Canonical Gradients (OPCG) in a forth coming
#' paper Quach and Li (2021).
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
#'
#' @return A 'pxd' matrix that estimates a basis for the central subspace based on 
#' the estimated local gradients
#'  
#' @export
#' 
#' 
opcg <- function(x_matrix, y_matrix, d, bw, ytype='continuous',
                 method="newton", parallelize=F, r_mat = NULL,
                 control_list=list()) {
  
  beta=opcg_wrap(x_matrix, y_matrix, d, bw, ytype,
                     method, parallelize, r_mat,
                     control_list)$opcg 
  
  return(beta)
  
 
  
}



# A wrapper for Refined OPCG

ropcg=function(x_matrix, y_matrix, d, bw, ytype='continuous',
               method="newton", parallelize=F, r_mat = NULL,
               control_list=list()) {
  
  est0=opcg(x_matrix, y_matrix, d, bw, ytype,
           method, parallelize, r_mat,
           control_list)$opcg 
  
  refined_est0 = opcg(x_matrix, y_matrix, d, bw/4, ytype,
                     method, parallelize, r_mat = est0,
                     control_list)$opcg 
  
  for(iter in 1:25) {
    refined_est1 = opcg(x_matrix, y_matrix, d, bw/4, ytype,
                        method, parallelize, r_mat = refined_est0,
                        control_list)$opcg 
    
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
