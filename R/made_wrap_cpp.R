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

# OPCG-MADE Wrapper ####


############################## MADE wrappers ##################################


# x_matrix=X; y_matrix=Y; h=1.19; ytype="multinomial"; B_mat = B_1;
# method="newton"; parallelize = T; r_mat=NULL; control_list=list();
# ahat_list = aDhat$ahat;  Dhat_list = aDhat$Dhat;

made_iter <- function(x_matrix, y_matrix, bw,
                      ahat_list, Dhat_list, B_mat,
                      ytype='continuous',
                      method="newton", parallelize=F, r_mat=NULL,
                      control_list=list()) {
  # y.matrix should be m x n, with m depending on ytype

  # Supported ytypes are
  # continuous: Defaults to Squared Loss
  # multinomial: Multinomial-Categorical GLM with categorical link
  # ordinal: Multinomial-Ordinal GLM with ordinal link
  # other: Custom Loss Functions - to be done at a later date



  p=dim(x_matrix)[1]; n=dim(x_matrix)[2];
  B=B_mat; d=dim(B)[2];

  c_0 = as.vector(t(B));

  # Setting up the response and gradient functions ----
  if (ytype=='continuous'){
    # y.matrix should be m x n, m can be >= 1
    # This is just OPG, with the WLS as the exact solution per j


    mv_Y=y_matrix; m=dim(y_matrix)[1];


    B_list_j <- function(j, test=F) {
      Xj=matcenter_cpp(x_matrix, index=j,x0=NULL);

      # Kernel Weights at j
      if (is.null(r_mat)) {
        Wj=gauss_kern_cpp(Xj,bw)
      } else {
        Wj=gauss_kern_cpp( t(r_mat)%*%Xj,bw)
      }

      score_j=nm_score_j_made(c_0, Xj, mv_Y, Wj,
                              ahat=(ahat_list[[j]]), Dhat = Dhat_list[[j]]);
      info_j=nm_info_j_made(c_0, Xj, mv_Y, Wj,
                            ahat=(ahat_list[[j]]), Dhat = Dhat_list[[j]]);

      return( list( score=score_j, info=info_j) )

    }

  } else if (ytype %in% c( "multinomial", "ordinal")) {
    # y.matrix should be 1 x n

    m_classes=as.numeric(levels(as.factor(y_matrix)));
    m=length(m_classes);
    mv_Y = mnY_to_mvY( y_matrix, m_classes, ytype);
    if (ytype=="multinomial" ) {
      linktype="expit";

      k_vec = colSums(mv_Y);
      mv_Y=matrix(mv_Y[1:(m-1),], m-1, n)

      # Empirical Logit Transform of the response
      # link_mv_y=emp_logit( mv_Y, k_vec, tune=0.05 ) ;


    } else if (ytype=="ordinal" ) {
      linktype="culmit";

      k_vec = as.vector(y_matrix);
      mv_Y=matrix(mv_Y[1:(m-1),], m-1, n)

      # Empirical Culmit Transform of the reponse
      # link_mv_y=emp_culmit( mv_Y, k_vec, tune=0.05 );

    }

    # Don't need Empirical Logit Transform of the reponse
    # Empirical Logit Transform of the response
    # logit_Y=emp_logit( mv_Y, k_vec, tune=0.05 ) ;

    # j=10; test=T
    B_list_j <- function(j, test=F) {
      # centering data at obs j
      Xj=matcenter_cpp(x_matrix, index=j,x0=NULL);

      # Kernel Weights at j
      if (is.null(r_mat)) {
        Wj=gauss_kern_cpp(Xj,bw)
      } else {
        Wj=gauss_kern_cpp( t(r_mat)%*%Xj,bw)
      }

      # Loss functions and related functions
      loss_j = mn_loss_j_made(c_0, Xj, mv_Y, Wj, 
                       ahat=ahat_list[[j]], Dhat=Dhat_list[[j]],
                       link=linktype, k=k_vec);
      
      score_j=mn_score_j_made(c_0, Xj, mv_Y, Wj,
                              ahat=ahat_list[[j]], Dhat = Dhat_list[[j]],
                              link=linktype, k=k_vec);
      info_j=mn_info_j_made(c_0, Xj, mv_Y, Wj,
                            ahat=(ahat_list[[j]]), Dhat = Dhat_list[[j]],
                            link=linktype, k=k_vec) ;

      return( list( score=score_j, info=info_j, loss=loss_j) )

    } # B_list_j(1,T) #j=1;

  }; # else if {} for alternative objective functions; future work?


  # Computing the mean losses over index j ----
  # Use version with no parallel
  # and one with parallel

  # if ("x" %in% ls(envir = .GlobalEnv)) {
  #   get("x", envir = .GlobalEnv)
  # } else {
  #   x
  # }

  if (parallelize) {

    # Compute estimates
    B_list = foreach::foreach(j = iterators::icount(n),
                              .packages = "linearsdr") %dopar% {
                                # Release results
                                return(B_list_j(j))
                              }

  } else {

    B_list = list();

    for(j in 1:n) {
      B_list[[j]] = B_list_j(j);
    };

  }


  # loss <- Reduce('+', lapply( 1:n , FUN=function(j) B_list[[j]][[1]] ) )/n
  score=lapply( 1:n , FUN=function(j) B_list[[j]][[1]] );
  info=lapply( 1:n , FUN=function(j) B_list[[j]][[2]] );
  loss=lapply( 1:n , FUN=function(j) B_list[[j]][[3]] );
  
  c_1 = vecB_hat(c_0, score, info);
  # c_0  - solve(info)%*%score
  
  # loss_c0 = Reduce("+", loss)

  # return(list(estimate=c_1 # , loss=loss_0
  #             ))
  # return(list( est=c_1, loss_c0 = loss_c0));
  return(c_1)
}


#### MADE Block update ----
# aD_list=aDhat; print_B_iter=T
made_update = function(x_matrix, y_matrix, d, bw, aD_list ,B_mat,  ytype="continuous",
                       method, parallelize=F, r_mat=NULL,
                       control_list=list()) {
  # Control Parameter Defaults
  control_args=control_list; control_names=names(control_args);
  tol_val=if ( "tol_val" %in% control_names ) control_args$tol_val else 1e-7;
  max_B_iter=if ( "max_B_iter" %in% control_names ) {
    control_args$max_B_iter
  } else {
    25
  };
  print_B_iter=if ( "print_B_iter" %in% control_names ) control_args$print_B_iter else F;
  
  
  c_0=made_iter(x_matrix, y_matrix, bw,
                   ahat_list=aD_list$ahat, Dhat_list=aD_list$Dhat,
                   B_mat, ytype, method, parallelize, r_mat,
                   control_list);
  
  p=dim(x_matrix)[1];
  for(iter in 1:max_B_iter){
    B_0=t(matrix(c_0, nrow=d, ncol=p));
    
    # # aD-block update
    # aDhat=opcg_made(x_matrix, y_matrix, h, B_mat=B_0, ytype,
    #                 method, parallelize, r_mat=NULL,
    #                 control_list);
    # B-block update
    
    c_1=made_iter(x_matrix, y_matrix, bw,
                  ahat_list = aD_list$ahat, Dhat_list=aD_list$Dhat, 
                  B_0,
                  ytype, method, parallelize, r_mat=NULL,
                  control_list); 
    # B_1=t(matrix(c_1, nrow=d, ncol=p));
    
    # A Matrix distance of B_1 to B_0;
    # subspace_dist=mat_dist(B_1, B_0);
    
    # Vector Distance
    euc_dist=euc_norm_cpp(c_0 - c_1)/euc_norm_cpp(c_0);
    
    # Loss Distance 
    # loss_dist = loss_prev - loss_0;
    
    
    if(print_B_iter) print(c("B Update: euc dist is", euc_dist, iter));
    if( euc_dist < tol_val ) break();
    
    # The new B_0 for next iteration
    # B_0 = B_1;
    c_0=c_1; #loss_prev=loss_0;
    
    
  }
  
  return(c_1);
} 

#### MADE; An alternating Block Co-ordinate Descent Approach ----
# 
# x_matrix=X; y_matrix=Y; d; bw; B_mat=B_1; ytype="ordinal";
# method="newton"; parallelize=T; r_mat=NULL;
# control_list=list(print_iter=T, max_iter_made=5)

made <- function(x_matrix, y_matrix, d, bw, B_mat, ytype="continuous",
                 method="newton", parallelize=F, r_mat=NULL,
                 control_list=list()) {

  # Force one step iteration of OPCG or aD-block;
  # But proper Block Co-ordinate Descent requires fully minimizing at each step.
  # Maybe we can get away with alternating one step improvements if we 
  # enforce the descent property on the loss function, so that we are always 
  # descending
  # control_list$max_iter=1;

  # Control Parameter Defaults
  control_args=control_list; control_names=names(control_args);
  tol_val=if ( "tol_val" %in% control_names ) control_args$tol_val else 1e-7;
  max_iter_made=if ( "max_iter_made" %in% control_names ) {
    control_args$max_iter_made
  } else {
    25
  };
  print_iter=if ( "print_iter" %in% control_names ) control_args$print_iter else F;


  # Block step for aD parameter
  aDhat=opcg_made(x_matrix, y_matrix, bw, B_mat, ytype,
                  method, parallelize, r_mat=NULL,
                  control_list);

  # Block step for B parameter
  
  c_0=made_update(x_matrix, y_matrix, d, bw,
                  aD_list = aDhat,
                  B_mat, ytype, method, parallelize, r_mat,
                  control_list);

  p=dim(x_matrix)[1];
  for(iter in 1:max_iter_made){
    B_0=t(matrix(c_0, nrow=d, ncol=p));

    # aD-block update
    aDhat1=opcg_made(x_matrix, y_matrix, bw, B_mat=B_0, ytype,
                    method, parallelize, r_mat=NULL,
                    control_list);
    # B-block update
    
    c_1=made_update(x_matrix, y_matrix, d,bw,
                    aD_list = aDhat1, B_0,
                    ytype, method, parallelize, r_mat=NULL,
                    control_list);


    # B_1=t(matrix(c_1, nrow=d, ncol=p));
    # A Matrix distance of B_1 to B_0;
    subspace_dist=mat_dist(t(matrix(c_1, nrow=d, ncol=p)), B_0);
    # euc_dist=euc_norm_cpp(c_0 - c_1)/euc_norm_cpp(c_0);
    
    # Loss Distance 
    # loss_dist = loss(c_0) - loss(c_1);
    # subspace_dist; euc_dist;
    if(print_iter) print(c("MADE: subspace dist is", subspace_dist, iter));
    if( subspace_dist < tol_val ) break();
    if(iter==max_iter_made) print("0 - non-convergence")
    # The new B_0 for next iteration
    # B_0 = B_1;
    c_0=c_1;

  }

  B_hat_made = t(matrix(c_1, nrow=d, ncol=p));
  # B_hat_made = matrix(c_1, nrow=p, ncol=d)
  return(B_hat_made)
  # return(list( estimate=B_hat_made, loss=B_hat_made_loss))
};



