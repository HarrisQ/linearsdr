###############################################################################
################ RADE and Regularized OPCG - Ridge and LASSO ##################
###############################################################################

library(glmnet)

rade_alt<-function(x_matrix, y_matrix, d, bw, ytype='continuous', 
                   method="newton", parallelize=F, 
                   l2_pen=0, l1_pen=0, control_list=list() ){
  
  # Control Parameter Defaults
  control_args=control_list; control_names=names(control_args);
  tol_val=if ( "tol_val" %in% control_names ) control_args$tol_val else 1e-7; 
  max_iter=if ( "max_iter" %in% control_names ) control_args$max_iter else 25 ; 
  print_iter=if ( "print_iter" %in% control_names ) control_args$print_iter else F ; 
  
  # Setting parameters
  n = length(grad_list); d = dim(init_est)[2];
  eye_d = diag(1,d,d);
  
  # l2_pen=1; l1_pen=.05;
  
  # Given initial Gamma, compute nu_s (i.e. Riesz Representers I think)
  nu_est0 = foreach (j=1:n, .packages="linearsdr") %dopar% { # Can parallelize
    return( ( (1+l2_pen)^(-1) )*
              solve_cpp( t(init_est)%*%init_est, t(init_est)%*%grad_list[[j]]) ) ;
  }
  # nu_est0 = list()
  # for (j in 1:n) { # Can parallelize
  #   nu_est0[[j]] = solve_cpp( t(init_est)%*%init_est, 
  #                             t(init_est)%*%grad_list[[j]] ) ;
  # }
  
  
  # Given nu estimates, compute gamma estimate 
  vvt0=list_sum(nu_est0,nu_est0);
  # inv_term0=kronecker( vvt0 ,sigx) + diag(n*l2_pen, d*p,d*p) 
  # ninv_term0 = matrix(c( sigx%*%list_sum(grad_list,nu_est0 ) ) )
  # gamma_hat0 = solve_cpp(inv_term0, ninv_term0);
  # Gamma_hat0 = matrix(gamma_hat0, p, d);
  
  
  # What about estimate Gamma row by row??? This would allow for l1 and l2 reg like
  # in Lexin Li's papers.
  # Gamma_hat0 = matrix(0, p, d)
  # for (i in 1:p) { # Can parallelize
  #   gamma_i = solve_cpp(vvt0 + diag(l2_pen, d, d), 
  #                       matrix(list_sum(grad_list,nu_est0 )[i,] ) );
  #   Gamma_hat0[i,] = gamma_i;
  # }
  Gamma_hat0=foreach(i=1:p, .combine='rbind', .packages='linearsdr') %dopar% {
    gamma_i = solve_cpp(vvt0, # + diag(l2_pen, d, d),
                        matrix(list_sum(grad_list,nu_est0 )[i,] ) );
    return(t( ((1+l2_pen)^(-1))* gamma_i));
  }
  
  
  # LASSO update for Gamma
  
  if(l1_pen > 0){
    # i=1
    vvt0_rt=matpower_cpp(vvt0, .5);
    x_star = lapply(1:p, function(i) rbind( vvt0_rt, sqrt(l2_pen)*eye_d ))
    y_star = lapply(1:p, function(i) c( vvt0_rt%*%Gamma_hat0[i,], rep(0,d)) )
    
    Gamma_hat0_l0=foreach(i=1:p, .combine = 'rbind', .packages =  "glmnet" ) %dopar% {
      gamma_hat0_l0_i=as.vector(glmnet(x_star[[i]], y_star[[i]], intercept=F,
                                       family="gaussian", alpha=1, lambda=l1_pen)$beta) 
      
      return(gamma_hat0_l0_i);
    }
    Gamma_hat0 = Gamma_hat0_l0;
  }
  
  
  loss0 = foreach(j=1:n, .combine='+') %dopar% { # Can parallelize
    diff_term = c(grad_list[[j]]  -  Gamma_hat0%*%nu_est0[[j]]);
    return(t(diff_term)%*%diff_term/n); 
  }
  
  
  Gam0_norml2 = Reduce("+", lapply(1:n, function(j) t(c(Gamma_hat0%*%nu_est0[[j]]))%*%
                                     c(Gamma_hat0%*%nu_est0[[j]]) ) ) ;
  Gam0_norml1 = t(c(Gamma_hat0))%*%c(Gamma_hat0);
  # loss0 = loss0 + l2_pen*Gam0_norm + l1_pen*sqrt(Gam0_norm) 
  loss0 = loss0 + l2_pen*Gam0_norml2 + l1_pen*sqrt(Gam0_norml1)
  
  # grad_list[[1]]  -  Gamma_hat0%*%nu_est0[[1]]
  # for (j in 1:n) { # Can parallelize
  #   loss0 = loss0 +
  #     ( t(c( grad_list[[j]]  -  Gamma_hat0%*%nu_est0[[j]] ))%*%
  #         (c( grad_list[[j]]  - Gamma_hat0%*%nu_est0[[j]] ))
  #     )/n 
  # }
  
  
  for(iter in 1:max_iter){
    
    # update nu_est
    nu_est1 = foreach (j=1:n, .packages="linearsdr") %dopar% { # Can parallelize
      return(((1+l2_pen)^(-1))*
               solve_cpp( t(Gamma_hat0)%*%Gamma_hat0, t(Gamma_hat0)%*%grad_list[[j]]) ) ;
    }
    
    # nu_est1 = list()
    # for (j in 1:n) {
    #   nu_est1[[j]] = solve(t(Gamma_hat0)%*%Gamma_hat0)%*%
    #     (t(Gamma_hat0)%*%grad_list[[j]]);
    # }
    
    # update Gamma est
    
    vvt1=list_sum(nu_est1, nu_est1 );
    # inv_term1=kronecker( vvt1 ,sigx) + diag(n*l2_pen, d*p,d*p)
    # ninv_term1 = matrix(c( sigx%*%list_sum(grad_list,nu_est1 ) ) )
    # 
    # gamma_hat1 = solve_cpp(inv_term1, ninv_term1);
    # Gamma_hat1 = matrix(gamma_hat1, p, d)
    
    Gamma_hat1=foreach(i=1:p, .combine='rbind', .packages='linearsdr') %dopar% {
      gamma_i = solve_cpp(vvt1 + l2_pen*eye_d, matrix(list_sum(grad_list,nu_est0 )[i,] ) );
      return(t( ((1+l2_pen)^(-1))*gamma_i ));
    }
    # Gamma_hat1 = matrix(0, p, d)
    # for (i in 1:p) {
    #   gamma_i = solve(vvt1 + diag(l2_pen, d, d) )%*%list_sum(grad_list,nu_est1 )[i,]
    #   Gamma_hat1[i,] = gamma_i
    # }
    
    # LASSO update for Gamma
    
    # library(glmnet)
    if(l1_pen > 0) {
      # i=1
      vvt1_rt = matpower_cpp(vvt1, .5);
      x_star = lapply(1:p, function(i) rbind( vvt1_rt, sqrt(l2_pen)*eye_d ))
      y_star = lapply(1:p, function(i) c( vvt1_rt%*%Gamma_hat1[i,], rep(0,d)) )
      
      Gamma_hat1_l0=foreach(i=1:p, .combine = 'rbind', .packages = "glmnet" ) %dopar% {
        gamma_hat1_l0_i=as.vector(glmnet(x_star[[i]], y_star[[i]], intercept=F,
                                         family="gaussian", alpha=1, lambda=l1_pen)$beta)
        # gamma_hat1_l0_i=as.vector(l1ce(y_star[[i]] ~ x_star[[i]] - 1, 
        #                                sweep.out=NULL, standardize=F, bound=l1_pen, absolute.t=T)$coef)
        return(gamma_hat1_l0_i)
      }
      Gamma_hat1 = Gamma_hat1_l0;
    }
    
    loss1 = foreach(j=1:n, .combine='+') %dopar% { # Can parallelize
      diff_term = c(grad_list[[j]]  -  Gamma_hat1%*%nu_est1[[j]]);
      return(t(diff_term)%*%diff_term/n); 
    }
    
    Gam1_norml2 = Reduce("+", lapply(1:n, function(j) t(c(Gamma_hat1%*%nu_est1[[j]]))%*%
                                       c(Gamma_hat1%*%nu_est1[[j]]) ) )
    Gam1_norml1 = t(c(Gamma_hat1))%*%c(Gamma_hat1);
    
    loss1 = loss1 + l2_pen*Gam1_norml2 + l1_pen*sqrt(Gam1_norml1) 
    
    # Gam1_norm=t( c(Gamma_hat1))%*%c(Gamma_hat1)
    # loss1 = l2_pen*Gam1_norm + l1_pen*sqrt( Gam1_norm )
    # for (j in 1:n) { # cov(t(X))%*% # sig_root = matpower_cpp(cov(t(X)), 0.5);
    #   loss1 = loss1 +
    #     ( t(c( grad_list[[j]]  -  Gamma_hat1%*%nu_est1[[j]] ))%*%
    #         (c( grad_list[[j]]  -  Gamma_hat1%*%nu_est1[[j]] )) )/n 
    #   
    # }
    
    loss_dist = loss0 - loss1;
    subspace_dist = mat_dist(Gamma_hat0, Gamma_hat1);
    
    if(print_iter) print(c(loss_dist, subspace_dist, iter))
    if(loss_dist < tol_val) break;
    
    Gamma_hat0=Gamma_hat1; loss0 = loss1;
    
  }
  
  Gamma_hat0=apply(Gamma_hat0, 2, normalize_cpp)
  Gamma_hat0[,1]=normalize_cpp(Gamma_hat0[,1])
  
  return(Gamma_hat0)
} 






