###---------------------------------------------------------------
###| A Generalized Eigenproblem Solver                           |
###---------------------------------------------------------------

#* Lexin Li's approach requires taking the root of M and G, 
#* which can be ill-posed when p >> n, so we need a way around it

beta_hat = function(mat, mat_A, d, lambda2b, metric_mat=NULL) {
  # mat=cand_mat$opcg_mat; 
  # mat_A = init_mat;
  # d=3; lambda2b = 1; metric_mat = NULL;
  
  # Setting some parameters
  p = dim(mat)[1];
  
  # Defining the metric matrix
  G = if (is.null(metric_mat)) diag(1,p,p) else metric_mat;
  
  mat_root = matpower_cpp(mat, 0.5);
  m_star = rbind( mat_root, sqrt(lambda2b)*matpower_cpp(G, 0.5)  );

  est_list=lapply(1:d, function(j) {
    #j=1
    
    u_star = matrix( c(mat_root%*%mat_A[,j], rep(0,p)), nrow=1, ncol=2*p );

    theta_hat = wls_cpp(x_matrix =  t(m_star), y_matrix = u_star,
                        weights = rep(1,2*p), reg = 0 )
    
    return(theta_hat);
  }) 
  
  return(Reduce('cbind', est_list))
}  
  
beta_hat2 = function(mat, mat_A, d, lambda2b, metric_mat=NULL) {
  # mat=cand_mat$opcg_mat; 
  # mat_A = init_mat;
  # d=3; lambda2b = 1; metric_mat = NULL;
  
  # Setting some parameters
  p = dim(mat)[1];
  
  # Defining the metric matrix
  G = if (is.null(metric_mat)) diag(1,p,p) else metric_mat;
  
  mat_root = matpower_cpp(mat, 0.5);
  m_star = rbind( mat_root, sqrt(lambda2b)*matpower_cpp(G, 0.5)  );
  
  est_list=lapply(1:d, function(j) {
    #j=1
    
    u_star = matrix( c(mat_root%*%mat_A[,j], rep(0,p)), nrow=1, ncol=2*p );
    
    theta_hat = wls_cpp(x_matrix =  t(m_star), y_matrix = u_star,
                        weights = rep(1,2*p), reg = 0 )
    
    return(theta_hat);
  }) 
  
  return(Reduce('cbind', est_list))
}  

alpha_hat = function(mat, mat_B, d, metric_mat=NULL) {
  # mat_B=
  
  p = dim(mat)[1];
  G = if (is.null(metric_mat)) diag(1,p,p) else metric_mat;
  
  G_root = matpower_cpp(G, -0.5)
  alpha_svd=svd(G_root%*%mat%*%mat_B)
  
  a_hat = G_root%*%alpha_svd$u%*%t( alpha_svd$v )
  
  return(a_hat)
}

reg_eigen = function(mat, init_mat, d, lambda2b, metric_mat=NULL, 
                     max_iter=25, tol_val=10e-7, print_iter=F) {
  
  #' mat is the positive matrix you want to decompose
  #' d is the number of eigenvectors to recover
  #' metric_mat is the matrix defining the inner product
  #' to consider. The default is NULL, which corresponds to the regular
  #' euclidean inner product
  
  # mat=cand_mat$opcg_mat; 
  # init_mat; d; lambda2b=20; metric_mat=cov(t(x_matrix))+diag(.05,p,p); 
  # max_iter=25; tol_val=10e-7; print_iter=T
  
  
  G= if (is.null(metric_mat)) diag(1,p,p) else metric_mat ;
  M= mat 
  # Given init_mat as mat_A to start
  b_hat0 = beta_hat(M, mat_A = init_mat, d, lambda2b, metric_mat=G)
  
  # Given this beta_hat
  a_hat0 = alpha_hat(M, mat_B = b_hat0, d, metric_mat=G)
  
  for(iter in 1:max_iter) {
    
    b_hat1 = beta_hat(M, mat_A = a_hat0, d, lambda2b, metric_mat=G);
    
    beta_dist = mat_dist(b_hat0, b_hat1);
    
    if(print_iter) print(c("B Update: loss dist is", beta_dist, iter));
    if( beta_dist < tol_val ) break();
    
    a_hat1 = alpha_hat(M, mat_B = b_hat1, d, metric_mat=G);
    # The new B_0 for next iteration
    # B_0 = B_1;
    a_hat0 = a_hat1; b_hat0 = b_hat1; 
  }
  
  b_hat = Reduce('cbind',
                 lapply(1:d, function(j) 
                   c(b_hat1[,j])/sqrt( t(c(b_hat1[,j]) )%*%c(b_hat1[,j]) ) ));
  
  # return the eigenvectors
  return(b_hat)
  
}



reg_eigen2 = function(mat, init_mat, d, gev_reg, metric_mat=NULL, 
                      max_iter=25, tol_val=10e-7, print_iter=F) {
  
  #' mat is the positive matrix you want to decompose
  #' d is the number of eigenvectors to recover
  #' metric_mat is the matrix defining the inner product
  #' to consider. The default is NULL, which corresponds to the regular
  #' euclidean inner product
  
  # mat=opcg_obj[["cand_mat"]][["opcg_mat"]];
  # init_mat = eigen_cpp( opcg_obj[["cand_mat"]][["opcg_mat"]])$vec[,1:d];
  # d; gev_reg = 10;
  # metric_mat = cov(t(X)) + diag(reg, p, p)
  # max_iter=25; tol_val=10e-7; print_iter=T
  
 
  G=if (is.null(metric_mat)) diag(1,p,p) else metric_mat ;
 
  # Given init_mat as mat_A to start
  b_hat0_list = lapply(1:d, function(j) ginv(mat + gev_reg*G)%*%mat%*%init_mat[,j] )
  b_hat0 = Reduce('cbind', b_hat0_list)
  # b_hat0 = beta_hat(mat, mat_A = init_mat, d, lambda2b, metric_mat=G)
  
  # Given this beta_hat
  a_hat0 = alpha_hat(mat, mat_B = b_hat0, d, metric_mat=G ) 
  
  
  for(iter in 1:max_iter) {
    # iter=1
    b_hat1_list = lapply(1:d, function(j) ginv(mat + gev_reg*G)%*%mat%*%a_hat0[,j] )
    b_hat1 = Reduce('cbind', b_hat1_list)
    # b_hat1 = beta_hat(mat, mat_A = a_hat0, d, lambda2b, metric_mat=G)

    beta_dist = mat_dist(b_hat0, b_hat1);
    
    if(print_iter) print(c("B Update: loss dist is", beta_dist, iter));
    if( beta_dist < tol_val ) break();
    
    a_hat1 = alpha_hat(mat, mat_B = b_hat1, d, metric_mat=G ) 
    
    # The new B_0 for next iteration
    # B_0 = B_1;
    a_hat0 = a_hat1; b_hat0 = b_hat1; 
  }
  
  b_hat = Reduce('cbind',
                 lapply(1:d, function(j) 
                   c(b_hat1[,j])/sqrt( t(c(b_hat1[,j]) )%*%c(b_hat1[,j]) ) ));
  
  # return the eigenvectors
  return(b_hat)
  
}

 