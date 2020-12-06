#tester function

# Create a wrapper for the function
gauss_mean<- function(k){
  Sys.sleep(5) # Forces a 5 second wait 
  data<-rnorm(100000,mean=10*k,sd=sqrt(5)) # Draw 100k samples from N(mu,5)
  # pdf(paste("Histogram",filetitle,".pdf",sep="")) # Create Histogram PDF
  # hist(data, main =paste("Histogram: mean=",10*filetitle,sep=""))
  # dev.off()
  return(mean(data))
  
}

gauss_mean_cpp<- function(k){
  Sys.sleep(5) # Forces a 5 second wait 
  data<-rnorm(100000,mean=10*k,sd=sqrt(5)) # Draw 100k samples from N(mu,5)
  # pdf(paste("Histogram",filetitle,".pdf",sep="")) # Create Histogram PDF
  # hist(data, main =paste("Histogram: mean=",10*filetitle,sep=""))
  # dev.off()
  return(mean(center_cpp(data)))
  
}

test_function = function(t) {
  mean_list=foreach::foreach(k=1:t, .packages = "linearsdr") %dopar% {gauss_mean(k)}
  return(Reduce('+',mean_list)/t)
}

test_function_cpp = function(t) {
  mean_list=foreach::foreach(k=1:t, .packages = "linearsdr") %dopar% {gauss_mean_cpp(k)}
  return(Reduce('+',mean_list)/t)
}

############################ 

# Intercept and Gradient estimation

aD_j = function(j, mv_y_matrix, x_matrix, h, B,
                link="continuous", mvlink_y=NULL, k_vec=NULL,
                r_mat=NULL, tol_val= 1e-07, max_iter=25, test=F) {
  
  n = dim(x_matrix)[2]; d <- dim(B)[2];
  m=dim(mv_y_matrix)[1] + 1;
  
  # centering data at obs j 
  Xj = matcenter_cpp(x_matrix, index=j,x0=NULL); #(x_matrix - x_matrix[,j])
  B_Xj=t(B)%*%Xj;
  Vj=rbind(rep(1,n), B_Xj);
  
  # Kernel Weights at j
  if (is.null(r_mat)) {
    Wj=gauss_kern_cpp(Xj,h) 
  } else { 
    rXj = t(r_mat)%*%Xj; 
    Wj=gauss_kern_cpp(rXj,h)
  }
  
  if(link=="continuous"){
    c_j_ls=wls_cpp(Vj,mv_Y,Wj);
    
    
    # Don't need to return the wls starting values for this  
    ## for least squares, we undo the vec to get (m-1)x(d+1), which is t(Aj)
    tA_hatj_ls <- matrix(c_j_ls, nrow = m, ncol = d+1) 
    
    a_hatj_ls <- tA_hatj_ls[,1]  
    D_hatj_ls <- t( tA_hatj_ls[,2:d] )  
    
    
    if (m > 1) {
      return( list( ahat=a_hatj_ls,  
                    Dhat=D_hatj_ls ) );
    } else {
      return( list( ahat=a_hatj_ls,  
                    Dhat=t(D_hatj_ls) ) );
    }
  } else if (link %in% c("expit", "culmit") ) {
    # Initial Value
    # WLS gives a (d+1)x(m-1) matrix; 
    # We want its transpose, a (m-1)x(d+1) matrix 
    c_j_ls=as.vector(t(wls_cpp(Vj,mvlink_y,Wj)));
    
    # Run Newton-Raphson
    c_j_1=aD_j_newton(c_j_ls, Vj, mv_y_matrix, Wj, link, k_vec, 
                      tol_val,max_iter, test);
    
    # Don't need to return the wls starting values for this  
    ## for least squares, we undo the vec to get (m-1)x(d+1), which is t(Aj)
    tA_hatj=matrix(c_j_1, nrow = m-1, ncol = d+1);
    a_hatj=tA_hatj[,1]; D_hatj=t( tA_hatj[,2:(d+1)] );
    D_hatj_ls=t( matrix(c_j_ls, nrow = m-1, ncol = d+1)[,2:(d+1)] )
    
    if (m > 2) {
      return( list( ahat=a_hatj, Dhat=D_hatj,
                    Dhat_ls=D_hatj_ls) );
    } else {
      return( list( ahat=a_hatj, Dhat=t(D_hatj),
                    Dhat_ls=t(D_hatj_ls) ) );
    }
    
  } 
  
} # aD_j(1,T)


aD_list_fn = function(k, Y, X, h, export=F){
  
  m_classes=as.numeric(levels(as.factor(Y)));
  mv_Y = mnY_to_mvY( Y, m_classes, ytype="multinomial");
  m=length(m_classes);
  k_vec = colSums(mv_Y);
  mv_Y=matrix(mv_Y[1:(m-1),], m-1, n)
  # Empirical Logit Transform of the response
  logit_Y=emp_logit( mv_Y, k_vec, tune=0.05 ) ;
  
  
  B = diag(1,p,p) 
  if (export) {
    foreach::foreach(j =1:k, 
                     .packages = "linearsdr",
                     .export = c("mv_Y", "X", "bw","B","logit_Y","k_vec")) %dopar%{ 
                       
                       
                       
                       aD_k=aD_j(j, mv_Y, X, h, B, 
                                 link="expit", mvlink_y=logit_Y, k_vec=k_vec,
                                 r_mat=NULL, tol_val= 1e-07, max_iter=25, test=F)
                       
                       # Release results
                       return(aD_k);
                     }
  } else {
    foreach::foreach(j = 1:k, 
                     .packages = "linearsdr" ) %dopar%{ 
         
                       aD_k=aD_j(j, mv_Y, X, h, B, 
                                 link="expit", mvlink_y=logit_Y, k_vec=k_vec,
                                 r_mat=NULL, tol_val= 1e-07, max_iter=25, test=F)
                       
                       # Release results
                       return(aD_k);
                     }
  }

}
