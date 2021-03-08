#tester function

# Testing out hCG
library('Rcpp')
library('RcppArmadillo')
library('inline')

p=10; n=1000
beta=matrix(1:30, 10, 3)
c_init = c(matrix(rep(1,30), 10,3))


X = matrix(rnorm(n*p), nrow=p, ncol=n)
Y = t(beta)%*%(X)# + rnorm(n)

sq_loss(c_init, X, Y)
sq_grad(c_init, X, Y)

control_list1=list(tol_val=1e-7, max_iter=50, 
                   init_stepsize=rep(1,50), beta_bt=.5,
                   c_ag=1e-3, c_ag2=9e-1, c_wolfe=9e-1,
                   max_iter_line=100)
test=T

# c_star=
aD_j_cg_test(c_init, X, Y, wj = 1, link='hi', k=1,  control_list1, test)

t(matrix(aD_j_cg_test(c_init, X, Y, wj = 1, link='hi', k=1,  control_list1, test), 3,10))

# This works 
c_star=aD_j_cg_test2(c_init, X, Y, control_list1, test)
t(matrix(c_star, 3,10))

aD_j_cg_test2 =function(init,
                  vj,
                  y_datta,
                  control_list,
                  test ) {
  # init=c_init
  # vj=X
  # y_datta=Y
  # control_list=control_list1
  
  tol_val = control_list$tol_val;
  max_iter = control_list$max_iter;
  
  s=control_list$init_stepsize; 
  beta_bt=control_list$beta_bt; #//backtracking beta
  c_ag=control_list$c_ag; #//parameter for Armijo condition, delta in Hybrid-DY
  c_ag2=control_list$c_ag2; #//parameter for second Armijo bound in nesterov 
  c_wolfe=control_list$c_wolfe; #//parameter for curvature condition, sigma in H-DY
  max_iter_line=control_list$max_iter_line;
  
  
  #// Step 0: Set the initial value, and compute the loss and score
  #// Also set the initial p_0 to the gradient
  
  p= dim(vj)[1]; m=dim(y_datta)[1];n=dim(y_datta)[2];
  c_mat = matrix(init, p, m); 
  c_now = init;
  y_vec = c(y_datta);
  y_vec
  
  
  nll_now = t(y_vec - kronecker(t(vj), diag(1,m,m))%*%c_now )%*%
    (y_vec - kronecker(t(vj), diag(1,m,m))%*%c_now )/n; #mn_loss_j(c_now,vj,y_datta,wj,link,k);
  grad_now = -2*t(kronecker(t(vj),diag(1, m,m), FUN="*"))%*%
              (y_vec - kronecker(t(vj), diag(1,m,m))%*%c_now)/n;
  p_now = grad_now; 

  #arma::vec c_next;
  #arma::uword iter;
  for (iter in 1:max_iter) { 
    # iter=1
    s_now = s[iter]; 
    #// Step 1: Line search 
    m_ag=0;
    # m_cg;
    for (m_cg in 1:max_iter_line) {
      # m_cg=1
      # // the Armijo rule
      # // the computational effort is not just dependent on the m_cg search
      # // it also depends on the step-size; large s_now means the m_cg will have
      # // to search farther.
      # // So we would like to have the smallest s_now to allow for fastest m_cg find. 
      armijo_bound = (c_ag*(beta_bt^m_cg)*s_now*
                                         ( t(p_now)%*%grad_now));
      
      
      
      # // evaluation for armijo condition
      # arma::vec c_search; 
      c_search = c_now - (beta_bt^m_cg)*s_now*p_now ;
      
      # double suff_dec_ag;
      suff_dec_ag = t(y_vec - kronecker(t(vj), diag(1,m,m))%*%c_now )%*%
        (y_vec - kronecker(t(vj), diag(1,m,m))%*%c_now )/n - 
        t(y_vec - kronecker(t(vj), diag(1,m,m))%*%c_search )%*%
        (y_vec - kronecker(t(vj), diag(1,m,m))%*%c_search )/n;
      armijo_cond = as.numeric(suff_dec_ag >= armijo_bound);
      
      
      #int 
      armijo_cond2=0;
      if (c_ag2 > 0) {
        #// the second bound in armijo-goldstein in nesterovs intro to conv opt text
        armijo_bound2 = (c_ag2*(beta_bt^m_cg)*s_now*
                                           ( t(p_now)%*%grad_now));
        #// second sufficient descent bound uses the same suff_dec_ag
        armijo_cond2 =+ (suff_dec_ag <= armijo_bound2);
      }
      
      wolfe_cond=0;
      if (c_wolfe > 0) {
        # // the weak Wolfe condition
        wolfe_bound = c_wolfe*( t(p_now)%*%grad_now);
        
        # // evaluation for curvature in weak wolfe
        # double curv_wolfe;
        curv_wolfe = t(p_now)%*%(-2*t(kronecker(t(vj),diag(1, m,m), FUN="*"))%*%
                                   (y_vec - kronecker(t(vj), diag(1,m,m))%*%c_search)/n) ;
        
        wolfe_cond = as.numeric(curv_wolfe <= wolfe_bound);
      }  
      
      
      if ( armijo_cond + armijo_cond2== 2 ) {# //+ wolfe_cond 
        m_ag = m_cg;
        break;
      }
      
    }
    
    h_now =  (beta_bt^m_ag)*s_now ;
    c_next = c_now - h_now*p_now;
    
    
    #// #Step 2a: Compute Loss;
    # arma::mat nll_next(1,1); 
    nll_next=t(y_vec - kronecker(t(vj), diag(1,m,m))%*%c_next )%*%
      (y_vec - kronecker(t(vj), diag(1,m,m))%*%c_next )/n;    
    
    # double nll_dist; 
    nll_dist = ( nll_now - nll_next);
    
    if (test) {
      # // Rprintf("Printing: iter %iter, ll Dist %ll_dist, eu Dist %eu_dist ", 
                 # //         iter, ll_dist, eu_dist);
      # Rcout << "Printing nll_dist, iter: " << nll_dist<< ", "  << iter << "\n";
      print(paste("Printing nll_dist: ", as.character(nll_dist), ", iter",as.character(iter)) )
    }
    
    if( nll_dist < tol_val) {
      break;
    } else {
      
      # // #Step 2b: Compute gradient;
      # arma::vec
      grad_next = (-2*t(kronecker(t(vj),diag(1, m,m), FUN="*"))%*%
                     (y_vec - kronecker(t(vj), diag(1,m,m))%*%c_next)/n)

      # // Step 3: Compute the coeffiecient
      # // Fletcher-Reeves
      # // double beta_cg_fr = as_scalar( ( grad_next.t()*grad_next )/
      # //                                ( grad_now.t()*grad_now ) );

      # // Dai-Yuan
      beta_cg_dy = -( ( t(grad_next)%*%grad_next )/
                                         ( t(p_now)%*%(grad_next-grad_now) ) );

      # // Hestenes-Stiefel
      beta_cg_hs = -( ( t(grad_next)%*%(grad_next-grad_now) )/
                                        ( t(p_now)%*%(grad_next-grad_now) ) );

      # // Hybird
        # arma::vec beta2(2); beta2(0) = beta_cg_dy; beta2(1)=beta_cg_hs;
        # arma::vec beta3(2); beta3(0) = 0; beta3(1)=min(beta2);
      beta_cg_hybrid=max(0, min(beta_cg_dy, beta_cg_hs));

      # // Step 4: Update p
      p_next = grad_next - beta_cg_hybrid*p_now;

      #  // Update all inputs for next iteration
      c_now=c_next;
      nll_now=nll_next;
      grad_now=grad_next;
      p_now=p_next;
    }
  }
  
  return(c_next);
  
}

#########################################################################
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
