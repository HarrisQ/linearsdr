#######################################################################
#           Boxcox transformation to multivariate normality
#                             Modified from 
#           Sufficient Dimension Reduction with R (Li, 2018)
#######################################################################

# Not ready for use

#######################################################################
#                  function 1: evaluate likelihood
#######################################################################


#'
#' @keywords internal
#' @noRd
#' 
#' 
likelihood = function(x,lam,eps, xmar, paralleize=F){ 
  n=dim(x)[1]; 
  # xmar=standmar(x,eps); # moved to a higher level to avoid redundancy
  xlam=bocotranmat(x,lam,eps); 
  mu=colMeans(xlam); 
  sum_lam_log_xmar=rowSums(t((lam - 1) * t( log(xmar) ) ) ) 
  sig=var(xlam); 
  inv_sig=ginv(sig);
  constant_1 <- (-1/2) * log( abs( det(sig)) ); 
  ll_func <- function(i) { # i <- 1
    (-1/2)*(xlam[i,]-mu)%*%inv_sig%*%(xlam[i,]-mu) + constant_1 + sum_lam_log_xmar[i]  
  }
  if ( !paralleize ) {
    loglik <- sum ( sapply(1:n, ll_func ) )
  } else if (paralleize){   
    
    loglik <- foreach(i=1:n, .combine='+', .packages="linearsdr") %dopar% {
      return(ll_func(i))  
    }
  }
  return( loglik )
}


#######################################################################
#                  function 2: standardize a matrix
#                              treating each row as a random vector
#                              in an iid sample
#######################################################################
#'
#' @keywords internal
#' @noRd
#' 
standmat = function(x){
  mu = colMeans(x) 
  signrt = matpower_cpp(var(x),-1/2)
  return(t(t(x) - mu)%*%signrt)
}

#######################################################################
#                  function 3: marginally standardize a matrix
#                              treating each row as a random vector
#                              in an iid sample, and then transform 
#                              it to be positive
#######################################################################
#'
#' @keywords internal
#' @noRd
#' 
standmar = function(x,eps){
  mu = colMeans(x)
  # sig = diag(diag(var(x)))
  signrt = matpower_cpp(var(x),-1/2)
  x1 = t(t(x) - mu)%*%signrt
  x2 = t(t(x1) - apply(x1,2,min))+eps
  return(x2)
}
#######################################################################
#                  function 4: standardize a vector
#######################################################################
#'
#' @keywords internal
#' @noRd
#' 
standvec = function(x) (x - mean(x))/sd(x); 

#######################################################################
#                  function 5: boxcox transform a vector
#                              x is n dim vector
#                              lam is box cox parameter
#                              eps is distance above 0
#######################################################################
#'
#' @keywords internal
#' @noRd
#' 
bocotranvec = function(x,lam,eps){
  n = length(x)
  x1 = standvec(x)
  x2 = x1 - min(x1) + eps
  if(abs(lam)< 10^(-10)) 
  x3 = log(x2) else
  x3 = (x2^lam - 1)/lam
  return(x3)
}
#######################################################################
#                  function 6: boxcox transform a matrix
#                              x is n x p matrix
#                              lam is p vector
#                              eps is 0.5 say
#######################################################################
#'
#' @keywords internal
#' @noRd
#' 
bocotranmat = function(x,lam,eps){
  n = dim(x)[1]
  p = dim(x)[2]
  xlam=sapply(1:p, FUN=function(i) bocotranvec(x[,i],lam[i],eps) ) 
  return(xlam)
}
#######################################################################
#                  function 7: find maximizer
#######################################################################
#'
#' @keywords internal
#' @noRd
#' 
argmax = function(x,y) x[ order(y)[length(x)] ]; 
#######################################################################
#       function 8: maximizing over one lambda in gauss seidel
#######################################################################
#'
#' @keywords internal
#' @noRd
#' 
gaussonestep=function(x,lam,eps,ilam,mlam, xmar, lam_bd, parallelize){
# ilam <- 1
# Stand the Margins here to avoid doing it mlam times
# xmar=standmar(x,eps);
  onelam=seq(from=-lam_bd,to=lam_bd,length=mlam)
  loglik <- sapply(1:mlam, FUN =function(k) { # k <- 1 # replaces/updates an element in lam with value from onelam;
    lam[ilam]=onelam[k];  
    return(likelihood(x, lam, eps, xmar, parallelize) )
    } )
  lamopt=argmax(onelam,loglik)
  lamout=lam 
  lamout[ilam]=lamopt
  return(lamout)
}
#######################################################################
#                function 9:  gauss seidel iterations
#######################################################################
#'
#' @keywords internal
#' @noRd
#' 
gauss=function(x,lam,eps,mlam,lam_bd = 2, n_iter=5, parallelize=F){
  # x; lam=rep(1,p); eps=.5; mlam=20; lam.bd = 2; n.iter=5; mc=FALSE
  # p = length(lam); lam=rep(1,p);
  xmar=standmar(x,eps);
  p = length(lam); 
  for(igauss in 1:n_iter) {
    for(i in 1:p){ # i <- 1
      lam1 = gaussonestep(x,lam,eps,i,mlam,xmar, lam_bd, parallelize) 
      lam = lam1 # Can't paralleize actually; this is iterative
    }
  }
  return(lam)
}

