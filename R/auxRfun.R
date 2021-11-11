################################################################
#          Auxiliary R Functions
################################################################

################################################################
#          General Matrix Power Functions
################################################################
#' This is an internal function used for computing matrix powers
#' 
#'
#' @param a a square matrix 
#' @param alpha power of matrix to be taken 
#' @param lead number of leading eigenvalues to consider; default is all
#' @param ignore decimal place to ignore; default is 10^(-15)
#'
#' @return a matrix to the power alpha
#'
#' @keywords internal
#' 
#' @noRd
#' 
#' @export
#' 
matpower=function(a,alpha, lead=NULL, ignore=10^(-15)){
  a=(a+t(a))/2
  tmp=eigen(a,sym=T); p = length(tmp$values);
  eig_vals = tmp$values;
  
  # Incase smallest eigenvalues are negative by computational floating error 
  m0 = length(eig_vals[abs(eig_vals)>ignore])
  # tmp = evec[,1:m]%*%diag(eval[1:m]^power)%*%t(evec[,1:m])
  
  
  if (is.null(lead)){
    # ideally, m0=p here
    
    return(tmp$vectors[,1:m0]%*%
             diag((eig_vals[1:m0])^(alpha),m0,m0)%*%
             t(tmp$vectors[,1:m0]))
  } else {
    m = min(lead, m0)
    return(tmp$vectors[,1:m]%*%
             diag((eig_vals[1:m])^(alpha),m,m)%*%
             t(tmp$vectors[,1:m]))
  }
  
}


################################################################
#          distance between subspaces
################################################################
#' This is an internal function used for comparing the subspaces
#' spanned by two bases.
#' 
#' Taken from Li (2018).
#' 
#'
#' @param v1 matrix whose columns are a basis for the subspace A 
#' @param v2 matrix whose columns are a basis for the subspace B 
#'
#' @return a distance between subspace A and B
#'
#' @keywords internal
#' 
#' @noRd
#' 
mat_dist=function(v1,v2){
  v1=as.matrix(v1);v2=as.matrix(v2)
  if(dim(v1)[1]>1){
    p1 <- v1%*%matpower_cpp(t(v1)%*%v1,-1)%*%t(v1)
    p2 <- v2%*%matpower_cpp(t(v2)%*%v2,-1)%*%t(v2)
    if(dim(v1)[1]==1){
      p1=v1%*%t(v1)/c(t(v1)%*%v1)
      p2=v2%*%t(v2)/c(t(v2)%*%v2)}
    d <- sqrt(sum((p1-p2)*(p1-p2)))}
  return(d)
}