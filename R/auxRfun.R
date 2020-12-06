################################################################
#          Auxiliary R Functions
################################################################


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