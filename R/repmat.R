#' Replicate and tile a matrix
#'
#' This function creates a large matrix consisting of an m-by-n tiling of copies of X.
#' The dimensions of the returned matrix are nrow(X)*m x ncol(X)*n.
#' This is the equivalent of the repmat MATLAB function.
#' @param X matrix to be replicated
#' @param m row dimension of replication
#' @param n column dimension of replication
#' @return A replicated matrix of X with dimensions nrow(X)*m x ncol(X)*n.
#' @author Olivier Gimenez <olivier.gimenez@cefe.cnrs.fr>
#' @keywords package
#' @export

repmat = function(X,m,n){
mx = dim(X)[1]
nx = dim(X)[2]
res = matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
res
}
