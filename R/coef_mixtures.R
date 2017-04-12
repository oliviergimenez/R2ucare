#' Estimation of multinomial mixture distributions parameters
#'
#' This function performs maximum likelihood inference for multinomial mixture distributions.
#' @param Mp a matrix of mixtures (a row matrix if a vector)
#' @param Np a matrix of bases (a row matrix if a vector)
#' @return This function returns a list of maximum likelihood estimates for the cells of a mixture distribution:
#' @return P matrix of cell probabilities estimates for mixtures
#' @return PI matrix of mixture probabilities
#' @return GAM matrix of cell probabilities estimates for bases
#' @author Olivier Gimenez <olivier.gimenez@cefe.cnrs.fr>, Roger Pradel, Rémi Choquet
#' @keywords package
#' @references Yantis, S., Meyer, D. E., and Smith, J. E. K. (1991). Analyses of multinomial mixture distributions: New tests for stochastic models of cognition and action. Psychological Bulletin 110, 350–374.
#' @export

coef_mixtures <- function(Mp,Np){

# various quantities to be defined
M = Mp
N = Np
s = nrow(N) # s = nb of bases
n = ncol(N)
nbmel = nrow(M) # nb of mixtures
dim(M) = c(1,nbmel*n)
dim(N) = c(1,s*n)


# run 15 times the optimization procedure to try and handle with local minima
# initial values
x = matrix(stats::runif(nbmel*(s-1)+s*(n-1)),nrow=nbmel*(s-1)+s*(n-1),ncol=1)
# Minimization
tmpmin = stats::optim(x,deviance_mixture,NULL,hessian=FALSE,M,N,s,n,nbmel,method="BFGS",control=list(trace=0, reltol=.0000001,abstol=.000001))
for (i in 1:14){
  x2 = matrix(stats::runif(nbmel*(s-1)+s*(n-1)),nrow=nbmel*(s-1)+s*(n-1),ncol=1)
  # Minimization
  tmpmin2 = stats::optim(x2,deviance_mixture,NULL,hessian=FALSE,M,N,s,n,nbmel,method="BFGS",control=list(trace=0, reltol=.0000001,abstol=.000001))
  #tmpmin2$value
  if (tmpmin2$value < tmpmin$value) {tmpmin = tmpmin2}
  }
x <- tmpmin$par

# reconstruction of Gam,Pi,P
rec_param = reconstitution(x,s,n,nbmel)
P = rec_param$P
Gam = rec_param$Gam
Pi = rec_param$Pi

# outputs
res = list(P=P, PI=Pi, GAM=Gam)
res

}
