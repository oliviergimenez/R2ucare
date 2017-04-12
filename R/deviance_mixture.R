#' Deviance of multinomial mixture distributions
#'
#' This function calculates the deviance of multinomial mixture distributions.
#' @param x value to which the deviance is to be evaluated
#' @param M a vector of mixtures (see coef_mixtures.R)
#' @param N a vector of bases (see coef_mixtures.R)
#' @param s number of bases
#' @param n number of cell probabilities
#' @param nbmel number of mixtures
#' @return This function returns the value of the deviance for mixture distributions.
#' @author Olivier Gimenez <olivier.gimenez@cefe.cnrs.fr>, Roger Pradel, Rémi Choquet
#' @keywords package
#' @references Yantis, S., Meyer, D. E., and Smith, J. E. K. (1991). Analyses of multinomial mixture distributions: New tests for stochastic models of cognition and action. Psychological Bulletin 110, 350–374.
#' @export

deviance_mixture <- function(x,M,N,s,n,nbmel){

res = reconstitution(x,s,n,nbmel)
P = res$P
Gam = res$Gam
# if Gam a vector, make it a matrix
if (is.null(dim(Gam))) Gam = matrix(Gam,nrow=1)
Pi = res$Pi

# transform in vectors
dim(P) = c(nrow(P)*ncol(P),1)
dim(Gam) = c(nrow(Gam)*ncol(Gam),1)
dim(Pi) = c(nrow(Pi)*ncol(Pi),1)
P = as.vector(P)
Gam = as.vector(Gam)
Pi = as.vector(Pi)
N = as.vector(N)
M = as.vector(M)

# compute deviance and return the value
f = sum(M * log(P))
f = f + sum(N * log(Pi))
f = -2 * f
f

}
