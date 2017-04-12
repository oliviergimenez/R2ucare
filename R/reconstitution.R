#' Reformat outputs of multinomial mixture distributions parameters
#'
#' This function reformat the outputs of multinomial mixture distributions parameters.
#' @param x vector with cell probabilities estimates for mixtures and bases, along with mixture probilities
#' @param s number of bases
#' @param n number of cell probabilities
#' @param nbmel number of mixtures
#' @return This function returns a list of maximum likelihood estimates for the cells of a mixture distribution with:
#' @return P matrix of cell probabilities estimates for mixtures
#' @return PI matrix of mixture probabilities
#' @return GAM matrix of cell probabilities estimates for bases
#' @author Olivier Gimenez <olivier.gimenez@cefe.cnrs.fr>, Roger Pradel, Rémi Choquet
#' @keywords package
#' @export

reconstitution <- function(x,s,n,nbmel){

if (s>1){
    # parameters linked to Gamma
    y=x[1:(nbmel*(s-1))]
    # reconstruction of Gamma
    dim(y) = c(nbmel,s-1)
    # inverse generalised logit link
    for (i in 1:nbmel){
        y[i,] = t(inv_logit_gen(t(y[i,])))
    }
    compy = 1 - apply(y,1,sum)
    Gam = cbind(compy,y)
} else {
    Gam = rep(1,nbmel)
}

# parameters linked to Pi
z = x[(nbmel*(s-1)+1):length(x)]

# reconstruction of Pi
dim(z) = c(s,n-1)
# inverse generalised logit link
for (i in 1:s){
    z[i,] = t(inv_logit_gen(t(z[i,])))
}
compz = 1-apply(z,1,sum)
Pi = cbind(compz,z)

# reconstitution of parameters Pi and P
P = Gam%*%Pi

res = list(P=P,Gam=Gam,Pi=Pi)
res
}
