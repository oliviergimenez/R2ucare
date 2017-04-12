#' Inverse generalized logit link
#'
#' This function computes the inverse (or reciprocal) of the generalized logit link function.
#' @param petitv vector of values to be transformed
#' @return ev vector of transformed values
#' @author Olivier Gimenez <olivier.gimenez@cefe.cnrs.fr>, Roger Pradel, Rémi Choquet
#' @keywords package
#' @export

inv_logit_gen <- function(petitv){
bound = rep(10,length(petitv))
petitv = pmin(petitv, bound)
petitv = pmax(petitv,-bound)
ev = exp(petitv)
dead = 1 / (1 + apply(ev,1,sum))
dead = as.matrix(dead)
ev = ev * diag(dead)
ev
}
