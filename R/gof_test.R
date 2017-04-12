#' Goodness-of-fit test for contingency tables
#'
#' This function carries out goodness-of-fit tests for contingency tables from the power-divergence family.
#' @param lambda parameter defining the statistic to be used: lambda = -0.5 is for the Freeman-Tuckey statistic, lambda = 0 for the G2 statistic, lambda = 2/3 for the Cressie-Read statistic and lambda = 1 for the classical Chi-square statistic
#' @param observes vector of observed probabilities
#' @param theoriques vector of theoretical/expected probabilities
#' @return This function returns the value of the goodness-of-fit statistic.
#' @author Olivier Gimenez <olivier.gimenez@cefe.cnrs.fr>, Roger Pradel, Rémi Choquet
#' @keywords package
#' @export

gof_test <- function(lambda, observes, theoriques) {

  cache1 = (observes == 0)
  if (sum(cache1) > 0) {
    observes = observes[!cache1]
    theoriques = theoriques[!cache1]
  }

  cache2 = (theoriques < 0.01)
  if (sum(cache2) > 0) {
    observes = observes[!cache2]
    theoriques = theoriques[!cache2]
  }

  if (lambda == 0) {
    # the G2 statistic is obtained as the limit of the statistic when lambda tends to 0
    stat = 2 * sum(observes * log((observes / theoriques)))
  } else {
    stat = 2 / (lambda * (lambda + 1)) * sum(observes * (((
      observes / theoriques
    ) ^ lambda) - 1))
  }
  stat
}
