#' Expected values in a contingency table
#'
#' This function calculates expected values for a rxc contingency table.
#' @param M a matrix of observed probabilities
#' @return A matrix of expected values.
#' @author Olivier Gimenez <olivier.gimenez@cefe.cnrs.fr>, Roger Pradel, Rémi Choquet
#' @keywords package
#' @export

expval_table <- function(M) {

  nblines = nrow(M)
  nbcolumns = ncol(M)
  ML = as.matrix(apply(M, 1, sum))
  MC = t(as.matrix(apply(M, 2, sum)))
  tot = sum(ML)
  if (tot != 0) {
    RES = repmat(ML, 1, nbcolumns) * repmat(MC, nblines, 1) / tot
  } else {
    RES = M
  }
  RES
}
