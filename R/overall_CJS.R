#' Overall goodness-of-fit test for the Cormack-Jolly-Seber model
#'
#' This function performs the overall goodness-of-fit test for the Cormack-Jolly-Seber model.
#' It is obtained as the sum of the 4 components Test3.SR, Test3.SM, Test2.CT and Test2.CL.
#' @param X is a matrix of encounter histories
#' @param freq is a vector of the number of individuals with the corresponding encounter history
#' @param rounding is the level of rounding for outputs; default is 3
#' @return This function returns a data.frame with the value of the test statistic, the degrees of freedom and the p-value of the test.
#' @author Olivier Gimenez <olivier.gimenez@cefe.cnrs.fr>,Jean-Dominique Lebreton, Rémi Choquet, Roger Pradel
#' @keywords package
#' @export
#' @examples
#' # read in the classical dipper dataset
#' dipper = system.file("extdata", "ed.inp", package = "R2ucare")
#' dipper = read_inp(dipper,group.df=data.frame(sex=c('Male','Female')))
#'
#' # Get encounter histories, counts and groups:
#' dip.hist = dipper$encounter_histories
#' dip.freq = dipper$sample_size
#' dip.group = dipper$groups
#'
#' # split the dataset in males/females
#' mask = (dip.group == 'Female')
#' dip.fem.hist = dip.hist[mask,]
#' dip.fem.freq = dip.freq[mask]
#' mask = (dip.group == 'Male')
#' dip.mal.hist = dip.hist[mask,]
#' dip.mal.freq = dip.freq[mask]
#'
#' # for females
#' overall_CJS(dip.fem.hist, dip.fem.freq)

overall_CJS <- function(X,freq,rounding=3){

# compute each component
res_test3sr = test3sr(X, freq)
res_test3sm = test3sm(X, freq)
res_test2ct = test2ct(X, freq)
res_test2cl = test2cl(X, freq)

# compute overall test as the sum
stat = round(res_test3sr$test3sr[1] + res_test3sm$test3sm[1] + res_test2ct$test2ct[1] + res_test2cl$test2cl[1],rounding)
dof = res_test3sr$test3sr[2] + res_test3sm$test3sm[2] + res_test2ct$test2ct[2] + res_test2cl$test2cl[2]
pval = 1 - stats::pchisq(stat,dof)
pval = round(pval,rounding)

#cat('chi2, degree of freedom and p-value')
res = data.frame(chi2 = stat,degree_of_freedom = dof,p_value = pval)
row.names(res) = 'Gof test for CJS model:'
res
}
