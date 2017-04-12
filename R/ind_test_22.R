#' Test of independence for 2x2 contingency tables
#'
#' This function tests independence in 2x2 contingency tables
#' @param M is a 2x2 contingency table
#' @param threshold is a threshold for low expected numbers; default is 2
#' @param rounding is the level of rounding for outputs; default is 3
#' @return This function returns a vector with statistic of quadratic chi2 or inv chi2 corresponding to pvalue of Fisher test, p-value of quadratic chi2 test or Fisher test for low numbers, signed test and test performed (Chi-square, Fisher or None).
#' @author Olivier Gimenez <olivier.gimenez@cefe.cnrs.fr>,Jean-Dominique Lebreton, Rémi Choquet, Roger Pradel
#' @keywords package
#' @export

ind_test_22 <- function(M,threshold=2,rounding=3){

# calculate margins, total, expected table
MC = apply(M,2,sum)
ML = apply(M,1,sum)
N = sum(ML)

# take care of empty rows and/or columns
rw = which(ML>0)
cl = which(MC>0)
if (((length(rw)==0)+(length(cl)==0)) > 0) {
   df = 0
} else {
   df = (length(rw)-1)*(length(cl)-1)
}

res = rep(0,4)
if (df>0){ # perform test
   M = M[rw,cl] # keeps non empty rows and columns BEHAVIOR OF chisq/Fisher.test TO BE CHECKED IN R
   TT = (1/N) * ML[rw] %*% t(MC[cl]) # calculate expected values on this subtable
   D = M - TT
   test_low = (sum(TT<threshold)>0)
   res[3] = test_low
   if (test_low) {
   res[2] = stats::fisher.test(M)$p.value
   res[1] = stats::qchisq(1-res[2],1)
   res[1] = round(res[1],rounding)
   res[2] = round(res[2],rounding)
   res[3] = round(sign(D[1,1]) * sqrt(res[1]),rounding)
   res[4] = 'Fisher'}
   else {
   	old.warn <- options()$warn # to suppress the warning messages
   	options(warn = -1)
   	res.tempo = stats::chisq.test(M,correct=F)
   options(warn = old.warn)
   res[1] = res.tempo$statistic
   res[1] = round(res[1],rounding)
   res[2] = res.tempo$p.value
   res[2] = round(res[2],rounding)
   res[3] = round(sign(D[1,1]) * sqrt(res[1]),rounding)
   res[4] = 'Chi-square'}
} # if df>0
res
}
