#' Test of independence for rxc contingency tables
#'
#' This function tests independence in rxc contingency tables
#' @param M is an r by c table of non-negative integers
#' @param threshold is a threshold for low expected numbers; default is 2
#' @param rounding is the level of rounding for outputs; default is 3
#' @return This function returns a vector with statistic of quadratic chi2 or inv chi2 corresponding to pvalue of Fisher test, p-value of quadratic chi2 test or Fisher test for low numbers, degree of freedom and test performed (Chi-square, Fisher or None).
#' @author Olivier Gimenez <olivier.gimenez@cefe.cnrs.fr>,Jean-Dominique Lebreton, Rémi Choquet, Roger Pradel
#' @keywords package
#' @export

ind_test_rc <- function(M,threshold=2,rounding=3){

# calculate margins
ML = apply(M,1,sum)
MC = apply(M,2,sum)
N = sum(MC)

# take care of empty rows and/or columns
rw = which(ML>0)
cl = which(MC>0)
if (((length(rw)==0) + (length(cl)==0))>0){ # no rows or no columns
   df = 0
} else {
   df = (length(rw)-1)*(length(cl)-1) # takes account of empty rows and columns
}

res = rep(0,4) # initialize table of results
res[3] = df
# calculation for a non-trivial table
if (df>0){
   M = M[rw,cl] # keeps non empty rows and columns
   TT = (1/N) * ML[rw] %*% t(MC[cl]) # calculate expected values on this subtable
   test_low = (sum(TT<threshold)>0)
   if (test_low) { # low numbers
      res[2] = as.numeric(stats::fisher.test(M)$p.value) # zero machine issue
      if (1-res[2] < 0.000000000000001) res[2] = 1 # zero machine issue
      res[1] = stats::qchisq(1-res[2],df)
	  res[1] = round(res[1],rounding)
      res[2] = round(res[2],rounding)
      res[4] = 'Fisher'
 }  else { # no low numbers
   	old.warn <- options()$warn # to suppress the warning messages
   	options(warn = -1)
   	res.tempo = stats::chisq.test(M,correct=F)
   options(warn = old.warn)
   res[1] = as.numeric(res.tempo$statistic)
   res[1] = round(res[1],rounding)
   res[2] = res.tempo$p.value
   res[2] = round(res[2],rounding)
   res[4] = 'Chi-square'}
} # if (df>0)
res
}

