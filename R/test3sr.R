#' Test3.SR
#'
#' This function performs Test3.SR
#' @param X is a matrix of encounter histories with K occasions
#' @param freq is a vector of the number of individuals with the corresponding encounter history
#' @param verbose controls the level of the details in the outputs; default is TRUE for all details
#' @param rounding is the level of rounding for outputs; default is 3
#' @return This function returns a list with first component the overall test and second component a data.frame with 4 columns for components i (2:K-1) (in rows) of test3.sri: component, statistic of the test, p-value, signed test, test performed.
#' @author Olivier Gimenez <olivier.gimenez@cefe.cnrs.fr>, Jean-Dominique Lebreton, Rémi Choquet, Roger Pradel
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
#' # Test3SR for males
#' res.males = test3sr(dip.mal.hist, dip.mal.freq)
#' res.males

test3sr <- function(X,freq,verbose=TRUE,rounding=3){

n = dim(X)[1]
K = dim(X)[2]
result = data.frame(component = rep(NA,K-2),stat = rep(NA,K-2), p_val = rep(NA,K-2), signed_test = rep(NA,K-2), test_perf = rep(FALSE,K-2))
it3 = 2:(K-1) # occasions of capture 2 to K-1
before = rep(0,n) # nr of captures before 1
after = apply(X,1,sum) - X[,1] # nr of captures after 1
i = 0
for (j in it3){ # scan occasions of capture from 2 to K-1
   U = c(0,0,0,'None')
   i = i+1 # index from 1 to K-1 (rows of result)
   result[i,1] = j
   before = before + X[,j-1] # nr of captures before j
   after = after - X[,j] # nr of captures after j
   # if no capture and effj < 0 do not count (handles individuals not released)
   select = X[,j] & (!((after == 0) & (freq < 0)))
   # indices and number of individuals to keep in test3.sri
   ind = which(select)
   nj = length(ind)
   if (nj > 0){
   	newORold = matrix(0,nrow = nj,ncol = 2)
   	recORnev = newORold
    freqj = freq[ind]
    afreqj = abs(freqj)
   	newORold[,1] = (before[ind]>0)*afreqj # numbers seen before by RH
   	newORold[,2] = (before[ind]==0)*afreqj # numbers of new individuals by rh
    recORnev[,1] = (after[ind]!=0) # recaptured later
    recORnev[,2] = (after[ind]==0) # never recaptured
    # cont_tab 2x2 contingency table
    cont_tab = t(newORold) %*% recORnev
   	ML = apply(cont_tab,1,sum)
    MC = apply(cont_tab,2,sum)
    df = all(ML!=0) & all(MC!=0) # are row and column margins non-zero?
    if (df == TRUE) U = ind_test_22(cont_tab,2)
    } # if (nj > 0)
   #result[i,2] = df # df (1 ou 0)
   result[i,2] = U[1] # stat chi-square (also if Fisher performed)
   result[i,3] = U[2] # p-value of chi-square/Fisher
   result[i,4] = U[3] # signed test
   result[i,5] = U[4] # chi-square/Fisher
} # for (j in it3)
# compute overall test:
stat = sum(as.numeric(result[,2]))
stat = round(stat,rounding)
dof = sum(result$test_perf != 'None')
pval = 1 - stats::pchisq(stat,dof)
pval = round(pval,rounding)
tot_signed = round(sum(as.numeric(result$signed_test))/sqrt(length(result$signed_test)),rounding)
# if user specifies all outputs
if (verbose==TRUE) return(list(test3sr=c(stat=stat,df=dof,p_val=pval,sign_test = tot_signed),details=result))
# otherwise
if (verbose==FALSE) return(list(test3sr=c(stat=stat,df=dof,p_val=pval,sign_test = tot_signed)))

}
