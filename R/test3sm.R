#' Test3.SM
#'
#' This function performs Test3.SM
#' @param X is a matrix of encounter histories with K occasions
#' @param freq is a vector of the number of individuals with the corresponding encounter history
#' @param verbose controls the level of the details in the outputs; default is TRUE for all details
#' @param rounding is the level of rounding for outputs; default is 3
#' @return This function returns a list with first component the overall test and second component a data.frame with 5 columns for components i (2:K-1) (in rows) of test3.smi: component, degree of freedom, statistic of the test, p-value, test performed.
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
#' # for females
#' res.females = test3sm(dip.fem.hist, dip.fem.freq)
#' res.females

test3sm <- function(X,freq,verbose=TRUE,rounding=3){

n = dim(X)[1]
K = dim(X)[2]
result = data.frame(component = rep(NA,K-2),stat = rep(NA,K-2), df = rep(NA,K-2), p_val = rep(NA,K-2), test_perf = rep(FALSE,K-2))
it3 = 2:(K-1) # occasions of capture 2:K-1
before = rep(0,n) # nr of captures before 1
after = apply(X,1,sum)-X[,1] # nr of captures after 1
i = 0
for (j in it3){ # scan occasions of capture 2:K-1
   i = i+1
   result[i,1] = j
   before = before + X[,j-1] # nr of captures before j
   after = after - X[,j] # nr of captures after j
   select = X[,j] & (after>0) # individuals recaptured at least once after j
   ind = which(select) # rows of these individuals
   nj = length(ind)
   rest = (j+1):K # remaining occasions
   restmj = 1:(K-j)
   df = 0
   drapeau = 0
   U = rep(0,4)
   if (nj > 0){
      newORold = matrix(0,nrow=nj,ncol=2)
      if (length(rest)==1){ when = matrix(X[ind,rest],nrow=length(X[ind,rest]))}
      if (length(ind)==1){ when = matrix(X[ind,rest],ncol=length(X[ind,rest]))}
      if ((length(rest)!=1) & (length(ind)!=1)) {when = as.matrix(X[ind,rest])} # later recaptures of those captured in j
      deja = rep(0,nj)
      # restrict to next recaptured
      for (ii in 1:nj){
         for (jj in restmj){
            deja[ii] = deja[ii] + when[ii,jj] # number of reobservations
            if (deja[ii]>1) when[ii,jj] = 0 # if already seen again, neglect observation
          } # loop on jj (remaining occasions)
      } # loop on ii (selected individuals)
      # end calculation of next recapture
      effj = freq[ind] # works if eff is a columnn
      aeffj = abs(effj) # handles negative numbers (if seen again after j)
   	  newORold[,1] = (before[ind]>0)*aeffj # numbers of old by rh (nj x 2)
   	  newORold[,2] = (before[ind]==0)*aeffj # numbers of new by rh
      cont_tab = t(newORold) %*% when # (2 x nj by nj x K-j i.e. 2 by K-j
      cont_tab = pool2K(cont_tab,2) # pool if needed
      ML = apply(cont_tab,2,sum)
      MC = apply(cont_tab,1,sum)
      # calculate df
      if (any(as.logical(cont_tab))){ # non empty table
         df = (sum(ML>0)-1)*(sum(MC>0)-1) # takes account of empty rows and columns
      } else {
         df=0 # empty table
      }
      # end of df calculation
      if (df>0) {
      	U = ind_test_rc(cont_tab,2)
    } else {
      	U = c(0,0,0,'None')
      	}
   } # if (nj > 0)
   result[i,2] = U[1] # stat chi-square (also if Fisher performed)
   result[i,3] = U[3] # degree of freedom
   result[i,4] = U[2] # p-value of chi-square/Fisher
   result[i,5] = U[4] # chi-square/Fisher
} # for (j in it3)
# compute overall test:
stat = sum(as.numeric(result[,2]))
stat = round(stat,rounding)
dof = sum(as.numeric(result[,3]))
pval = 1 - stats::pchisq(stat,dof)
pval = round(pval,rounding)
# if user specifies all outputs
if (verbose==TRUE) return(list(test3sm=c(stat=stat,df=dof,p_val=pval),details=result))
# otherwise
if (verbose==FALSE) return(list(test3sm=c(stat=stat,df=dof,p_val=pval)))
}
