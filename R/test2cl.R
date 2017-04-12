#' Test2.CL
#'
#' This function performs Test2.CL
#' @param X is a matrix of encounter histories with K occasions
#' @param freq is a vector of the number of individuals with the corresponding encounter history
#' @param verbose controls the level of the details in the outputs; default is TRUE for all details
#' @param rounding is the level of rounding for outputs; default is 3
#' @return This function returns a list with first component the overall test and second component a data.frame with 5 columns for components i (2:K-3) (in rows) of test2.cli following Pradel 1993 (in Lebreton and North, Birkhauser Verlag): component, degree of freedom, statistic of the test, p-value, test performed.
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
#' # for males
#' X = dip.mal.hist
#' freq = dip.mal.freq
#' res.males = test2cl(X,freq)
#' res.males

test2cl <- function(X,freq,verbose=TRUE,rounding=3){

# build m-array
m_array = marray(X,freq)
m = m_array$m[,,] # marray returns an array
km1 = dim(m)[1]
km4 = km1 - 3
result = data.frame(component=rep(NA, km4),dof=rep(NA, km4),stat=rep(NA, km4),p_val=rep(NA, km4),test_perf=rep(TRUE, km4))
for (indice in 1:km4){
     i = indice + 1
     result[indice,1] = i
     rw = 1:(i-1)
     cl = (i+1):km1
     cont_tab = matrix(0,nrow=2,ncol=length(cl))
     if (i>2){
        cont_tab[1,] = apply(m[rw,cl],2,sum)
     } else {
        cont_tab[1,] = m[rw,cl]
     }
     cont_tab[2,] = m[i,cl]
     cont_tab = pool2K(cont_tab,2) # pool if needed
     MC = apply(cont_tab,2,sum)
     ML = apply(cont_tab,1,sum)
# calculate df
if (any(as.logical(cont_tab))){ # non empty table
	df=(sum(ML>0)-1)*(sum(MC>0)-1) # account for empty rows and columns
} else {
	df=0 # empty table
}
# end of df calculation
if (df>0){
	U = ind_test_rc(cont_tab) # test of ind
} else {
	U = c(0,0,0,'None')
}

     result[indice,2] = U[3] # df = 1 or 0
     result[indice,3] = U[1]
     result[indice,4] = U[2]
     result[indice,5] = U[4]
}
# compute overall test:
stat = sum(as.numeric(result[,3]))
stat = round(stat,rounding)
dof = sum(as.numeric(result[,2]))
pval = 1 - stats::pchisq(stat,dof)
pval = round(pval,rounding)
# if user specifies all outputs
if (verbose==TRUE) return(list(test2cl=c(stat=stat,df=dof,p_val=pval),details=result))
# otherwise
if (verbose==FALSE) return(list(test2cl=c(stat=stat,df=dof,p_val=pval)))

}
