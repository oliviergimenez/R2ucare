#' Test2.CT
#'
#' This function performs Test2.CT
#' @param m is a triangular m-array (in a K-1 by K-1 array)
#' @param verbose controls the level of the details in the outputs; default is TRUE for all details
#' @param rounding is the level of rounding for outputs; default is 3
#' @return This function returns a list with first component the overall test and second component a data.frame with 5 columns for components i (2:K-2) (in rows) of test2.Cti: component, degree of freedom, statistic of the test, p-value, test performed.
#' @author Olivier Gimenez <olivier.gimenez@@cefe.cnrs.fr>, Rémi Choquet, Jean-Dominique Lebreton, Anne-Marie Reboulet, Roger Pradel
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
#' X = dip.fem.hist
#' freq = dip.fem.freq
#' res.fem = marray(X,freq)
#' m = res.fem$m[,,] # marray returns an array
#' res.females = test2ct(m)
#' res.females

test2ct <- function(m,verbose=TRUE,rounding=3){

km1 = dim(m)[1]
km3 = km1 - 2
cont_tab = matrix(0,nrow=2,ncol=2)
result = data.frame(component=rep(NA,km3),dof=rep(NA,km3),stat=rep(NA,km3),p_val=rep(NA,km3),test_perf=rep(TRUE,km3))
for (indice in 1:km3){
     i = indice + 1
     result[indice,1] = i
     rw = 1:(i-1)
     cl = (i+1):km1
     # make test2ct table
     HG = m[rw,i]
     HD = m[rw,cl]
     BG = m[i,i]
     BD = m[i,cl]
     cont_tab[1,1] = sum(HG)
     cont_tab[1,2] = sum(sum(HD))
     cont_tab[2,1] = BG
     cont_tab[2,2] = sum(BD)
     # margins 
     ML = apply(cont_tab,1,sum)
     MC = apply(cont_tab,2,sum)
     df = all(ML!=0) & all(MC!=0)
     # statistics
     if (df==TRUE){ 
     	U = ind_test_22(cont_tab,2)
    }else {
    		U = c(0,0,'None')
    }
     result[indice,2] = as.numeric(df) # df = 1 or 0
     result[indice,3] = U[1]
     result[indice,4] = U[2]     
     result[indice,5] = U[3]     
} # for (indice in 1:km3)

# compute overall test:
stat = sum(as.numeric(result[,3]))
stat = round(stat,rounding)
dof = sum(result[,2])
pval = 1 - pchisq(stat,dof)
pval = round(pval,rounding)
# if user specifies all outputs
if (verbose==TRUE) return(list(test2ct=c(stat=stat,df=dof,p_val=pval),details=result))
# otherwise
if (verbose==FALSE) return(list(test2ct=c(stat=stat,df=dof,p_val=pval)))

}