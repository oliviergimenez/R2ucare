#' m-array: table of first recaptures
#'
#' This function calculates the m-array, the number of released and never seen again individuals; deals with more than 1 group
#' @param X a matrix of encounter histories over K occasions
#' @param freq is a vector with the number of individuals having the corresponding encounter history
#' @return This function returns a list with R the number of released individuals (K-1 x g matrix), m the m-array (K-1 x K-1 x g array) with upper triangle filled only and never the number of individuals never recaptured (K-1 x g matrix).
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
#' # get female data
#' mask = (dip.group == 'Female')
#' dip.fem.hist = dip.hist[mask,]
#' dip.fem.freq = dip.freq[mask]
#'
#' # get number of released individuals (R), 
#' # the m-array (m) and 
#' # the number of individuals never seen again (never)
#' marray(dip.fem.hist,dip.fem.freq)

marray <- function(X,freq){

# INITIALIZATION
n = dim(X)[1]
k = dim(X)[2]
if (is.vector(freq)){
	g = 1
	freq = as.matrix(freq)
}
if (is.matrix(freq)) g = dim(freq)[2]
km1 = k-1
m = array(0,dim=c(km1,km1,g))
R = matrix(0,nrow=km1,ncol=g)
r = R
for (i in 1:n){ # loop on encounter histories
   h = X[i,] # pick up history in row i
   e = freq[i,] # get number of individuals with this ch in the g groups
   ae = abs(e) # absolute nvamue of numbers
   dates = which(h>0) # occasions of capture
   ncap=length(dates) # find number of captures
   if (ncap>0){ # more than 1 capture
      d1 = dates[1] # first capture
      dnext = d1 - 1
      for (j in 2:ncap){ # loop on dates of capture
         if (ncap<2) break # manage reverse loop
         dlast = dnext+1 # date of last capture
         dnext = dates[j] - 1 # date of next capture shifted left (to fill m)
         for (gg in 1:g){ # loop on groups
            R[dlast,gg] = R[dlast,gg] + ae[gg]
            m[dlast,dnext,gg] = m[dlast,dnext,gg] + ae[gg]
         }
      }
      # treat last release : nobody in m
      dlast = dnext + 1
      for (gg in 1:g){
         if ((e[gg]>0) & (dlast<k)) R[dlast,gg] = R[dlast,gg] + e[gg]
      }
   } # if ncap>0
} # for (i in 1:n)

for (d in 1:km1){
   for (gg in 1:g){
      v = m[d,,gg]
      r[d,gg] = r[d,gg] + sum(v)
   }
}
never = R-r
res = list(R=R,m=m,never=never)
res
}
