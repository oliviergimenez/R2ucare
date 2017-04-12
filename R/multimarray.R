#' Multistate m-array
#'
#' This function calculates the m-array for multistate capture-recapture data, the number of released and never seen again individuals.
#' @param X a matrix of encounter histories over K occasions
#' @param freq is a vector with the number of individuals having the corresponding encounter history
#' @return This function returns a matrix in which R the number of released individuals is in the first column, the number of individuals never recaptured (K-1) is in the last column and m the m-array (K-1 x K-1) with upper triangle filled only is in sandwich between these two vectors.
#' @author Olivier Gimenez <olivier.gimenez@cefe.cnrs.fr>,Jean-Dominique Lebreton, Rémi Choquet, Roger Pradel
#' @keywords package
#' @export
#' @examples
#' # Read in Geese dataset:
#' geese = system.file("extdata", "geese.inp", package = "R2ucare")
#' geese = read_inp(geese)
#'
#' # Get encounter histories and number of individuals with corresponding histories
#' geese.hist = geese$encounter_histories
#' geese.freq = geese$sample_size
#'
#' # build m-array
#' multimarray(geese.hist, geese.freq)

multimarray <- function(X,freq){

  n = dim(X)[1] # nb of encounter histories
  s = max(X) # nb of sites
  k = dim(X)[2]
  km1 = k - 1 # nb of recapture occasions
  effx = matrix(freq)
  g = ncol(effx)

  #get.first <- function(x) min(which(x!=0))
  #e = apply(X, 1, get.first) # beginning of encounter histories in e
  #f = e
  #for (i in 1:nh) f[i] = chmatx[e[i],i] # site of initial release

  R = matrix(0,nrow=s,ncol=km1)
  m = array(0,dim=c(s,s,km1,km1))
  rprime = R
  # loop on encounter histories
  for (i in 1:n){
    h = X[i,] # pick up history in row i
    e = effx[i,] # get number of individuals with this ch in the g groups
    ae = abs(e) # absolute value of numbers
    dates = which(h>0) # occasions of capture
    ncap = length(dates) # find number of captures
    if (ncap>0){ # more than 1 capture
      if (ncap>1){
        dlast = dates[1:(ncap-1)]
        dnext = dates[2:ncap]
        sd = unlist(h[dlast])
        sa = unlist(h[dnext])
        dnext = dnext-1
        for (j in 1:(ncap-1)){
          R[sd[j],dlast[j]] = R[sd[j],dlast[j]] + ae
          m[sd[j],sa[j],dlast[j],dnext[j]] = m[sd[j],sa[j],dlast[j],dnext[j]] + ae
        }
      } # if (ncap>1)
      dl = dates[ncap]
      sdl = unlist(h[dl])
      if (dl<k) R[sdl,dl] = R[sdl,dl] + 0.5 * (e + ae)
    } # if (ncap>0)
  } # for (i in 1:n)

  v = m
  for (ns in 1:s){
    for (d in 1:km1){
      v = m[ns,,d,] # sum over strata of arrival and occasions
      rprime[ns,d] = rprime[ns,d] + sum(sum(v))
    }
  }
  never = R - rprime

  # format output
  released = c(R)
  never = c(never)

  marray = NULL
  for (i in 1:km1){
    for (j in 1:km1){
      marray = cbind(marray,m[,,i,j])
    }
  }
  deb = seq(1,ncol(marray),by=s*km1)
  fin = seq(s*km1,ncol(marray),by=s*km1)
  marray2 = NULL
  for (j in 1:km1){
    marray2 = rbind(marray2,marray[,deb[j]:fin[j]])
  }

  marray = cbind(released,marray2,never)
  marray

}
