#' Pooling algorithm
#'
#' This function pools columns of a 2xK contingency table (if needed, ie if low numbers present)
#' @param M is a 2 by K contingency table (or a K by 2 table)
#' @param low is a threshold for low expected numbers; default is 2 (if this argument is big enough, the table is pooled down to 2 x 2; if this argument is 0, the table is not pooled)
#' @return This function returns a matrix with the pooled contingency table.
#' @author Olivier Gimenez <olivier.gimenez@cefe.cnrs.fr>, Jean-Dominique Lebreton, Rémi Choquet, Roger Pradel
#' @keywords package
#' @export

pool2K <- function(M,low=2){

# get dimensions
rw = dim(M)[1]
cl = dim(M)[2]

# check dimensions, and transpose if needed
transp = 0
if (rw>2){
	M = t(M)
	transp = 1
	rw = dim(M)[1]
	cl = dim(M)[2]
}

# pool from right to left
# flip when OK from the right and not OK from the left
if (cl>2){ # more than two columns
  ML = apply(M,1,sum)
  MC = apply(M,2,sum)
  N = sum(MC)
  if (N>0){ # M is not full of zeros
     unsurN = 1/N
     TT = unsurN * ML %*% t(MC)
     flip = 0
     while (any(TT<low)){
        if (cl==2) break
        cm1 = cl - 1
        M[,cm1] = M[,cm1] + M[,cl]
        M = M[,1:cm1]
        cl = cm1
        MC = apply(M,2,sum)
        TT = unsurN * ML %*% t(MC)
        # restart from left if needed
        if (MC[1]<MC[cl]){
           M = M[,c(ncol(M):1),drop = FALSE]
           flip = 1-flip
        }
    }
     # transfer de-flipped or de-transposed in Mpooled
     if (flip==1){
        Mpooled = M[,c(ncol(M):1),drop = FALSE]
     } else {
        Mpooled = M
     }

     if (transp==1) Mpooled = t(Mpooled)
}  else { # M is full of zeros
     Mpooled = matrix(0,nrow=2,ncol=2)
  }
} else { # M is already 2 by 2
   Mpooled = M
} # if (cl>2)
Mpooled
}


