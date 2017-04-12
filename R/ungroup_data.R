#' Ungroup encounter capture-recapture data in individual histories
#'
#' This function splits encounter histories in as many individual histories as required.
#' @param X matrix of encounter capture-recapture histories
#' @param effX vector with numbers of individuals with that particular encounter history
#' @return matrix with ungrouped capture-recapture histories and counts in the last column (should be 1s)
#' @author Olivier Gimenez <olivier.gimenez@cefe.cnrs.fr>, Roger Pradel, Rémi Choquet
#' @keywords package
#' @export
#' @examples
#' # Generate fake capture-recapture dataset
#' X = matrix(round(runif(9)),nrow=3)
#' freq=c(4,3,-8)
#' cbind(X,freq)
#' ungroup_data(X,freq)

ungroup_data <- function(X,effX){

effX = as.matrix(effX)
g = ncol(effX)
bloc = matrix(0,nrow=1,ncol=g)
i = 1
ma = 0
s = ncol(X)

while (i <= nrow(X)){
	bloc = abs(effX[i,])
	if ((i < nrow(X))&&(sum(X[i+1,]==X[i,])==s)){
		bloc = bloc + abs(effX[i+1,])
		i = i + 1
	}
	ma = ma + max(bloc)
	i = i + 1
}

effY = matrix(nrow=ma,ncol=g)
Y = matrix(0,nrow=ma,ncol=ncol(X))
comp_init = 1
comp_fin = 1
i = 1

while (i <= nrow(X)){
	i0 = i
	for (gr in 1:g){
		comp = comp_init
		i = i0
		e = effX[i,gr]
		for (j in 1:abs(e)){
			Y[comp,] = X[i,]
			effY[comp,gr] = sign(e)
			comp=comp+1
		}
		if ((i < nrow(X))&&(sum(X[i+1,]==X[i,])==s)){
			i = i0 + 1
			e = effX[i,gr]
			for (j in 1:abs(e)){
				Y[comp,] = X[i,]
				effY[comp,gr] = sign(e)
				comp = comp + 1
			}
		}
		comp_fin = max(comp_fin,comp)
	}
	comp_init = comp_fin
	i = i + 1
} # while
return(cbind(Y,effY))
} # function
