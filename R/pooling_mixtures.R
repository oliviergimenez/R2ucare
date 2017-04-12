#' Pooling algorithm (multisite goodness-of-fit tests)
#'
#' This function pools rows and columns of a rxc bases and mixture table according to Pradel et al. (2003).
#' It provides the components of TestM in the multisite goodness-of-fit tests.
#' @param nk number of mixtures
#' @param nj number of bases
#' @param a number of sites/states
#' @param mixandbases matrix with mixtures and bases
#' @return This function returns a matrix with the pooled table.
#' @author Olivier Gimenez <olivier.gimenez@@cefe.cnrs.fr>, Rémi Choquet, Jean-Dominique Lebreton, Anne-Marie Reboulet, Roger Pradel
#' @keywords package
#' @references Pradel R., Wintrebert C.M.A. and Gimenez O. (2003). A proposal for a goodness-of-fit test to the Arnason-Schwarz multisite capture-recapture model. Biometrics 59: 43-53.
#' @export

pooling_mixtures <- function(nk,nj,a,mixandbases){

# pooling of mixture tables
# components i of TestM, denoted TestMi

expvalmixandbas = expval_table(mixandbases) # compute expected values
ind1 = apply(expvalmixandbas,2,min)
ind2 = apply(expvalmixandbas,2,which.min)
ncol = which.min(ind1) # coord col (ncol) with smallest expected value (ind3)
nline = ind2[ncol] # coord row (nline) with smallest expected value
vec_direction = c(ncol,nline)

# if smallest exp value is in bases, then pooling occurs only in columns (testmixorbas=1)
# if smallest exp value is in mixtures, no constraint (testmixorbas=0)
if (nline > (nk-nj)){
    testmixorbas = 1
} else {
    testmixorbas = 0
}

# pooling condition; while exp val < 2, proceed
while (expvalmixandbas[nline,ncol]<2){
    sumli = sum(mixandbases[nline,])/ncol(mixandbases) # row mean
    sumco = sum(mixandbases[,ncol])/nrow(mixandbases) # column mean
    if (testmixorbas==1){
        sumli = sumco + 1
    }
    if (sumli > sumco){
        pooldim = 1 # pool column
    } else {
        pooldim = 2
    }
    mindim = c(nj+1,nj+1)
    if (dim(mixandbases)[3-pooldim]==mindim[3-pooldim]){
        if (testmixorbas==1){
            expvalmixandbas[(nrow(expvalmixandbas)-a+1):nrow(expvalmixandbas),] = 2
            ind1 = apply(expvalmixandbas,2,min)
            ind2 = apply(expvalmixandbas,2,which.min)
            ncol = which.min(ind1)
            nline = ind2[ncol]
            vec_direction = c(ncol,nline)

            if (nline > nk-nj){
                # if mix, testmixorbas = 0 otherwise testmixorbas = 1
                testmixorbas = 1
            } else {
                testmixorbas = 0
            }
            next
        } else {
            if (dim(mixandbases)[pooldim]==mindim[pooldim]){
                break
            } else {
                pooldim = 3 - pooldim
            }
        }
    } # if

    # look for row/col with smallest count
    marge = apply(mixandbases,3-pooldim,sum)
    marge[vec_direction[pooldim]] = max(marge) + 1
    # set the row/col with smallest exp val to biggest obs
    if (pooldim==2){ # protect bases
        marge[length(marge)-nj+1:length(marge)] = max(marge) + 1
    }

    # ind2 = index of row/col with smallest count
    if (is.null(dim(marge))) marge = matrix(marge,ncol=1)
    ind1 = apply(marge,2,min)
    ind2 = apply(marge,2,which.min)
    # test pooldim to decide whether row or col should be pooled
    if (pooldim==1){
        mixandbases[,ncol] = mixandbases[,ncol] + mixandbases[,ind2]
        mixandbases = mixandbases[,-ind2]
    } else {
        mixandbases[nline,] = mixandbases[nline,] + mixandbases[ind2,]
        mixandbases = mixandbases[-ind2,]
    }

    expvalmixandbas = expval_table(mixandbases) #table of expected values
    # spot cell with smallest exp val (nline,ncol)
    ind1 = apply(expvalmixandbas,2,min)
    ind2 = apply(expvalmixandbas,2,which.min)
    ncol = which.min(ind1) # coord col (ncol) with smallest exp val (ind3)
    nline = ind2[ncol] # coord row (nline) with smallest exp val
    vec_direction = c(ncol,nline)

    nk = nrow(mixandbases)
    if (nline > nk-nj){
        # if mix, testmixorbas = 0 otherwise testmixorbas = 1
        testmixorbas = 1
    } else {
        testmixorbas = 0
    }


} # while
mixandbases
} # function
