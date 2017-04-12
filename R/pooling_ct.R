#' Pooling algorithm (multisite goodness-of-fit tests)
#'
#' This function pools rows and columns of a rxc contingency table according to Pradel et al. (2003).
#' @param table is a rxc contingency table
#' @return This function returns a matrix with the pooled contingency table.
#' @author Olivier Gimenez <olivier.gimenez@cefe.cnrs.fr>, Jean-Dominique Lebreton, Rémi Choquet, Roger Pradel
#' @keywords package
#' @references Pradel R., Wintrebert C.M.A. and Gimenez O. (2003). A proposal for a goodness-of-fit test to the Arnason-Schwarz multisite capture-recapture model. Biometrics 59: 43-53.
#' @export

pooling_ct <- function(table){

expvaltable = expval_table(table) # table of expected values
# spot the smallest expected value of coordinate (nline,ncol)
ind1 = apply(expvaltable,2,min)
ind2 = apply(expvaltable,2,which.min)
ncol = which.min(ind1)
nline = ind2[ncol]
vec_direction = c(ncol,nline)

while (expvaltable[nline,ncol]<2){
    if (sum(table[nline,])/ncol(table) > sum(table[,ncol])/nrow(table)){
        pooldim = 1 # poole column of smallest expected value
    } else {
        pooldim = 2
    }

# test the dimensions of the table we just pooled
# if none equals 2, then continue
# if one dim equals 2, pool according to the other dim
# if both equal 2, stop pooling and render Fisher test if still exp val < 2
    flag = which(dim(table)==2)
    if (length(flag)==1){
         pooldim = flag
    } else if (length(flag)==2){
         fisheroupas = 1
         break
    }
    # look for row or column with smallest count
    marge = apply(t(table),pooldim,sum)
    marge[vec_direction[pooldim]] = max(marge) + 1 # set row/column with smallest expected value to biggest obs
    ind1 = min(marge)
    ind2 = which.min(marge) # ind2 = index of row/column with smallest count

    # test pooldim to decide whether rows or columns should be pooled
    if (pooldim==1){
        table[,ncol] = table[,ncol] + table[,ind2]
        table = table[,-ind2]
    } else {
        table[nline,] = table[nline,] + table[ind2,]
        table = table[-ind2,]
    }
    expvaltable = expval_table(table)
    # spot the smallest expected value of coordinate (nline,ncol)
    ind1 = apply(expvaltable,2,min)
    ind2 = apply(expvaltable,2,which.min)
    ncol = which.min(ind1)
    nline = ind2[ncol]
    vec_direction = c(ncol,nline)

}
table = as.matrix(table)
table
}

