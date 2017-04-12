#' Test3G.WBWA
#'
#' This function performs Test3G.WBWA
#' @param X is a matrix of encounter histories with K occasions
#' @param freq is a vector of the number of individuals with the corresponding encounter history
#' @param verbose controls the level of the details in the outputs; default is TRUE for all details
#' @param rounding is the level of rounding for outputs; default is 3
#' @return This function returns a list with first component the overall test and second component a data.frame with occasion, site, the value of the test statistic, degree of freedom, p-value and test performed (chi-square, Fisher or none).
#' @author Olivier Gimenez <olivier.gimenez@@cefe.cnrs.fr>, Roger Pradel, Rémi Choquet
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
#' # perform Test.3GWBWA
#' test3Gwbwa(geese.hist,geese.freq)

test3Gwbwa <- function(X,freq,verbose=TRUE,rounding=3){

# various quantities to define
k = ncol(X)
res = group_data(X,freq)
his = res[,1:k]
eff = res[,k+1]
nh = nrow(his)
a = max(his)
kplusun = k + 1

# initialization
table_wbwa = data.frame(occasion = rep(NA,a*(k-2)),site = rep(NA,a*(k-2)),stat = rep(NA,a*(k-2)), df = rep(NA,a*(k-2)), p_val = rep(NA,a*(k-2)), test_perf = rep(FALSE,a*(k-2)))
where_in_table_wbwa = 0
#nsitereel=sum(filtre);

# sort batches
for (i in 2:(k-1)){ # loop on date
    for (l in 1:a){ # loop on sites
    	where_in_table_wbwa = where_in_table_wbwa + 1
        #if nsitereel>1 & filtre(l)==1
            fisheroupas = 0
            fisherWBWA = 0
            masque = (his[,i]==l)
            batch = his[masque,] # select encounter histories containing l in column i
            if (length(batch)==0){ # if no release at date i on site l, no test
                table_wbwa[where_in_table_wbwa,1] = i
                table_wbwa[where_in_table_wbwa,2] = l
                table_wbwa[where_in_table_wbwa,3] = 0
                table_wbwa[where_in_table_wbwa,4] = 0
                table_wbwa[where_in_table_wbwa,5] = 0
                table_wbwa[where_in_table_wbwa,6] = 'None'
                next
            }
            batcheff = eff[masque] # select counts corresponding to encounter histories with l in column i
            res = group_data_gen(batch,batcheff,(i+1):k) # sort according to columns i+1,...,k
            batchpost = res[,1:ncol(res)-1]
            batcheffpost = res[,ncol(res)]
            # look site on which previous obs occurred
            if (i!=2){
            		tt = t(apply(batchpost[,1:(i-1)],1,rev))
                eante = apply(tt!=0,1,which.max)
                eante = i - eante
            } else {
                eante = rep(1,nrow(batchpost))
            }
            # on cherche le site d'observation suivant
            if (i!=(k-1)){
                epost = apply(batchpost[,(i+1):k]!=0,1,which.max)
            } else {
                epost = rep(1,nrow(batchpost))
            }
            # build table of obs at date i on site l in rows
            # according to site of previous obs and in columns
            # according to site and date of next obs
            # PS: on first row, never seen before
            #     on last column, never seen again
            j = 0
            ind = (k-i)*a+1
            table = matrix(0,nrow=a+1,ncol=ind)
            cpt = matrix(0,nrow=a+1,ncol=1)
            while (j<nrow(batchpost)){ # go through encounter histories
                j=j+1
                date = epost[j]
                site = batchpost[j,i+date]
                if (site==0){
                    col = ind # never seen again
                } else{
                    col = (date-1)*a+site # seen again at date and site
                }
                if (j==1){
                    cold = col
                } else if (col!=cold){
                    table[,cold] = cpt
                    cpt = matrix(0,nrow=a+1,ncol=1)
                    cold = col
                }
                tempo = batchpost[j,eante[j]] + 1 # site of previous obs + 1
                if (site==0){
                    cpt[tempo] = cpt[tempo] + batcheffpost[j] * (batcheffpost[j]>0)
                } else {
                    cpt[tempo]=cpt[tempo]+abs(batcheffpost[j])
                }
            }
            table[,cold] = cpt #

            # ===================================================
            # test WBWA
            # ===================================================

            compoWBWA = table
            compoWBWA = compoWBWA[-1,]
            colWBWA = ncol(compoWBWA)
			compoWBWA = compoWBWA[,-colWBWA]

            for (j in 1:((colWBWA-1)/a-1)){
            	    if (((colWBWA-1)/a-1)<1) break
                compoWBWA[,1:a] = compoWBWA[,1:a] + compoWBWA[,(a+1):(2*a)]
                compoWBWA = compoWBWA[,-((a+1):(2*a))]
            }

            expvalWBWA = expval_table(compoWBWA) # table of expected values
            # reperer la plus faible valeur attendue de coord (nline,ncol)
            ind1 = apply(expvalWBWA,2,min)
            ind2 = apply(expvalWBWA,2,which.min)
            ncol = which.min(ind1)
            nline = ind2[ncol]
            vec_direction = c(ncol,nline)
            while (expvalWBWA[nline,ncol]<2){
                if (sum(compoWBWA[nline,])/ncol(compoWBWA) > sum(compoWBWA[,ncol])/nrow(compoWBWA)){
                    pooldim = 1 # on poole la colonne de la plus faible valeur attendue
                } else {
                    pooldim = 2
                }
                # on teste les dimensions de la table fraichement poolée
                # - soit aucune des dim ne vaut 2, alors on continue
                # - soit l'une des 2 dimensions vaut 2, on poole alors selon l'autre
                # - soit les 2 dimensions valent 2 et alors on s'arrete de pooler mais on rend test fisher
                # si ya encore des valeurs attendues < 2
                flag = which(dim(compoWBWA)==2)
                if (length(flag)==1){
                    pooldim = flag
                } else if (length(flag)==2){
                    fisherWBWA = 1
                    break # on sort du while et on garde la table 2 x 2 avec des valeurs attendues < 2
                }
                # je cherche la colonne ou ligne d'effectif le plus faible
                marge = apply(t(compoWBWA),pooldim,sum)
                marge[vec_direction[pooldim]] = max(marge) + 1 # on fixe la ligne ou la colonne qui contient
                # la valeur attendue la plus faible au max des obs.
                ind1 = min(marge)
                ind2 = which.min(marge) # ind2 = indice de la colonne ou de la ligne d'efectif la plus faible
                # il faut tester pooldim pour savoir si on poole les lignes ou les colonnes
                if (pooldim==1){
                    compoWBWA[,ncol] = compoWBWA[,ncol] + compoWBWA[,ind2]
                    compoWBWA = compoWBWA[,-ind2]
                } else {
                    compoWBWA[nline,] = compoWBWA[nline,] + compoWBWA[ind2,]
                    compoWBWA = compoWBWA[-ind2,]
                }
                expvalWBWA = expval_table(compoWBWA) # table des valeurs attendues
                # reperer la plus faible valeur attendue de coord (nline,ncol)
	            ind1 = apply(expvalWBWA,2,min)
    	        ind2 = apply(expvalWBWA,2,which.min)
        	    ncol = which.min(ind1)
            	nline = ind2[ncol]
                vec_direction = c(ncol,nline)
             }



#            if verbosity>=3
#                if verbosity>=3
#         strtable=[ strtable {strcat('component testWBWA, occasion: ',num2str(i),' site: ',num2str(l)) } ];
#         strtable=[ strtable {'-W-Before vs W-after- Next seen'} ];
#         for kk=1:size(compoWBWA,1)
#            strtable=[ strtable {strcat('Last seen_',num2str(kk),' | ',num2str(compoWBWA(kk,:)))}];
#         end
#                end
#            end



            # WBWA
            if (fisherWBWA == 1){
                fish = stats::fisher.test(compoWBWA)
                pvalfish = fish$p.value
                zeros_rows = (apply(compoWBWA,1,sum)==0)
                zeros_cols = (apply(compoWBWA,2,sum)==0)
                if (sum(!zeros_rows)+sum(!zeros_cols)==0){
                  dffish = 0
                } else {
                  dffish = (sum(!zeros_rows)-1)*(sum(!zeros_cols)-1)
                }
        				stafish = stats::qchisq(1-pvalfish, dffish)
                table_wbwa[where_in_table_wbwa,1] = i
                table_wbwa[where_in_table_wbwa,2] = l
                table_wbwa[where_in_table_wbwa,3] = stafish
                table_wbwa[where_in_table_wbwa,4] = dffish
                table_wbwa[where_in_table_wbwa,5] = pvalfish
                table_wbwa[where_in_table_wbwa,6] = 'Fisher'
            } else {
            	   	old.warn <- options()$warn # to suppress the warning messages
            	   	options(warn = -1)
            	   	chi2 = stats::chisq.test(compoWBWA,correct=F)
            	   	options(warn = old.warn)
                pvalchi2 = chi2$p.value
                dfchi2 = chi2$parameter
				stachi2 = chi2$statistic
                table_wbwa[where_in_table_wbwa,1] = i
                table_wbwa[where_in_table_wbwa,2] = l
                table_wbwa[where_in_table_wbwa,3] = stachi2
                table_wbwa[where_in_table_wbwa,4] = dfchi2
                table_wbwa[where_in_table_wbwa,5] = pvalchi2
                table_wbwa[where_in_table_wbwa,6] = 'Chi-square'
            }
        }
}
# compute overall test:
stat = sum(as.numeric(table_wbwa[,3]))
stat = round(stat,rounding)
dof = sum(as.numeric(table_wbwa[,4]))
pval = 1 - stats::pchisq(stat,dof)
pval = round(pval,rounding)
# if user specifies all outputs
if (verbose==TRUE) return(list(test3Gwbwa=c(stat=stat,df=dof,p_val=pval),details=table_wbwa))
# otherwise
if (verbose==FALSE) return(list(test3Gwbwa=c(stat=stat,df=dof,p_val=pval)))

}

