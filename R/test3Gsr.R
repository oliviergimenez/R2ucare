#' Test3G.SR
#'
#' This function performs Test3G.SR
#' @param X is a matrix of encounter histories with K occasions
#' @param freq is a vector of the number of individuals with the corresponding encounter history
#' @param verbose controls the level of the details in the outputs; default is TRUE for all details
#' @param rounding is the level of rounding for outputs; default is 3
#' @return This function returns a list with first component the overall test and second component a data.frame with occasion, site, the value of the test statistic, degree of freedom, p-value and test performed (chi-square, Fisher or none).
#' @author Olivier Gimenez <olivier.gimenez@cefe.cnrs.fr>, Rémi Choquet, Roger Pradel
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
#' # perform Test3.GSR
#' test3Gsr(geese.hist,geese.freq)

test3Gsr <- function(X,freq,verbose=TRUE,rounding=3){

# various quantities to define
k = ncol(X)
res = group_data(X,freq)
his = res[,1:k]
eff = res[,k+1]
nh = nrow(his)
a = max(his)
kplusun = k + 1

# initialization
table_multi_3sr = data.frame(occasion = rep(NA,a*(k-2)),site = rep(NA,a*(k-2)),stat = rep(NA,a*(k-2)), df = rep(NA,a*(k-2)), p_val = rep(NA,a*(k-2)), test_perf = rep(FALSE,a*(k-2)))
where_in_table_3sr = 0
#nsitereel=sum(filtre);

# sort batches
for (i in 2:(k-1)){ # loop on date
    for (l in 1:a){ # loop on sites
    	where_in_table_3sr = where_in_table_3sr + 1
        #if filtre(l)
            fisheroupas=0

            masque = (his[,i]==l)
            batch = his[masque,] # select encounter histories containing l in column i
            if (length(batch)==0){ # if no release at date i on site l, no test
                table_multi_3sr[where_in_table_3sr,1] = i
                table_multi_3sr[where_in_table_3sr,2] = l
                table_multi_3sr[where_in_table_3sr,3] = 0
                table_multi_3sr[where_in_table_3sr,4] = 0
                table_multi_3sr[where_in_table_3sr,5] = 0
                table_multi_3sr[where_in_table_3sr,6] = 'None'
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
            # test3G SR, not pooled (is not tested)
            # ===================================================

            compo3GSSR = cbind(apply(table[,1:ncol(table)-1],1,sum),table[,ncol(table)])

            # ===================================================
            # test3G SR (sub-batch, releases)
            # ===================================================
            compo3GSR=rbind(compo3GSSR[1,],apply(compo3GSSR[2:(2+a-1),],2,sum))
            #if verbosity>=3
            #    strtable=[ strtable {strcat('component test3G SR (transients), occasion: ',num2str(i),' site: ',num2str(l)) }];
            #    strtable=[ strtable {'----------------- | Seen again - Never seen again' }];
            #    strtable=[ strtable {strcat('---- Newly Marked | ',num2str(compo3GSR(1,:))) }];
            #    strtable=[ strtable {strcat('Previously Marked | ',num2str(compo3GSR(2,:))) }];
            #end
            fisher3GSR = 0

            if (any(expval_table(compo3GSR)<2)){
                fisher3GSR = 1
            }

            #
            # calcul des tests pour  3G_SR
            #
            if (fisher3GSR == 1){
                fish = stats::fisher.test(compo3GSR)
                pvalfish = fish$p.value
                dffish = (nrow(compo3GSR)-1)*(ncol(compo3GSR)-1)
				stafish = stats::qchisq(1-pvalfish, dffish)
                table_multi_3sr[where_in_table_3sr,1] = i
                table_multi_3sr[where_in_table_3sr,2] = l
                table_multi_3sr[where_in_table_3sr,3] = stafish
                table_multi_3sr[where_in_table_3sr,4] = dffish
                table_multi_3sr[where_in_table_3sr,5] = pvalfish
                table_multi_3sr[where_in_table_3sr,6] = 'Fisher'
            } else {
            	   	old.warn <- options()$warn # to suppress the warning messages
            	   	options(warn = -1)
            	   	chi2 = stats::chisq.test(compo3GSR,correct=F)
            	   	options(warn = old.warn)
                pvalchi2 = chi2$p.value
                dfchi2 = chi2$parameter
				stachi2 = chi2$statistic
                table_multi_3sr[where_in_table_3sr,1] = i
                table_multi_3sr[where_in_table_3sr,2] = l
                table_multi_3sr[where_in_table_3sr,3] = stachi2
                table_multi_3sr[where_in_table_3sr,4] = dfchi2
                table_multi_3sr[where_in_table_3sr,5] = pvalchi2
                table_multi_3sr[where_in_table_3sr,6] = 'Chi-square'
            }

        }
    }
# compute overall test:
stat = sum(as.numeric(table_multi_3sr[,3]))
stat = round(stat,rounding)
dof = sum(as.numeric(table_multi_3sr[,4]))
pval = 1 - stats::pchisq(stat,dof)
pval = round(pval,rounding)
# if user specifies all outputs
if (verbose==TRUE) return(list(test3Gsr=c(stat=stat,df=dof,p_val=pval),details=table_multi_3sr))
# otherwise
if (verbose==FALSE) return(list(test3Gsr=c(stat=stat,df=dof,p_val=pval)))

}

