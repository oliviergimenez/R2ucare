#' Test3G.SM
#'
#' This function performs Test3G.SM
#' @param X is a matrix of encounter histories with K occasions
#' @param freq is a vector of the number of individuals with the corresponding encounter history
#' @param verbose controls the level of the details in the outputs; default is TRUE for all details
#' @param rounding is the level of rounding for outputs; default is 3
#' @return This function returns a list with first component the overall test and second component a data.frame with occasion, site, the value of the test statistic, degree of freedom, p-value and test performed (chi-square, Fisher or none).
#' @author Olivier Gimenez <olivier.gimenez@cefe.cnrs.fr>, Roger Pradel, Rémi Choquet
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
#' # perform Test.3.GSm
#' test3Gsm(geese.hist,geese.freq)

test3Gsm <- function(X,freq,verbose=TRUE,rounding=3){

# various quantities to define
k = ncol(X)
res = group_data(X,freq)
his = res[,1:k]
eff = res[,k+1]
nh = nrow(his)
a = max(his)
ns = a
kplusun = k + 1

# initialization
table_multi_3sm = data.frame(occasion = rep(NA,a*(k-2)),site = rep(NA,a*(k-2)),stat = rep(NA,a*(k-2)), df = rep(NA,a*(k-2)), p_val = rep(NA,a*(k-2)), test_perf = rep(FALSE,a*(k-2)))
where_in_table_3sm = 0
#nsitereel=sum(filtre)
stattotal = NULL

# sort batches
for (i in 2:(k-1)){ # loop on date
    for (l in 1:a){ # loop on sites
    	where_in_table_3sm = where_in_table_3sm + 1
        #if filtre(l)
            fisheroupas=0

            masque = (his[,i]==l)
            batch = his[masque,] # select encounter histories containing l in column i
            if (length(batch)==0){ # if no release at date i on site l, no test
                table_multi_3sm[where_in_table_3sm,1] = i
                table_multi_3sm[where_in_table_3sm,2] = l
                table_multi_3sm[where_in_table_3sm,3] = 0
                table_multi_3sm[where_in_table_3sm,4] = 0
                table_multi_3sm[where_in_table_3sm,5] = 0
                table_multi_3sm[where_in_table_3sm,6] = 'None'
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
            # Build table for test3G SM with pooling
            # No distinction between bases and mixtures, pooling is performed
            # on rows/columns indistinctively
            # ===================================================

            table1 = table[1,1:(ncol(table)-1)]
            table2 = table[1,ncol(table)]
            tablerevu = table[2:(1+ns),1:(ncol(table)-1)]
            table4 = table[2:(1+ns),ncol(table)]

            #   affichage des tables au fur et a mesure de leur creation
            table = rbind(table1,apply(tablerevu,2,sum))
            if (is.null(dim(table))) {table=as.matrix(table)} # if table is a vector
            if ((nrow(table)*ncol(table))>=4){
                table = pooling_ct(table)
                #if verbosity>=3
                #    strtable=[ strtable {'============================================================================='} ];
                #    strtable=[ strtable {strcat('component test3G SM, occasion: ',num2str(i),' site: ',num2str(l))} ];
                #    strtable=[ strtable {strcat('New from everywhere --------- |',num2str(table(1,:))) }];
                #    strtable=[ strtable {strcat('Has been seen early somewhere |',num2str(table(2,:))) }];
                #    %             table
                #end

            if (any(expval_table(table)<2)){
                fish = stats::fisher.test(table)
                pvalfish = fish$p.value
                zeros_rows = (apply(table,1,sum)==0)
                zeros_cols = (apply(table,2,sum)==0)
                if (sum(!zeros_rows)+sum(!zeros_cols)==0){
                  dffish = 0
                } else {
                  dffish = (sum(!zeros_rows)-1)*(sum(!zeros_cols)-1)
                }
				        stafish = stats::qchisq(1-pvalfish, dffish)
                table_multi_3sm[where_in_table_3sm,1] = i
                table_multi_3sm[where_in_table_3sm,2] = l
                table_multi_3sm[where_in_table_3sm,3] = stafish
                table_multi_3sm[where_in_table_3sm,4] = dffish
                table_multi_3sm[where_in_table_3sm,5] = pvalfish
                table_multi_3sm[where_in_table_3sm,6] = 'Fisher'
            } else {
            	 old.warn <- options()$warn # to suppress the warning messages
            	 options(warn = -1)
            	 chi2 = stats::chisq.test(table,correct=F)
            	 options(warn = old.warn)
               pvalchi2 = chi2$p.value
               dfchi2 = chi2$parameter
				       stachi2 = chi2$statistic
               table_multi_3sm[where_in_table_3sm,1] = i
               table_multi_3sm[where_in_table_3sm,2] = l
               table_multi_3sm[where_in_table_3sm,3] = stachi2
               table_multi_3sm[where_in_table_3sm,4] = dfchi2
               table_multi_3sm[where_in_table_3sm,5] = pvalchi2
               table_multi_3sm[where_in_table_3sm,6] = 'Chi-square'
            }
                #if verbosity>=3
                #    strtable=[ strtable {strcat('Associated test of the last table :',num2str(stat(1:3))) } ];
                #end
            }

            if (ns>1){ # not defined for single site
                table = cbind(apply(tablerevu,1,sum),table4)
                 if ((nrow(table)*ncol(table))>=4){
                    table = pooling_ct(table)
             #       if verbosity>=3
             #           strtable=[ strtable {'------------------------------------------------------------'} ];
             #           strtable=[ strtable {'------------------- Are seen again (and) not seen again'} ];
             #           for kk=1:size(table,1)
             #               strtable=[ strtable {strcat('previously seen at site_',num2str(kk),' | ',num2str(table(kk,:)))} ];
             #           end
             #           %             table
             #       end

            if (any(expval_table(table)<2)){
                fish = stats::fisher.test(table)
                pvalfish = fish$p.value
                zeros_rows = (apply(table,1,sum)==0)
                zeros_cols = (apply(table,2,sum)==0)
                if (sum(!zeros_rows)+sum(!zeros_cols)==0){
                  dffish = 0
                } else {
                  dffish = (sum(!zeros_rows)-1)*(sum(!zeros_cols)-1)
                }
                stafish = stats::qchisq(1-pvalfish, dffish)
                table_multi_3sm[where_in_table_3sm,1] = i
                table_multi_3sm[where_in_table_3sm,2] = l
                table_multi_3sm[where_in_table_3sm,3] = stafish + table_multi_3sm[where_in_table_3sm,3]
                table_multi_3sm[where_in_table_3sm,4] = dffish + table_multi_3sm[where_in_table_3sm,4]
                table_multi_3sm[where_in_table_3sm,5] = pvalfish
                table_multi_3sm[where_in_table_3sm,6] = 'Fisher'
                #           if verbosity>=3
             #               strtable=[ strtable {strcat('Associated test of the last table :',num2str([stafish pvalfish dffish]))} ];
             #           end
            } else {
            	   	old.warn <- options()$warn # to suppress the warning messages
            	   	options(warn = -1)
            	   	chi2 = stats::chisq.test(table,correct=F)
            	   	options(warn = old.warn)
                  pvalchi2 = chi2$p.value
                  dfchi2 = chi2$parameter
				        stachi2 = chi2$statistic
                table_multi_3sm[where_in_table_3sm,1] = i
                table_multi_3sm[where_in_table_3sm,2] = l
                table_multi_3sm[where_in_table_3sm,3] = stachi2 + table_multi_3sm[where_in_table_3sm,3]
                table_multi_3sm[where_in_table_3sm,4] = dfchi2 + table_multi_3sm[where_in_table_3sm,4]
                table_multi_3sm[where_in_table_3sm,5] = pvalchi2 + table_multi_3sm[where_in_table_3sm,5]
                table_multi_3sm[where_in_table_3sm,6] = 'Chi-square'
             #           if verbosity>=3
             #               strtable=[ strtable {strcat('Associated test of the last table :',num2str(chi2(table)))} ];
             #           end
            }
         }

                #if verbosity>=3
                #    strtable=[ strtable {'------------------------------------------------------------'} ];
                #    strtable=[ strtable {strcat('For this table, only next reencounters on site ',num2str(l) ,' are kept')} ];
                #end

                for (j in 1:ns){
                    table = tablerevu[,seq(j,ncol(tablerevu),by=ns)]
                    if (is.null(dim(table))) {table=as.matrix(table)} # if table is a vector
                   if ((nrow(table)*ncol(table))>=4){
                        table = pooling_ct(table)
               #         if verbosity>=3
               #             strtable=[ strtable {strcat('------------------- date of next reencounter on site ',num2str(j) ,' (after pooling)')} ];
               #             for kk=1:size(table,1)
               #                 strtable=[ strtable {strcat('previously seen at site_',num2str(kk),' | ',num2str(table(kk,:)))} ];
               #             end
               #         end

                       if (any(expval_table(table)<2)){
                fish = stats::fisher.test(table)
                pvalfish = fish$p.value
                zeros_rows = (apply(table,1,sum)==0)
                zeros_cols = (apply(table,2,sum)==0)
                if (sum(!zeros_rows)+sum(!zeros_cols)==0){
                  dffish = 0
                } else {
                  dffish = (sum(!zeros_rows)-1)*(sum(!zeros_cols)-1)
                }
    		        stafish = stats::qchisq(1-pvalfish, dffish)
                table_multi_3sm[where_in_table_3sm,1] = i
                table_multi_3sm[where_in_table_3sm,2] = l
                table_multi_3sm[where_in_table_3sm,3] = stafish + table_multi_3sm[where_in_table_3sm,3]
                table_multi_3sm[where_in_table_3sm,4] = dffish + table_multi_3sm[where_in_table_3sm,4]
                table_multi_3sm[where_in_table_3sm,5] = pvalfish
                table_multi_3sm[where_in_table_3sm,6] = 'Fisher'
                #            if verbosity>=3
                #                strtable=[ strtable {strcat('Associated test of the last table :',num2str([stafish pvalfish dffish]))}];
                #            end
                       } else {
                old.warn <- options()$warn # to suppress the warning messages
                options(warn = -1)
                chi2 = stats::chisq.test(table,correct=F)
                options(warn = old.warn)
                pvalchi2 = chi2$p.value
                dfchi2 = chi2$parameter
                stachi2 = chi2$statistic
                table_multi_3sm[where_in_table_3sm,1] = i
                table_multi_3sm[where_in_table_3sm,2] = l
                table_multi_3sm[where_in_table_3sm,3] = stachi2 + table_multi_3sm[where_in_table_3sm,3]
                table_multi_3sm[where_in_table_3sm,4] = dfchi2 + table_multi_3sm[where_in_table_3sm,4]
                table_multi_3sm[where_in_table_3sm,5] = pvalchi2 + table_multi_3sm[where_in_table_3sm,5]
                table_multi_3sm[where_in_table_3sm,6] = 'Chi-square'
                #            if verbosity>=3
                #                strtable=[ strtable {strcat('Associated test of the last table :',num2str(chi2(table)))}];
                #            end
                #        end
               }

                }
             }
                table_multi_3sm[where_in_table_3sm,5] = 1-stats::pchisq(table_multi_3sm[where_in_table_3sm,3],table_multi_3sm[where_in_table_3sm,4])
      } # if ns>1
            #stattotal = cbind(stattotal,table_multi_3sm[where_in_table_3sm,])
     }
     }
# compute overall test:
stat = sum(as.numeric(table_multi_3sm[,3]))
stat = round(stat,rounding)
dof = sum(as.numeric(table_multi_3sm[,4]))
pval = 1 - stats::pchisq(stat,dof)
pval = round(pval,rounding)
# if user specifies all outputs
if (verbose==TRUE) return(list(test3Gsm=c(stat=stat,df=dof,p_val=pval),details=table_multi_3sm))
# otherwise
if (verbose==FALSE) return(list(test3Gsm=c(stat=stat,df=dof,p_val=pval)))

}
