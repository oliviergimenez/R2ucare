#' TestM.LTEC
#'
#' This function performs TestM.LTEC
#' @param X is a matrix of encounter histories with K occasions
#' @param freq is a vector of the number of individuals with the corresponding encounter history
#' @param verbose controls the level of the details in the outputs; default is TRUE for all details
#' @param rounding is the level of rounding for outputs; default is 3
#' @return This function returns a list with first component the overall test and second component a data.frame with occasion, the value of the test statistic, degree of freedom, p-value and test performed (chi-square, Fisher or none).
#' @author Olivier Gimenez <olivier.gimenez@cefe.cnrs.fr>, Roger Pradel, Rémi Choquet
#' @keywords package
#' @export
#' @examples
#' \dontrun{
#' # Read in Geese dataset:
#' geese = system.file("extdata", "geese.inp", package = "R2ucare")
#' geese = read_inp(geese)
#'
#' # Get encounter histories and number of individuals with corresponding histories
#' geese.hist = geese$encounter_histories
#' geese.freq = geese$sample_size
#'
#' # perform TestM.LTEC
#' testMltec(geese.hist, geese.freq)
#' }

testMltec <- function(X,freq,verbose=TRUE,rounding=3){


#derocc=1-min(filtre);
#filtre=find(filtre);

# various quantities to define
k = ncol(X)
his = X
a = max(his)

# initialization
table_multi_litec = data.frame(occasion = rep(NA,k-4),stat = rep(NA,k-4), df = rep(NA,k-4), p_val = rep(NA,k-4), test_perf = rep(FALSE,k-4))

marray = multimarray(X,freq)

debutligne = 1
finligne = seq(a,a*k,by=a)
debutcolonne = seq(1,a*(k-1),by=a)
fincolonne = a * (k-1)

datat = marray[,2:(ncol(marray)-1)] # extrait du m-array avec les revus, sans les relaches ni les jamais revus

#for i=2:(k-3)+derocc % boucle sur les occasions
for (i in 2:(k-3)){ # boucle sur les occasions

        mixandbases = datat[debutligne:finligne[i],debutcolonne[i]:fincolonne]

        for (j in 1:(i-2)){
            if ((i-2)<1) break
            mixandbases[1:a,] = mixandbases[1:a,] + mixandbases[(a+1):(2*a),]
            mixandbases = mixandbases[-((a+1):(2*a)),]
        }

    mixandbases = mixandbases[,-(1:a)]
        nk = nrow(mixandbases) # nb de melanges + nb de bases
        tabtemp = mixandbases[(nk-a+1):nk,] # les bases
        #tabtemp = tabtemp[filtre,] # les bases filtrees
        mixandbases = rbind(mixandbases[1:(nk-a),],tabtemp)
        nj = nrow(tabtemp)
        nk = nrow(mixandbases)

        if (any(apply(tabtemp,1,sum)==0) | (nj==0)){
            table_multi_litec[i-1,1] = i
            table_multi_litec[i-1,2] = 0
            table_multi_litec[i-1,3] = 0
            table_multi_litec[i-1,4] = 0
            table_multi_litec[i-1,5] = 'None'
            break
        }

        # pooling
        mixandbases = pooling_mixtures(nk,nj,a,mixandbases)

        # DEBUT DU CALCUL DU TEST DE MELANGE (traduit de Yantis et al)

        indic_mixbasis = c(matrix(0,nrow=nrow(mixandbases)-nj,ncol=1),rep(1,nj)) # indique 1 ou 0 suivant que base ou melange
        data = cbind(mixandbases,indic_mixbasis)
    #if verbosity>=3
    #    strtable=[ strtable {strcat('occasion :',num2str(i),', last column : 1 for basis and 0 for mixture ')}];
    #    taille=size(data,2);
    #    strtable=[ strtable {'--------------- Seen again at t+2 - Seen later '}];
    #    for kk=1:size(data,1)
    #        if data(kk,taille)==0
    #            strtable=[ strtable {strcat('When last released | ',num2str(data(kk,1:taille-1)))}];
    #        else
    #            strtable=[ strtable {strcat('Currently released | ',num2str(data(kk,1:taille-1)))}];
    #        end
    #    end
    #end

        nk = nrow(data)
        r = ncol(data)
        ni = r - 1
        data = t(data)
        nature = data[r,] # il s'agit de indic_mixanbasis'
        data = data[-r,] # on la supprime !!!
        tri = which(nature!=0) # coordonnees des bases
        nj = length(tri) # nombre de bases
        tri = c(tri,which(nature==0)) # ajout des coordonnees des melanges
        M = data[,tri] # on renumerote bases et melanges
        totk = apply(M,2,sum) # effectif des colonnes
        CoorMelVide = which(totk[(nj+1):nk]==0)
        if (!(length(CoorMelVide)==0)) M = M[,-CoorMelVide]
        nk = ncol(M)
        totk = apply(M,2,sum) # actualisation des effectifs des colonnes
        # Si aucune base n'est vide & si melanges
        if (nj!=nk){
        	# NEW definition des bases
        	Np = t(M[,1:nj])
        	# definition des melanges
        	Mp = t(M[,(nj+1):nk])
        	# calcul des coefficients du melanges
        	res = coef_mixtures(Mp,Np)
        	Q = res$P
        	P = res$PI
        	A = res$GAM
        	Q = rbind(P,Q)
        	# calcul des valeurs attendues
        	theoriques = matrix(rep(totk,ni),byrow=T,nrow=ni) * t(Q)
        	# calcul du nombre de degre de liberte
        	df = (nk-nj)*(ni-nj)
        	# test LR
        	tempchi2 = gof_test(1,c(M),c(theoriques))
        	table_multi_litec[i-1,1] = i
        	table_multi_litec[i-1,2] = tempchi2
        	table_multi_litec[i-1,3] = df
        	table_multi_litec[i-1,4] = 1 - stats::pchisq(tempchi2,df)
        	table_multi_litec[i-1,5] = 'Chi-square'
        	} else {
        	table_multi_litec[i-1,1] = i
        	table_multi_litec[i-1,2] = 0
        	table_multi_litec[i-1,3] = 0
        	table_multi_litec[i-1,4] = 0
        	table_multi_litec[i-1,5] = 'None'
        	}
}

# compute overall test:
stat = sum(as.numeric(table_multi_litec[,2]))
stat = round(stat,rounding)
dof = sum(as.numeric(table_multi_litec[,3]))
pval = 1 - stats::pchisq(stat,dof)
pval = round(pval,rounding)
# if user specifies all outputs
if (verbose==TRUE) return(list(testMltec=c(stat=stat,df=dof,p_val=pval),details=table_multi_litec))
# otherwise
if (verbose==FALSE) return(list(testMltec=c(stat=stat,df=dof,p_val=pval)))


} # function
