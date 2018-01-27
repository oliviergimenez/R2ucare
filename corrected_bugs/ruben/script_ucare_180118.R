##### Load the database ####
mat <- c()
mat <- read.table("/Users/oliviergimenez/Dropbox/OG/GITHUB/R2ucare/corrected_bugs/ruben/example_matrix", sep="\t",
                  header = T,
                  row.names = 1)

# Put the 4 as 0 to assess that all time steps have individuals and that all individuals had a capture
mat[mat==4]<-0
mat <- as.matrix(mat)
dim(mat)
# (1) Assess that there is at least one individual in any time step
check <- as.vector(apply(mat, 2, sum))
sum(check==0)
cat("0 indicates that all sampling occasions had at least 1 individual")
cat("I put it to assess that this was not the error")

# (2) Assess that all individuals were captured at least once
check_indiv <- apply(mat, 1, sum); min(check_indiv)
cat("The minimum value should be 1")
if (min(check_indiv>0)){
  cat("Correct")
} else {
  cat("error")
}

# Create a vector indicating the number of individuals with this history
# As I create one row per individual I have to put one a 1
n_history <- rep(1, times=nrow(mat))
length(n_history)==nrow(mat)

library(R2ucare)
# Put again the zeros as fours for U_CARE
#mat[mat==0]<-4
test3Gsm(mat,n_history)
test3Gsr(mat,n_history)
test3Gwbwa(mat,n_history)
testMitec(mat,n_history)
testMltec(mat,n_history)

overall_JMV(mat,n_history)
