## ---- message=FALSE, warning=FALSE---------------------------------------
if(!require(devtools)) install.packages("devtools")
library("devtools")
install_github('oliviergimenez/R2ucare')

## ---- message=FALSE, warning=FALSE---------------------------------------
if(!require(RMark)) install.packages("RMark")
library(RMark)
dipper = system.file("extdata", "ed.inp", package = "R2ucare")
dipper = convert.inp(dipper,group.df=data.frame(sex=c('Male','Female')))

## ---- message=FALSE, warning=FALSE---------------------------------------
dip.hist = matrix(as.numeric(unlist(strsplit(dipper$ch, ''))),nrow=nrow(dipper),byrow=T)

## ---- message=FALSE, warning=FALSE---------------------------------------
dip.freq = dipper$freq
dip.group = dipper$sex

## ---- message=FALSE, warning=FALSE---------------------------------------
mask = (dip.group == 'Female')
dip.fem.hist = dip.hist[mask,]
dip.fem.freq = dip.freq[mask]
mask = (dip.group == 'Male')
dip.mal.hist = dip.hist[mask,]
dip.mal.freq = dip.freq[mask]

## ---- message=FALSE, warning=FALSE---------------------------------------
library(R2ucare)

## ---- message=FALSE, warning=FALSE---------------------------------------
test3sr_females = test3sr(dip.fem.hist, dip.fem.freq)
test3sm_females = test3sm(dip.fem.hist, dip.fem.freq)
# we need the m-array to perform test2ct and test2cl
m = marray(dip.fem.hist,dip.fem.freq)$m[,,]
test2ct_females = test2ct(m)
test2cl_females = test2cl(m)
# display results:
test3sr_females
test3sm_females
test2ct_females
test2cl_females

## ---- message=FALSE, warning=FALSE---------------------------------------
overall_CJS(dip.fem.hist, dip.fem.freq)

## ---- message=FALSE, warning=FALSE---------------------------------------
library(RMark)
geese = system.file("extdata", "geese.inp", package = "R2ucare")
geese = convert.inp(geese)

## ---- message=FALSE, warning=FALSE---------------------------------------
geese.hist = matrix(as.numeric(unlist(strsplit(geese$ch, ''))),nrow=nrow(geese),byrow=T)
geese.freq = geese$freq

## ---- message=FALSE, warning=FALSE---------------------------------------
X = geese.hist
freq = geese.freq

## ---- message=FALSE, warning=FALSE---------------------------------------
library(R2ucare)

## ---- message=FALSE, warning=FALSE---------------------------------------
test3Gsr(X,freq)
test3Gsm(X,freq)
test3Gwbwa(X,freq)
testMitec(X,freq)
testMltec(X,freq)

## ---- message=FALSE, warning=FALSE---------------------------------------
overall_JMV(X,freq)

## ---- message=FALSE, warning=FALSE---------------------------------------
# Assuming the geese dataset has been read in R (see above):
X[X>1] = 1

## ---- message=FALSE, warning=FALSE---------------------------------------
# Assuming the geese dataset has been read in R (see above):
X[X==3]=2 # all 3s become 2s

## ---- message=FALSE, warning=FALSE,eval=FALSE----------------------------
#  # Assuming the male dipper dataset has been read in R (see above):
#  t(apply(X,1,rev))

## ---- message=FALSE, warning=FALSE,eval=FALSE----------------------------
#  # Assuming the male dipper dataset has been read in R (see above):
#  mask = (apply(X,1,sum)>0 & freq>0) # select non-empty histories, and histories with at least one individual
#  sum(!mask) # how many histories are to be dropped?
#  X[mask,] # drop these histories from dataset
#  freq[mask] # from counts

## ---- message=FALSE, warning=FALSE, eval=FALSE---------------------------
#  # Assuming the male dipper dataset has been read in R (see above):
#  X[,c(1,4,6)] # pick occasions 1, 4 and 6 (might be a good idea to clean the resulting dataset)
#  gather_146 = apply(X[,c(1,4,6)],1,max) # gather occasions 1, 4 and 6 by taking the max
#  X[,1] = gather_146 # replace occasion 1 by new occasion
#  X = X[,-c(4,6)] # drop occasions 4 and 6

## ---- message=FALSE, warning=FALSE, eval=FALSE---------------------------
#  # Assuming the geese dataset has been read in R (see above):
#  for (i in 1:nrow(X)){
#  occasion_first_encounter = min(which(X[i,]!=0)) # look for occasion of first encounter
#  X[i,occasion_first_encounter] = 0 # replace the first non zero by a zero
#  }
#  # delete empty histories from the new dataset
#  mask = (apply(X,1,sum)>0) # select non-empty histories
#  sum(!mask) # how many histories are to be dropped?
#  X[mask,] # drop these histories from dataset
#  freq[mask] # from counts

