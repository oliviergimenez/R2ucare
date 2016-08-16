# R package R2ucare

## What it does
Ths package contains R functions to perform goodness-of-fit tests for capture-recapture models (and various manipulations on these data). 
This is basically a Matlab to R translation of U-CARE (Choquet et al. 2009). 
For Cormack-Jolly-Seber models (single-state), we refer to Lebreton et al. (1992) and Pradel et al. (2005) for the theory. 
For Arnason-Schwarz models (multistate), have a look to Pradel et al. (2003). 
[Chapter 5 of the Gentle Introduction to MARK](http://www.phidot.org/software/mark/docs/book/pdf/chap5.pdf) also provides a good start for understanding goodness-of-fit test.

## To install the package
```R
install.packages("devtools")
library("devtools")
install_github('oliviergimenez/R2ucare')
```

## A session example

```R
# read in the classical dipper dataset using package RMark
library(RMark)
dipper = system.file("extdata", "ed.inp", package = "R2ucare")
dipper = convert.inp(dipper,group.df=data.frame(sex=c('Male','Female')))

# add spaces between columns
dip.hist = matrix(as.numeric(unlist(strsplit(dipper$ch, ''))),nrow=nrow(dipper),byrow=T)

# get counts and encounter histories
dip.freq = dipper$freq
dip.group = dipper$sex

# split tha dataset in females/males
mask = (dip.group == 'Female')
dip.fem.hist = dip.hist[mask,]
dip.fem.freq = dip.freq[mask]
mask = (dip.group == 'Male')
dip.mal.hist = dip.hist[mask,]
dip.mal.freq = dip.freq[mask]

# load R2ucare package
library(R2ucare)

# perform Test.3Sr, Test3.Sm, Test2.Ct and Test.Cl for females
test3sr_females = test3sr(dip.fem.hist, dip.fem.freq)
test3sm_females = test3sm(dip.fem.hist, dip.fem.freq)
X = dip.fem.hist
freq = dip.fem.freq
m = marray(X,freq)$m[,,]
test2ct_females = test2ct(m)
test2cl_females = test2cl(m)
test3sr_females
test3sm_females
test2ct_females
test2cl_females

# perform Test.3Sr, Test3.Sm, Test2.Ct and Test.Cl for males
test3sr_males = test3sr(dip.mal.hist, dip.mal.freq)
test3sm_males = test3sm(dip.mal.hist, dip.mal.freq)
X = dip.mal.hist
freq = dip.mal.freq
m = marray(X,freq)$m[,,]
test2ct_males = test2ct(m)
test2cl_males = test2cl(m)
test3sr_males
test3sm_males
test2ct_males
test2cl_males
```

## References 

* Choquet, R., Lebreton, J.-D., Gimenez, O., Reboulet, A.-M., and R. Pradel. (2009). [U-CARE: Utilities for performing goodness of fit tests and manipulating CApture-REcapture data](https://dl.dropboxusercontent.com/u/23160641/my-pubs/Choquetetal2009UCARE.pdf). Ecography. 32: 1071-1074.
* Lebreton, J.-D. et al. (1992). Modeling survival and testing biological hypotheses using marked animals: a unified approach with case studies. Ecol. Monogr. 62: 67-118.
* Pradel, R., Gimenez O. and J.-D. Lebreton (2005). [Principles and interest of GOF tests for multistate capture-recapture models](https://dl.dropboxusercontent.com/u/23160641/my-pubs/Pradeletal2005ABC.pdf). Animal Biodiversity and Conservation 28: 189–204.
* Pradel R., Wintrebert C.M.A. and Gimenez O. (2003). [A proposal for a goodness-of-fit test to the Arnason-Schwarz multisite capture-recapture model](https://dl.dropboxusercontent.com/u/23160641/my-pubs/Pradeletal2003Biometrics.pdf). Biometrics 59: 43-53.

## Yet to be done
* check the behavior of the R chisq/Fisher.test functions when cells are deleted and the contingency table is no longer a square or a rectangle
* add signed test
* add functions to manipulate data
* write vignette 
* add function to perform overall test (sum of all components)
* allow reading in files with Headed format
* add multistate gof
* add gof for mixture of recapture and recoveries (Rachel)
* add test for heterogeneity (Anita)
* plan a sequence of root tests
* add a series of datasets for end-users to manipulate package
