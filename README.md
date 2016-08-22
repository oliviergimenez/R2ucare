# R2ucare: An R package to test the goodness of fit of capture-recapture models

Author: [Olivier Gimenez](https://oliviergimenez.wordpress.com/)

Vignette: to come

Please email all comments/questions to olivier.gimenez [AT] cefe.cnrs.fr

Citation: to come

## What it does
Ths package contains R functions to perform goodness-of-fit tests for capture-recapture models. We will enrich it soon with various functions to manipulate capture-recapture data and visualize the tests. 

Despite what its name might suggest, **you do not need** to download and install U-CARE to run the R2ucare package. This package is basically a Matlab to R translation of U-CARE (Choquet et al. 2009). 

For Cormack-Jolly-Seber models (single-state), we refer to Lebreton et al. (1992) and Pradel et al. (2005) for the theory. For Arnason-Schwarz models (multistate), have a look to Pradel et al. (2003). [Chapter 5 of the Gentle Introduction to MARK](http://www.phidot.org/software/mark/docs/book/pdf/chap5.pdf) also provides a good start for understanding goodness-of-fit test. 

## To install the package

This repository hosts the development version of the package. It will also be available soon on CRAN (translation: I have to drastically reduce the to-do list below). For the time being, just use the working version:

```R
if(!require(devtools)) install.packages("devtools")
library("devtools")
install_github('oliviergimenez/R2ucare')
```

## Goodness-of-fit tests for the Cormack-Jolly-Seber model

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

## Goodness-of-fit tests for the Arnason-Schwarz model

```R
# read in the classical dipper dataset using package RMark
library(RMark)
geese = system.file("extdata", "geese.inp", package = "R2ucare")
geese = convert.inp(geese)

# add spaces between columns
geese.hist = matrix(as.numeric(unlist(strsplit(geese$ch, ''))),nrow=nrow(geese),byrow=T)
geese.freq = geese$freq

# encounter histories and number of individuals with corresponding histories
X = geese.hist
freq = geese.freq

# load R2ucare package
library(R2ucare)

# perform Test3.GSr, Test.3.GSm, Test3G.wbwa, TestM.ITEC and TestM.LTEC
test3Gsr(X,freq)
test3Gsm(X,freq)
test3Gwbwa(X,freq)
testMitec(X,freq)
testMltec(X,freq)
```

## References 

* Choquet, R., Lebreton, J.-D., Gimenez, O., Reboulet, A.-M., and R. Pradel. (2009). [U-CARE: Utilities for performing goodness of fit tests and manipulating CApture-REcapture data](https://dl.dropboxusercontent.com/u/23160641/my-pubs/Choquetetal2009UCARE.pdf). Ecography. 32: 1071-1074.
* Lebreton, J.-D. et al. (1992). Modeling survival and testing biological hypotheses using marked animals: a unified approach with case studies. Ecol. Monogr. 62: 67-118.
* Pradel, R., Gimenez O. and J.-D. Lebreton (2005). [Principles and interest of GOF tests for multistate capture-recapture models](https://dl.dropboxusercontent.com/u/23160641/my-pubs/Pradeletal2005ABC.pdf). Animal Biodiversity and Conservation 28: 189–204.
* Pradel R., Wintrebert C.M.A. and Gimenez O. (2003). [A proposal for a goodness-of-fit test to the Arnason-Schwarz multisite capture-recapture model](https://dl.dropboxusercontent.com/u/23160641/my-pubs/Pradeletal2003Biometrics.pdf). Biometrics 59: 43-53.

## Yet to be done

1. Stuff to write
    + add AS and JMV fitting to complete multisite gof test
    + add function to perform overall test (sum of all components)
    + add functions to manipulate data as in U-CARE
    + write vignette 
    + add signed tests
    + allow reading in files with Headed format
    + add gof for mixture of recapture and recoveries (Rachel)
    + add test for heterogeneity (Anita)
    + add verbose with the option to print tables
    + make readme similar to [stm](https://github.com/bstewart/stm)

2. Debugging
    + check df for test3Gsm
    + check all next/break conditions (unit tests)
    + check the behavior of the R chisq/Fisher.test functions when cells are deleted and the contingency table is no longer a square or a rectangle
    
3. Not so urgent
    + pass testMitec and testMltec in C++
    + plan a sequence of unit tests
    + add filtre like in MATLAB code (?)
    + add mosaic plot for wbwa, and other ways to visually represent contingency tables for the other tests
    + class ucare?
