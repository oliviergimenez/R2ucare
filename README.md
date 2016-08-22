# R2ucare: An R package to test the goodness of fit of capture-recapture models

Authors: [Olivier Gimenez](https://oliviergimenez.wordpress.com/), Roger Pradel, Jean-Dominique Lebreton, Rémi Choquet, Anne-Marie Reboulet

Vignette: to come

Please email all comments/questions to olivier.gimenez [AT] cefe.cnrs.fr

Citation: to come

## What it does

Ths package contains R functions to perform goodness-of-fit tests for capture-recapture models. We will enrich it soon with various functions to manipulate capture-recapture data and visualize the tests. Any other function you would like to see available in R2ucare, just email us.

Despite what its name might suggest, **you do not need** to download and install U-CARE to run the R2ucare package. This package is basically a Matlab to R translation of U-CARE (Choquet et al. 2009). 
For Cormack-Jolly-Seber models (single-state), we refer to Lebreton et al. (1992) and Pradel et al. (2005) for the theory. For Arnason-Schwarz models (multistate), have a look to Pradel et al. (2003). [Chapter 5 of the Gentle Introduction to MARK](http://www.phidot.org/software/mark/docs/book/pdf/chap5.pdf) also provides a good start for understanding goodness-of-fit test. 

**Warning**: to date, no goodness-of-fit test exists for models with individual covariates (unless you discretize them and use groups), individual time-varying covariates (unless you treat them as states) or temporal covariates; therefore, remove these covariates from your dataset before using it with R2ucare. For groups, just treat the group separately as in the Dipper example below. 

## To install the package

This repository hosts the development version of the package. It will also be available soon on CRAN (I have to drastically reduce the to-do list below before submitting it; any help welcome!). For the time being, just use the working version:

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

# or all tests at once:
overall_CJS(X,freq)

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

# or all tests at once:
overall_CJS(X,freq)

```

## Goodness-of-fit tests for the Arnason-Schwarz model

```R
# read in the geese dataset using package RMark
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

# or all tests at once:
overall_JMV(X,freq)
```

## Various functions that might be useful

Several U-CARE functions to manipulate and process capture-recapture data can be mimicked with R built-in functions. For example, recoding multisite/state encounter histories in single-site/state ones is easy:
```R
# Assuming the geese dataset has been read in R (see above):
geese.hist[geese.hist>1] = 1
```
Also, reversing time is not that tough:
```R
# Assuming the male dipper dataset has been read in R (see above):
t(apply(dip.mal.hist,1,rev))
```

However, several useful functions need a proper R equivalent. I have coded a few of them, see the list below and ?name-of-the-function for more details (more to follow hopefully, suggestions and help welcome). 

* `marray` build the m-array for single-site/state capture-recapture data;
* `multimarray` build the m-array for multi-site/state capture-recapture data;
* `group_data` pool together individuals with the same encounter capture-recapture history.

## References 

* Choquet, R., Lebreton, J.-D., Gimenez, O., Reboulet, A.-M., and R. Pradel. (2009). [U-CARE: Utilities for performing goodness of fit tests and manipulating CApture-REcapture data](https://dl.dropboxusercontent.com/u/23160641/my-pubs/Choquetetal2009UCARE.pdf). Ecography. 32: 1071-1074.
* Lebreton, J.-D. et al. (1992). Modeling survival and testing biological hypotheses using marked animals: a unified approach with case studies. Ecol. Monogr. 62: 67-118.
* Pradel, R., Gimenez O. and J.-D. Lebreton (2005). [Principles and interest of GOF tests for multistate capture-recapture models](https://dl.dropboxusercontent.com/u/23160641/my-pubs/Pradeletal2005ABC.pdf). Animal Biodiversity and Conservation 28: 189–204.
* Pradel R., Wintrebert C.M.A. and Gimenez O. (2003). [A proposal for a goodness-of-fit test to the Arnason-Schwarz multisite capture-recapture model](https://dl.dropboxusercontent.com/u/23160641/my-pubs/Pradeletal2003Biometrics.pdf). Biometrics 59: 43-53.

## Yet to be done

1. Stuff to write
    + add more functions to manipulate data as in U-CARE (clean data, suppress first encounter, ungroup, ...)
    + write vignette 
    + add signed tests
    + allow reading in files with Headed format (see [here](https://github.com/hadley/haven), [here](http://blog.analytixware.com/2015/02/reading-sas-into-r.html) and [here](http://rconvert.com/sas-vs-r-code-compare/5-ways-to-convert-sas-data-to-r/))
    + add gof for mixture of recapture and recoveries (Rachel)
    + add test for heterogeneity (Anita)
    + add verbose with the option to print tables
    + add AS and JMV fitting to complete multisite gof test (https://github.com/oliviergimenez/multievent_jags_R; translate in C++ using Rcpp to speed up fitting process)

2. Debugging
    + check df for test3Gsm
    + check all next/break conditions (unit tests)
    + check the behavior of the R chisq/Fisher.test functions when cells are deleted and the contingency table is no longer a square or a rectangle
    
3. Not so urgent
    + pass testMitec and testMltec in C++ using Rcpp
    + plan a sequence of unit tests
    + add filtre like in MATLAB code for ghost states
    + add mosaic plot for wbwa, and other ways to visually represent contingency tables for the other tests
    + class ucare?
