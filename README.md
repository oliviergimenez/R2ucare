# R2ucare: An R package to test the goodness of fit of capture-recapture models

Authors: [Olivier Gimenez](https://oliviergimenez.wordpress.com/), Jean-Dominique Lebreton, Rémi Choquet, Roger Pradel

Please email all comments/questions to olivier.gimenez [AT] cefe.cnrs.fr

Citation: to come

## What it does (and does not do)

Ths package contains R functions to perform goodness-of-fit tests for capture-recapture models. It also has various functions to manipulate capture-recapture data. Please email all suggestions for improvements, questions, comments and bugs to olivier.gimenez [AT] cefe.cnrs.fr.

For Cormack-Jolly-Seber models (single-state), we refer to Lebreton et al. (1992) and Pradel et al. (2005) for the theory. For Arnason-Schwarz models (multistate), have a look to Pradel et al. (2003). [Chapter 5 of the Gentle Introduction to MARK](http://www.phidot.org/software/mark/docs/book/pdf/chap5.pdf) also provides a good start for understanding goodness-of-fit test. 

**Warning**: to date, no goodness-of-fit test exists for models with individual covariates (unless you discretize them and use groups), individual time-varying covariates (unless you treat them as states) or temporal covariates; therefore, remove these covariates from your dataset before using it with R2ucare. For groups, just treat the group separately as in the Dipper example below. 

## To install the package

This repository hosts the development version of the package. It will also be available soon on CRAN (I have to drastically reduce the to-do list below before submitting it; any help welcome!). For the time being, just use the working version:

```R
if(!require(devtools)) install.packages("devtools")
library("devtools")
install_github('oliviergimenez/R2ucare')
```

Despite what its name might suggest, **you do not need** to download and install U-CARE to run the R2ucare package. This package is basically a Matlab to R translation of U-CARE (Choquet et al. 2009). 

## Getting started

The simplest way to get started is to have a look to the [R2ucare vignette](https://github.com/oliviergimenez/R2ucare/blob/master/inst/doc/vignette_R2ucare.Rmd).

## References 

* Choquet, R., Lebreton, J.-D., Gimenez, O., Reboulet, A.-M., and R. Pradel. (2009). [U-CARE: Utilities for performing goodness of fit tests and manipulating CApture-REcapture data](https://dl.dropboxusercontent.com/u/23160641/my-pubs/Choquetetal2009UCARE.pdf). Ecography. 32: 1071-1074.
* Lebreton, J.-D. et al. (1992). Modeling survival and testing biological hypotheses using marked animals: a unified approach with case studies. Ecol. Monogr. 62: 67-118.
* Pradel, R., Gimenez O. and J.-D. Lebreton (2005). [Principles and interest of GOF tests for multistate capture-recapture models](https://dl.dropboxusercontent.com/u/23160641/my-pubs/Pradeletal2005ABC.pdf). Animal Biodiversity and Conservation 28: 189–204.
* Pradel R., Wintrebert C.M.A. and Gimenez O. (2003). [A proposal for a goodness-of-fit test to the Arnason-Schwarz multisite capture-recapture model](https://dl.dropboxusercontent.com/u/23160641/my-pubs/Pradeletal2003Biometrics.pdf). Biometrics 59: 43-53.

## Yet to be done

1. Stuff to write
    + add gof for recoveries
    + add gof for mixture of recapture and recoveries (Rachel)
    + add test for heterogeneity (Anita)
    + add verbose with the option to print tables

2. Not so urgent
    + pass testMitec, testMltec and AS/JMV model fitting in C++ using Rcpp
    + plan a sequence of unit tests
    + add mosaic plot for wbwa, and other ways to visually represent contingency tables for the other tests
    + gof tests for closed pop models, and occupancy models
    + class ucare?
