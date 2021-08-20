# R2ucare <img src="man/figures/logo.png" align="right" width="190" height="200"/>

## What it does (and does not do)

This package contains `R` functions to perform goodness-of-fit tests for capture-recapture models. It also has various functions to manipulate capture-recapture data. Please email all suggestions for improvements, questions, comments and bugs to olivier.gimenez[at]cefe[dot]cnrs[dot]fr.

For Cormack-Jolly-Seber models (single-state), we refer to Lebreton et al. (1992) and Pradel et al. (2005) for the theory. 
For Arnason-Schwarz models (multistate), have a look to Pradel et al. (2003). [Chapter 5 of the Gentle Introduction to MARK](http://www.phidot.org/software/mark/docs/book/pdf/chap5.pdf) also provides a good start for understanding goodness-of-fit tests. 

**Warning**: to date, no goodness-of-fit test exists for models with individual covariates (unless you discretize them and use groups), individual time-varying covariates (unless you treat them as states) or temporal covariates; therefore, remove these covariates from your dataset before using it with `R2ucare`. For groups, just treat groups separately as in the Dipper example below. 

## To install the package

The latest stable version of the package can be downloaded from `CRAN` with the `R` command
```
install.packages("R2ucare")
```

The repository on `GitHub` https://github.com/oliviergimenez/R2ucare hosts the development version of the package, to install it:
```R
if(!require(devtools)) install.packages("devtools")
devtools::install_github("oliviergimenez/R2ucare")
```

Despite what its name might suggest, **you do not need** to download and install `U-CARE` to run the `R2ucare` package. 
This package is basically a `Matlab` to `R` translation of `U-CARE` (Choquet et al. 2009). 

## Getting started

The simplest way to get started is to have a look to the 
`R2ucare` [vignette](https://oliviergimenez.github.io/R2ucare/articles/vignette_R2ucare.html).

## How to cite R2ucare?

Please, cite our (free access) paper if you use `R2ucare`:

Gimenez, O, Lebreton, J-D, Choquet, R, Pradel, R. (2018) [R2ucare: An R package to perform goodness-of-fit tests for capture–recapture models](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13014). Methods in Ecology and Evolution 9: 1749-1754.

A bibtex entry is as follows:

```
@Article{,
  title = {R2ucare: An {R} package to perform goodness-of-fit tests for capture–recapture models},
  author = {{Gimenez} and {Olivier} and {Lebreton} and {Jean-Dominique} and {Choquet} and {Rémi} and {Pradel} and {Roger}},
  journal = {Methods in Ecology and Evolution},
  year = {2018},
  volume = {9},
  number = {7},
  pages = {1749-1754},
  url = {https://oliviergimenez.github.io/R2ucare/},
}
```

## To-do list

0. Fix bugs      
1. Add more tests (any help welcome)
    + gof test for recovery data
    + gof test for mixture of recapture and recoveries (based on Rachel McCrea's work)
    + gof tests for heterogeneity in detection and in transition propensity, test for underlying state-structure (based on Anita Jeyam's work)
    + Gof tests for closed pop models
    + Gof tests for occupancy models (relevant?)    
2. Make it more user-friendly
    + decision trees to suggest what to do when things go wrong 
    + pass testMitec, testMltec and AS/JMV model fitting in C++ using Rcpp
    + mosaic plot for wbwa, and other ways to visually represent contingency tables for the other tests
    + class ucare

## References 

* Choquet, R., Lebreton, J.-D., Gimenez, O., Reboulet, A.-M., and R. Pradel. (2009). [U-CARE: Utilities for performing goodness of fit tests and manipulating CApture-REcapture data](https://dl.dropboxusercontent.com/u/23160641/my-pubs/Choquetetal2009UCARE.pdf). Ecography. 32: 1071-1074.
* Lebreton, J.-D. et al. (1992). Modeling survival and testing biological hypotheses using marked animals: a unified approach with case studies. Ecol. Monogr. 62: 67-118.
* Pradel, R., Gimenez O. and J.-D. Lebreton (2005). [Principles and interest of GOF tests for multistate capture-recapture models](https://dl.dropboxusercontent.com/u/23160641/my-pubs/Pradeletal2005ABC.pdf). Animal Biodiversity and Conservation 28: 189–204.
* Pradel R., Wintrebert C.M.A. and Gimenez O. (2003). [A proposal for a goodness-of-fit test to the Arnason-Schwarz multisite capture-recapture model](https://dl.dropboxusercontent.com/u/23160641/my-pubs/Pradeletal2003Biometrics.pdf). Biometrics 59: 43-53.

