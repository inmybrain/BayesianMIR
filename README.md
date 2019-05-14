
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- To add a badge  -->
<!-- [![Travis-CI Build Status](https://travis-ci.org/geanders/countyweather.svg?branch=master)](https://travis-ci.org/geanders/countyweather) -->
<!-- [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/countyweather)](https://cran.r-project.org/package=countyweather) -->
BayesianMIR
===========

This is a README file of the R package *BayesianMIR*. In our paper, not yet available, we develop the Bayesian multi-instance regression (BMIR) model, applied to the multiple instance regression problem. More details about the structure of data and the Bayesian modeling, see our upcoming paper.

<!-- We assume the primary instance assumption used in @RayPage2005, that is, there is one primary instance explaining the bag-level response variable. -->
Installation of the package
---------------------------

To install our package, you may simply execute the following codes:

``` r
# install.packages("devtools")
devtools::install_github("inmybrain/BayesianMIR", subdir = "BayesianMIR") # don't forget to specify subdir!
```

If you come across a problem like [this](https://github.com/r-lib/remotes/issues/130), please refer to [this answer](https://github.com/r-lib/remotes/issues/130#issuecomment-423830669) in that issue.

Or you can install the source file using the command line after downloading it from [here](XXX) (under construction);

``` bash
R CMD INSTALL BayesianMIR_1.0.tar.gz
```

A basic example of using the package
------------------------------------

We give a toy example to apply the main function `BMIR_sampling` to do sampling from the posterior distribution.

<!-- First, generate $40$ samples from the Gaussian mixture model with two components. One of the component has zero mean vector, while the other a vector parallel to $(1, \cdots, 1)^{\rm T}$ with modulus $5$. We assume the heterogeneous covariance structure between them. -->
``` r
library("BayesianMIR")
# K <- 2
# p <- 2
# q <- 3
# U <- lapply(1:K, function(noarg) getCovariance(p, 0.3, "AR"))
# V <- lapply(1:K, function(noarg) getCovariance(q, 0.2, "CS"))
# Sigma <- Map(kronecker, U, V) # separable covariance matrix
# mu <- list(rep(0, p * q), 5 / sqrt(p*q) * rep(1, p * q)) # distinct mean vectors
# Y <- vector(mode = "list", length = K)
# set.seed(6) # to make results reproducible
# for(i in 1:K){
#  Y[[i]] <- mvtnorm::rmvnorm(n = 20, mean = mu[[i]], sigma = Sigma[[i]])
# }
```

Notes
-----

<!-- - For available covariance structures, see the help page; -->
<!-- ```{r, eval = FALSE} -->
<!-- ?Mclust_SEP_cpp -->
<!-- ``` -->
<!-- - As for initial assignment of cluster membership, each sample is assigned randomly to clusters. -->
Issues
------

We are happy to troubleshoot any issue with the package;

-   please contact to the maintainer by <seongohpark6@gmail.com>, or

-   please open an issue in the github repository.

<!-- ## Error and warning messages you may get -->
<!-- ## References  -->
