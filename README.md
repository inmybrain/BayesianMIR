
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

<!-- Or you can install the source file using the command line after downloading it from [here](XXX) (NOT AVAILABLE NOW); -->
<!-- ```{bash, eval = FALSE} -->
<!-- R CMD INSTALL BayesianMIR_1.0.tar.gz -->
<!-- ``` -->
A basic example of using the package
------------------------------------

### Generate the Monte Carlo Markov Chains (model estimation)

We give a toy example to apply the main function `BMIR_sampler` to do sampling from the posterior distribution.

To do it, we first generate multiple instance regression data.

``` r
rm(list = ls())

nsample <- 200
ninst <- 5
nprime <- 1 # do not change
nfeature <- 5

set.seed(6)
npoint <- (ninst - nprime) * nfeature
bag <- list()
for(i in 1:nsample){
  prime_inst <- matrix(runif(nfeature, -5, 5), ncol = nfeature)
  nonprime_inst <- (-1)^rbinom(npoint, 1, 1/2) * runif(npoint, min = -10, max = 10)
  nonprime_inst <- matrix(nonprime_inst, ncol = nfeature, byrow = F)
  
  # nonprime_inst <- matrix(c(rnorm(npoint, -5, 1),
  #                           rnorm(npoint, 5, 1)),
  #                         nrow = 2, byrow = T)
  # nonprime_inst <- nonprime_inst[cbind(sample(1:2, npoint, replace = T), 1:npoint)]
  # nonprime_inst <- matrix(nonprime_inst, ncol = nfeature, byrow = F)
  
  bag[[i]] <- as.data.frame(rbind(prime_inst, nonprime_inst))
}
beta_true <- rep(2, nfeature + 1)
label <- unlist(lapply(bag, function(feature){
  beta_true[1] + as.matrix(feature[1,,drop = FALSE]) %*% beta_true[-1]
}))
label <- label + rnorm(nsample, mean = 0, sd = 1)
```

For convenience, we use `Tidy_dataset` function to tidy data.

``` r
library(BayesianMIR)
#> Warning: replacing previous import 'ggplot2::margin' by
#> 'randomForest::margin' when loading 'BayesianMIR'
tidydata <- Tidy_dataset(label = label[1:100],
                         feature_inst = bag[1:100])
newtidydata <- Tidy_dataset(feature_inst = bag[-(1:100)])
```

We can obtain the MCMC samples using `BMIR_sampler`.

``` r
## BMIR model fitting
ntotal <- 10000
BMIR_fit <- BMIR_sampler(ntotal = ntotal, tidydata = tidydata)
#> =============================================================
#> Multiple Instance Bayesian Regression
#> Elapsed time for chain1=0.017 mins: MCMC sampling is done!
```

### Visualization

Using the fitted model, a scatter plot for multiple instance regression can be provided as follows.

``` r
MIScatterPlot(tidydata = tidydata, 
              bag_size = 5,
              true_primary = lapply(1:tidydata$nsample, function(x) rep(c(T,F), c(1, ninst - 1))), 
              pred_primary = lapply(split(BMIR_fit$pip[,1], tidydata$membership), function(x) rank(-x, ties.method = "min") <= 1)
)
```

![](README-unnamed-chunk-6-1.png)

Using slightl modified `ggmcmc::ggs_density` function, we can have the Bayesian inference.

``` r
# install.packages("ggmcmc")
library("ggmcmc")
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
#> Loading required package: tidyr
#> Loading required package: ggplot2
ggs_density <- function (D, ncol, family = NA, rug = FALSE, hpd = FALSE, greek = FALSE) 
  ## - ncol is added!
  ## - ci -> ggmcmc::ci
  ## - [Low, High] interval is commented
{
  if (!is.na(family)) {
    D <- get_family(D, family = family)
  }
  if (attributes(D)$nChains <= 1) {
    f <- ggplot(D, aes(x = value))
  }
  else {
    f <- ggplot(D, aes(x = value, colour = as.factor(Chain), 
                       fill = as.factor(Chain)))
  }
  f <- f + geom_density(alpha = 0.3) + scale_fill_discrete(name = "Chain") + 
    scale_colour_discrete(name = "Chain")
  if (!greek) {
    f <- f + facet_wrap(~Parameter, ncol = ncol, scales = "free")
  }
  else {
    f <- f + facet_wrap(~Parameter, ncol = ncol, scales = "free", 
                        labeller = label_parsed)
  }
  if (rug) 
    f <- f + geom_rug(alpha = 0.1)
  if (hpd) {
    ciD <- ggmcmc::ci(D)
    f <- f + geom_segment(data = ciD, size = 2, color = "blue", inherit.aes = FALSE, 
                          aes(x = low, xend = high, y = 0, yend = 0)) 
    # +geom_segment(
    #   data = ciD,
    #   size = 1,
    #   inherit.aes = FALSE,
    #   aes(
    #     x = Low,
    #     xend = High,
    #     y = 0,
    #     yend = 0
    #   )
    # )
  }
  return(f)
}
ggs_mcmc <- ggmcmc::ggs(BMIR_fit$mcmclist)
ggs_mcmc$Parameter <- factor(ggs_mcmc$Parameter, labels = c(paste0("coef", 1:(nfeature + 1)), "sig2_error"))
ggs_density(ggs_mcmc %>% 
              filter(Iteration > ntotal * 1 / 4), 
            ncol = 2,
            hpd = TRUE) + 
  geom_vline(data = data.frame(Parameter = c(paste0("coef", 1:(nfeature + 1)), "sig2_error"),
                               true_val = c(rep(2, 1 + nfeature), 1)),
             aes(xintercept = true_val), color = "red") +
  labs(x = "Value", y = "Density") + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
```

![](README-unnamed-chunk-7-1.png)

### Prediction in new bags

When new bags (i.e. without labels) are given, we can predict both labels and primary instances using `predict.BMIR`.

``` r
pred_fit <- predict.BMIR(BMIRchain = BMIR_fit$mcmclist$Chain1, 
                         pip = BMIR_fit$pip[,1], 
                         tidydata = tidydata, 
                         newtidydata = newtidydata, 
                         k = 1)
```

Let us see how prediction works.

``` r
ggplot(data = data.frame(pred = pred_fit$newtidydata$label, 
                         true = label[-(1:100)]), 
       mapping = aes(x = pred, y = true)) + 
  geom_point() + geom_abline(intercept = 0, slope = 1, color = "red")
```

![](README-unnamed-chunk-9-1.png)

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