% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MIEDA.R
\name{MIScatterPlot}
\alias{MIScatterPlot}
\title{Scatter Plot for Multiple Instance Regression}
\usage{
MIScatterPlot(tidydata, bag_size, true_primary = NULL, pred_primary = NULL)
}
\arguments{
\item{tidydata}{A list created from \code{\link{Tidy_dataset}}.}

\item{bag_size}{Maximum size of bags to be plotted. This parameter is for reducing the object size.}

\item{true_primary}{A list of logical vectors. Each logical vector denotes which instances are primary in a bag. In general, it is not available in real data.}

\item{pred_primary}{A list of logical vectors. Each logical vector is a fitted or predicted outcome that indicates which instances are primary in a bag.}
}
\value{
Return a ggplot object created from \code{ggplot}.
}
\description{
Visualize multiple instances and/or identification of primary instances.
}
\examples{
## Basic set-up
nsample <- 100
ninst <- 5
nprime <- 1 # do not change
nfeature <- 1

## Generate synthetic data
set.seed(6)
npoint <- (ninst - nprime) * nfeature
bag <- list()
for(i in 1:nsample){
  prime_inst <- matrix(runif(nfeature, -5, 5), ncol = nfeature)
  nonprime_inst <- matrix(c(rnorm(npoint, -5, 1),
                            rnorm(npoint, 5, 1)),
                          nrow = 2, byrow = TRUE)
  nonprime_inst <- nonprime_inst[cbind(sample(1:2, npoint, replace = TRUE), 1:npoint)]
  nonprime_inst <- matrix(nonprime_inst, ncol = nfeature, byrow = FALSE)
  bag[[i]] <- as.data.frame(rbind(prime_inst, nonprime_inst))
}
beta_true <- rep(2, nfeature + 1)
label <- unlist(lapply(bag, function(feature){
  beta_true[1] + as.matrix(feature[1,,drop = FALSE]) \%*\% beta_true[-1]
}))
label <- label + rnorm(nsample, mean = 0, sd = 1)

## Tidy data
tidydata <- Tidy_dataset(label = label,
                         feature_inst = bag)
                         
## BMIR model fitting
BMIR_fit <- BMIR_sampler(ntotal = 1000, tidydata = tidydata)
                         
## Visualize fitted results
tp <- lapply(1:tidydata$nsample, function(x) rep(c(TRUE, FALSE), c(1, ninst - 1)))
pp <- lapply(split(BMIR_fit$pip[,1], tidydata$membership), function(x) rank(-x, "min") <= 1)
MIScatterPlot(tidydata = tidydata, 
              bag_size = 5,
              true_primary = tp, 
              pred_primary = pp)
}
