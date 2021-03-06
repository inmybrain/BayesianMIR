% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MIEDA.R
\name{MISummarize}
\alias{MISummarize}
\title{Summarize Bags with Multiple Instances}
\usage{
MISummarize(tidydata, plot = FALSE, bag_size = 5)
}
\arguments{
\item{tidydata}{A list created from \code{\link{Tidy_dataset}}.}

\item{plot}{Logical. Draw a scatter plot using \code{\link{MIScatterPlot}}. Default is FALSE.}

\item{bag_size}{Maximum size of bags to be plotted. This parameter is for reducing the object size.}
}
\value{
No return value.
}
\description{
Summarize basic informations about instances in bags
}
\examples{
## Generate bags
tidydata <- Tidy_dataset(label = rnorm(3), feature_inst = replicate(3, matrix(rnorm(10), 2), simplify = FALSE))
MISummarize(tidydata)
}
\seealso{
MIScatterPlot
}
