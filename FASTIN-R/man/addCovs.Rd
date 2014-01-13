\name{addCovs}
\alias{addCovs}
\title{Add Covariates and Group indices for predators}
\usage{
addCovs(Groups = "", Covariates = "")
}
\arguments{
  \item{Groups}{Index of group membership for each
  predator, one (named) column per grouping variable}

  \item{Covariates}{Covariate values for each predator, one
  (named) column per covariate}
}
\description{
Files must be in .csv format
}
\details{
Use \code{\link{simulation}} to simulate and write these
files to inpect the file structure.
}
\author{
Philipp Neubauer
}
\references{
Neubauer,.P. and Jensen, O.P. (in prep)
}
\seealso{
\code{\link{addFA}},\code{\link{addSI}},\code{\link{run_MCMC}},\code{\link{simulation}}
}

