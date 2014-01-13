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
files to inspect the file structure.
}
\examples{
Cov_path <- system.file("extdata", "Simdata_Covariates.csv", package="FASTIN")
Group_path<- system.file("extdata", "Simdata_Groups.csv", package="FASTIN")
addCovs(Groups=Group_path,Covariates=Cov_path)
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

