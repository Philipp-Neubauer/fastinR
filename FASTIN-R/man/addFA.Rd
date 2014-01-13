\name{addFA}
\alias{addFA}
\title{Add Fatty Acid profile data for predators and prey items}
\usage{
addFA(FA.predators = NULL, FA.preys = NULL, fat.conts = "",
  Conv.Coeffs.mean = "", Conv.Coeffs.var = "", FC.mean = 1, FC.var = 1,
  CC.mean = 1, CC.var = 1, R.diag = 0.2, datas = NULL)
}
\arguments{
  \item{FA.predators}{Predator index/names (first column)
  and fatty acid profiles (1 row pey predator), with fatty
  acids named across the first row}

  \item{FA.preys}{Prey names/sample id (first column) and
  fatty acid profiles (1 row pey prey item), with fatty
  acids names across the first row}

  \item{fat.conts}{Prey fat contents, as (columnwise) mean
  and variance per prey species or specified for each prey
  sample for the main analysis, in that case the first
  column is the prey sample id id and the second column is
  the individual sample's fat content}

  \item{Conv.Coeffs.mean}{Prey specific conversion
  coefficient means: Prey names (first column) and an n x P
  matrix for n preys and P fatty acids}

  \item{Conv.Coeffs.var}{Prey specific conversion
  coefficient variances, dimensions as for the means}
}
\description{
Files must be in .csv format.
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
\code{\link{addSI}},\code{\link{addCovs}},\code{\link{selectvars}},\code{\link{run_MCMC}},\code{\link{simulation}}
}

