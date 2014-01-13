\name{addSI}
\alias{addSI}
\title{Add Stable Isotope data for predators and prey items}
\usage{
addSI(SI.predators = NULL, SI.preys = NULL, Frac.Coeffs.mean = "",
  Frac.Coeffs.var = "", FC.mean = 1, FC.var = 1, R.diag.SI = 0.2,
  datas = NULL)
}
\arguments{
  \item{SI.predators}{Predator index/names (first column)
  and Stable Isotopes (1 row pey predator), with Stable
  Isotope named across the first row}

  \item{SI.preys}{Prey names/sample id (first column) and
  fatty acid profiles (1 row pey prey item), ith Stable
  Isotope named across the first row}

  \item{Frac.Coeffs.mean}{Prey specific additive
  fractionation coefficient means: Prey names (first
  column) and an n x P matrix for n preys and P Stable
  Isotopes}

  \item{FracCoeffs.var}{Prey specific Fractionation
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
\code{\link{addFA}},\code{\link{addCovs}},\code{\link{run_MCMC}},\code{\link{simulation}}
}

