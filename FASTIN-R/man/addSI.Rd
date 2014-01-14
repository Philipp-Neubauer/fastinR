\name{addSI}
\alias{addSI}
\title{Add Stable Isotope data for predators and prey items}
\usage{
addSI(SI.predators = NULL, SI.preys = NULL, Frac.Coeffs.mean = "",
  Frac.Coeffs.var = "", FC.mean = 1, FC.var = 1, datas = NULL)
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

  \item{Frac.Coeffs.var}{Prey specific Fractionation
  coefficient variances, dimensions as for the means}

  \item{FC.mean}{optional - if no prey specific
  fractionation coefficiants are supplied via
  Frac.Coeffs.mean, FC mean can provide either a global
  (single) mean coefficient or fatty acid specific mean
  coefficients using R's c(FA_1,FA_2,...) notation for ALL
  fatty acids.}

  \item{FC.var}{optional - if no prey specific
  fractionation coefficiants are supplied via
  Frac.Coeffs.mean, FC var can provide either a global
  (single) coefficient variance or fatty acid specific
  coefficient variances using R's c(FA_1,FA_2,...) notation
  for ALL fatty acids.}

  \item{datas}{a data structure as produced by
  \code{\link{addSI}}, needed if fatty acids and stable
  isotopes are added sequentially.}
}
\description{
Files must be in .csv format.
}
\details{
Use \code{\link{simulation}} to simulate and write these
files to inspect the file structure.
}
\examples{
SI.predators <- system.file("extdata", "Simdata_SI_preds.csv", package="FASTIN")
SI.preys <- system.file("extdata", "Simdata_SI_preys.csv", package="FASTIN")
Frac.Coeffs.mean <- system.file("extdata", "Simdata_SI_fc_means.csv", package="FASTIN")
Frac.Coeffs.var <- system.file("extdata", "Simdata_SI_fc_var.csv", package="FASTIN")
dats <- addSI(SI.predators=SI.predators,SI.preys=SI.preys,Frac.Coeffs.mean=Frac.Coeffs.mean,Frac.Coeffs.var=Frac.Coeffs.var)
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

