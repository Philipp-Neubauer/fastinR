\name{add_SI}
\alias{add_SI}
\title{Add Stable Isotope data for predators and prey items}
\usage{
add_SI(SI.predators = NULL, SI.preys = NULL, Frac.Coeffs.mean = "",
  Frac.Coeffs.var = "", FC.mean = 1, FC.var = 1, datas = NULL)
}
\arguments{
  \item{SI.predators}{A data frame or a path to a csv file
  with the predator index/names (first column) and Stable
  Isotopes (1 row per predator), with Stable Isotope named
  across the first row.}

  \item{SI.preys}{A data frame or a path to a csv file with
  the prey names/sample id (first column) and SI
  measurements (1 row pey prey item), with Stable Isotope
  named across the first row. Samples from the same species
  should have the same prey name/sample id.}

  \item{Frac.Coeffs.mean}{A data frame or a path to a csv
  file with the prey specific additive fractionation
  coefficient means: Prey names (first column) and an n x P
  matrix for n preys and P Stable Isotopes}

  \item{Frac.Coeffs.var}{A data frame or a path to a csv
  file with the prey specific fractionation coefficient
  variances, dimensions as for the means}

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
  \code{\link{add_SI}}, needed if fatty acids and stable
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
SI.predators <- system.file("extdata", "Simdata_SI_preds.csv", package="fastinR")
SI.preys <- system.file("extdata", "Simdata_SI_preys.csv", package="fastinR")
Frac.Coeffs.mean <- system.file("extdata", "Simdata_SI_fc_means.csv", package="fastinR")
Frac.Coeffs.var <- system.file("extdata", "Simdata_SI_fc_var.csv", package="fastinR")
dats <- add_SI(SI.predators=SI.predators,SI.preys=SI.preys,Frac.Coeffs.mean=Frac.Coeffs.mean,Frac.Coeffs.var=Frac.Coeffs.var)
}
\author{
Philipp Neubauer
}
\references{
Neubauer,.P. and Jensen, O.P. (in prep)
}
\seealso{
\code{\link{add_FA}},\code{\link{add_Covs}},\code{\link{run_MCMC}},\code{\link{simulation}}
}

