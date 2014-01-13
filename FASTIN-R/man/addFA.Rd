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

  \item{FC.mean}{optional - if no prey or sample specific
  fat content means are supplied in a fat.conts file, prey
  specific coefficients can be entered here using R's
  c(FC_1,FC_2,...) notation.}

  \item{FC.var}{optional - if no prey or sample specific
  fat content variances are supplied in a fat.conts file,
  prey specific coefficients can be entered here using R's
  c(FC_1,FC_2,...) notation.}

  \item{CC.mean}{optional - if no prey specific
  fractionation coefficiants are supplied via
  Conv.Coeffs.mean, CC.mean can provide either a global
  (single) mean coefficient or fatty acid specific mean
  coefficients using R's c(FA_1,FA_2,...) notation for ALL
  fatty acids.}

  \item{CC.var}{optional - if no prey specific
  fractionation coefficiants are supplied via
  Conv.Coeffs.mean, CC.var can provide either a global
  (single) coefficient variance or fatty acid specific
  coefficient variances using R's c(FA_1,FA_2,...) notation
  for ALL fatty acids.}

  \item{datas}{a data structure as produced by
  \code{\link{addSI}}, needed if fatty acids and stable
  isotopes are added sequentially.}

  \item{R.diag}{diagnoal of the prior for the predator
  covariance matrix. Note that this parameter can
  significantly influence convergence, especially if there
  are few predator signatures to estiamte the covariance
  from - handle with care! If too small, the sampler will
  get stuck in local modes and extreme values, if too high
  it can produce nonsenseical estiamtes where all
  proportions are equal (posterior mean at
  \code{1/n.preys})}
}
\description{
Files must be in .csv format.
}
\details{
Use \code{\link{simulation}} to simulate and write these
files to inspect the file structure.
}
\examples{
FA.predators <- system.file("extdata", "Simdata_FA_preds.csv", package="FASTIN")
FA.preys <- system.file("extdata", "Simdata_FA_preys.csv", package="FASTIN")
Conv.Coeffs.mean <- system.file("extdata", "Simdata_FA_cc_means.csv", package="FASTIN")
Conv.Coeffs.var <- system.file("extdata", "Simdata_FA_cc_var.csv", package="FASTIN")
fat.conts <- system.file("extdata", "Simdata_fat_cont.csv", package="FASTIN")
dats <- addFA(FA.predators=FA.predators,FA.preys=FA.preys,fat.conts=fat.conts,Conv.Coeffs.mean=Conv.Coeffs.mean,Conv.Coeffs.var=Conv.Coeffs.var)
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

