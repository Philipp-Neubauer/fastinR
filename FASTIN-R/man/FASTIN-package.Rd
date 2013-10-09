\name{FASTIN-package}
\alias{FASTIN-package}
\alias{FASTIN}
\docType{package}
\title{
Fatty Acids and Stable Isotopes in R
}
\description{
This packages provides gui facilities to perfom Bayesian estiamtion of diet proportions based on fatty acid profiles, stable isotopes and their combination. Bayesian estiamtion is currently performed using the R package BRugs to call OpenBUGS from R. The package can handle covariates and grouped (e.g., regional) data as predictors of diet proportions, and will estiamte both individual and popualtion level diet protportions. Facilities for simulating, plotting and extracting/writing summary statistics to file are provided as well.
}
\details{
\tabular{ll}{
Package: \tab FASTIN\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-03-19\cr
License: \tab GPL (>= 2)\cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Philipp Neubauer

Maintainer: Philipp Neubauer <neubauer.phil@gmail.com>

}
\references{
Neubauer.P. and Jensen, O.P. (in prep)
}

\keyword{ package }
\seealso{
\code{\link{fastin}},\code{\link{simulation}},\code{\link{run_FASTIN}}
}
\examples{

simulation()
run_FASTIN()

}
