\name{DietProportionPlot}
\alias{DietProportionPlot}
\title{Plot density or denstrip plots of posterior diet proportions}
\arguments{
  \item{x}{MCMC output from \code{\link{run_MCMC}},
  containing diet proportion MCMC chains}

  \item{save}{Either a string to be used as prefix for
  saved plots, or FALSE for disabling saving to file.}

  \item{density}{If TRUE (default), density plots are
  drawn, if FALSE, denstrip plots drawn instead.}
}
\description{
A pairwise plot to inspect correlations between prey items
in the MCMC estiamtion, and a density/denstrip plot of diet
proportion posterior distributions. To be run with output
from \code{\link{run_MCMC}}
}
\details{
If plots are saved they are not drawn at the same time.
That may change in the future...
}
\author{
Philipp Neubauer
}
\references{
Neubauer.P. and Jensen, O.P. (in prep)
}
\seealso{
\code{\link{run_MCMC}},\code{\link{diags}}
}

