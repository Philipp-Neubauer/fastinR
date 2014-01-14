\name{multiplot}
\alias{multiplot}
\title{Plot density or denstrip plots from multiple MCMC runs}
\usage{
multiplot(MCMCouts, save = "fastinR_Multiplot", density = T)
}
\arguments{
  \item{MCMCouts}{A (named) list of objects produced by
  \code{\link{run_MCMC}}}

  \item{density}{If TRUE (default), density plots are
  drawn, if FALSE, denstrip plots drawn instead.}

  \item{save}{Either a string to be used as prefix for
  saved plots, or FALSE for disabling saving to file.}
}
\description{
A density/denstrip plot produced for multiple MCMC runs
from \code{\link{run_MCMC}}. Useful to compare across runs
with different data subsets or data types. The analyses
needs to be of the same type, for now only Analysis.Type =
Population.Proportions in \code{\link{run_MCMC}} works.
}
\details{
If plots are saved they are not drawn at the same time.
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

