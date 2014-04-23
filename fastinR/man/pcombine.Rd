\name{pcombine}
\alias{pcombine}
\title{Combine MCMC chains for prey items}
\usage{
pcombine(mcmc, index)
}
\arguments{
  \item{mcmc}{output from \code{\link{run_MCMC}}}

  \item{index}{A vector of integer indices or prey names to
  combine.}
}
\value{
An object of the same class as the input with combined MCMC
chains.
}
\description{
Prey items can be combined into groups post hoc by summing
their individual MCMC draws
}
\author{
Philipp Neubauer
}
\references{
Neubauer.P. and Jensen, O.P. (in prep)
}
\seealso{
\code{\link{run_MCMC}}
}

