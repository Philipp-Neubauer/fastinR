\name{diags}
\alias{diags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Diagnostics for MCMC runs
%%  ~~function to do ... ~~
}
\description{Raftery-Lewis and Gelman-Rubin diagnostics for MCMC chains
  from \code{\link{run_MCMC}}. Multiple chains are needed for the latter
  type of diagnostic.
}
\usage{
diags(MCMCout=NULL,accuracy=0.01,proba=0.95,quant=0.025)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{MCMCouts}{An object produced by \code{\link{run_MCMC}}}
  \item{accuracy}{Accuracy with which parameters are to be estiamted}
  \item{proba}{probability with which estiamtes are within the itnerval [quant,1-quant]}
  \item{quant}{interval within which to estiamte the parameter}
  
}
\details{

  The function allows to diagnose if a given MCMC run was liekly long
  enough to have produced reliable estiamtes. The number of iterations
  and a thinning interval are suggested - if these are significantly
  larger than the parameters used to produce MCMCout, then the chains
  should be re-run with the suggested parameters.

}
\references{

Neubauer.P. and Jensen, O.P. (in prep)

}
\author{

Philipp Neubauer

}
\seealso{
\code{\link{run_MCMC}},\code{\link{diags}}
}