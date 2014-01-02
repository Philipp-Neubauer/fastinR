\name{plot.cov_props}
\alias{plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot density or denstrip plots of posterior diet proportions
%%  ~~function to do ... ~~
}
\description{A pairwise plot to inspect correlations between prey items
  in the MCMC estiamtion, and a density/denstrip plot of diet proportion
  posterior distributions. To be run with output from \code{\link{run_MCMC}}
}
\usage{
\method{plot}{cov_props}(x,save="FASTIN_MCMC_",density=T,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{x}{MCMC output from \code{\link{run_MCMC}}, containing diet
  proportion MCMC chains}
  \item{save}{Either a string to be used as prefix for saved plots, or
    FALSE for disabling saving to file.}
  \item{density}{If TRUE (default), density plots are drawn, if FALSE,
  denstrip plots drawn instead.}
  
}
\details{

  If plots are saved they are not drawn at the same time.

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