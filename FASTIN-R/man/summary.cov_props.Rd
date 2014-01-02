\name{summary.cov_props}
\alias{summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Displays a summary (95\% credible interval and median) of posterior diet proportions
%%  ~~function to do ... ~~
}
\description{To be run with output from \code{\link{run_MCMC}}
}
\usage{
\method{summary}{cov_props}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{object}{MCMC output from \code{\link{run_MCMC}}, containing diet
  proportion MCMC chains}
   
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