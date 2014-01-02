\name{multiplot}
\alias{multiplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot density or denstrip plots from multiple MCMC runs
%%  ~~function to do ... ~~
}
\description{A density/denstrip plot produced for multiple MCMC runs
  from \code{\link{run_MCMC}}. Useful to compare across runs with
  different data subsets or data types. The analyses needs to be of the
  same type, for now only Analysis.Type = Population.Proportions in
  \code{\link{run_MCMC}} works.
}
\usage{
multiplot(MCMCouts,density=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{MCMCouts}{A (named) list of objects produced by  \code{\link{run_MCMC}}}
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