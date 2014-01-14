\name{selectvars}
\alias{selectvars}
\title{Select Fatty Acids from imported data (should really be a method for subset...)}
\usage{
selectvars(datas, ix = NULL)
}
\arguments{
  \item{datas}{A data structure produced by
  \code{\link{addSI}} and \code{\link{addFA}}}

  \item{ix}{Index of variables to retain, if not supplied,
  \code{\link{plotvarselect}} is called to itneractively
  select variables based on cosine distance and prey matrix
  condition.}
}
\value{
a data structure of the same form as datas, with Fatty
Acids selected by ix.
}
\description{
Select Fatty Acids from imported data (should really be a
method for subset...)
}
\examples{
\dontrun{
# load simulated example
data('Sim')
selectvars(datas,4:1)
# or select visualy
selectvars(datas)
}
}
\author{
Philipp Neubauer
}
\references{
Neubauer.P. and Jensen, O.P. (in prep)
}
\seealso{
\code{\link{addFA}},\code{\link{plotvarselect}},\code{\link{run_MCMC}},\code{\link{simulation}}
}

