\name{plotvarselect}
\alias{plotvarselect}
\title{Plot cosine distance and matrix condition number for variable selection}
\usage{
plotvarselect(prey_mat, prey.ix)
}
\arguments{
  \item{prey_mat}{(row-wise) dataframe of prey
  compositions}

  \item{prey.ix}{Indexes rows of \code{prey_mat}}
}
\value{
Two plots with 1) relative contributions of Fatty Acids to
source separation and 2) matrix condition number.
}
\description{
Cosine distance and matrix condition number are graphed as
a function of Fatty Acids, where fatty acids are ordered by
their relative contribution to Canonical axes in a CAP
(Canonical Analysis of Principal Coordinates).
}
\examples{
\dontrun{
# load simulated example
data('Sim')
preys <- datas$datas.FA$preys
prey.ix <- datas$prey.ix
plotvarselect(preys,prey.ix)
}
}
\author{
Philipp Neubauer
}
\references{
Neubauer.P. and Jensen, O.P. (in prep) #'
}
\seealso{
\code{\link{addFA}},\code{\link{selectvars}},\code{\link{run_MCMC}},\code{\link{simulation}}
}

