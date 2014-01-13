\name{dataplot}
\alias{dataplot}
\title{Plot Non-metric multidimensional scaling (NMDS) or (clr transformed) diet data}
\usage{
dataplot(datas = NULL)
}
\arguments{
  \item{datas}{Diet data produced by \code{\link{addFA}}
  and/or \code{\link{addSI}}}
}
\description{
Diet data produced by \code{\link{addFA}} and/or
\code{\link{addFA}} is projected onto 2 dimesions using
NMDS, where Fatty Acid data is first clr transformed and
then concatenated with SI data to produce the plot.
}
\author{
Philipp Neubauer
}
\references{
Neubauer,.P. and Jensen, O.P. (in prep)
}
\seealso{
\code{\link{addFA}},\code{\link{addSI}},\code{\link{selectvars}},\code{\link{plotvarselect}},\code{\link{run_MCMC}}
}

