\name{dataplot}
\alias{dataplot}
\title{Plot Non-metric multidimensional scaling (NMDS) or (ilr transformed) diet data}
\usage{
dataplot(datas = NULL)
}
\arguments{
  \item{datas}{Diet data produced by \code{\link{add_FA}}
  and/or \code{\link{add_SI}}}
}
\description{
Diet data produced by \code{\link{add_FA}} and/or
\code{\link{add_FA}} is projected onto 2 dimesions using
NMDS, where Fatty Acid data is first ilr transformed. When
SI and FA data are present, then Euclidean distances of
ilr(FA) and SI data are added to produce NMDS and the plot.
}
\author{
Philipp Neubauer
}
\references{
Neubauer,.P. and Jensen, O.P. (in prep)
}
\seealso{
\code{\link{add_FA}},\code{\link{add_SI}},\code{\link{select_vars}},\code{\link{var_select_plot}},\code{\link{run_MCMC}}
}

