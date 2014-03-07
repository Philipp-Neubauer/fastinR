\name{select_vars}
\alias{select_vars}
\title{Select Fatty Acids from imported data (should really be a method for subset...)}
\usage{
select_vars(datas, ix = NULL, plot = T)
}
\arguments{
  \item{datas}{A data structure produced by
  \code{\link{add_SI}} and \code{\link{add_FA}}}

  \item{ix}{Index of variables to retain, if not supplied,
  \code{\link{var_select_plot}} is called to itneractively
  select variables based on cosine distance and prey matrix
  condition.}

  \item{plot}{Logical (T/F): should selection be
  itneractive? Needs tcl-tk package(and X-quartz on Mac)}
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
select_vars(datas,4:1)
# or select visualy
select_vars(datas)
}
}
\author{
Philipp Neubauer
}
\references{
Neubauer.P. and Jensen, O.P. (in prep)
}
\seealso{
\code{\link{add_FA}},\code{\link{var_select_plot}},\code{\link{run_MCMC}},\code{\link{simulation}}
}

