\name{resetsc}
\alias{resetsc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Reset prey groupings
}
\description{Reset prey groupings as specified in  \code{\link{addFA}} and/or
  \code{\link{addSI}}
}
\usage{
resetsc(datas=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{datas}{A data object as produced by \code{\link{addFA}} and/or
  \code{\link{addSI}}. While it is possible to construct the data
  objects (see examples on github), it is advisable to use the above
  mentioned functions for convenience. IF no data object is supplied,
  the function will look for one in the gui namespace}

}
\details{
This function can be run from the gui or from this standalone R
  function. 
}
\value{
 Resets prey groupings - requires a new call to \code{\link{addFA}} or
 \code{\link{addSI}} to add data with the original/new prey groupings.
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