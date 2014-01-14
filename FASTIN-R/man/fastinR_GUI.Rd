\name{fastinR_GUI}
\alias{fastinR_GUI}
\title{GUI for the fastinR set of functions}
\usage{
fastinR_GUI()
}
\description{
Provides a gui to input data to be analysed with the
\code{\link{run_MCMC}} function. The gui rolls the
functionality of most of the functions contained within
fastinR into a single interface, but some options (plottin
options mainly) are not available through the gui.
}
\details{
Please consult the help individual functions to obtain help
on a particular topic
}
\section{Warning}{
  The Tcl/Tk gui interface is very unpredictable, leading
  to odd errors like internal functions not being found
  when calling FASTIN functions from the console after
  having used the GUI. Problems seem to be related to
  environments. It's usually best to completely restart R
  when this happens.
}
\examples{
\dontrun{fastin_GUI()}
}
\author{
Philipp Neubauer
}
\references{
Neubauer.P. and Jensen, O.P. (in prep)
}
\seealso{
\code{\link{simulation}},\code{\link{add_FA}},\code{\link{add_SI}},\code{\link{run_MCMC}}
}

