\name{fastin}
\alias{fastin}
\title{Gui for the FASTIN analysis}
\usage{
fastin()
}
\description{
Provides a gui to input data to be analysed with the
\code{\link{run_MCMC}} function. The gui rolls the
functionality of most of the functions contained within
FASTIN into a single interface, but some options (plottin
options mainly) are not available through the gui.
}
\details{
Please consult the help individual functions to obtain help
on a particular topic
}
\section{Warning the Tcl/Tk gui interface is very unpredictable, leading to odd errors like internal functions not being found when calling FASTIN functions from the console after having used the GUI. Problems seem to be related to environments. It's usually best to completely restart R when this happens.}{

}
\author{
Philipp Neubauer
}
\references{
Neubauer.P. and Jensen, O.P. (in prep)
}
\seealso{
\code{\link{addFA}},\code{\link{addSI}},\code{\link{run_MCMC}},\code{\link{diags}}
}

