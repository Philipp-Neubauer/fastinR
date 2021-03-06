% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/Grooming.R
\name{print_var_list}
\alias{print_var_list}
\title{Prints cosine distance and matrix condition number for variable selection}
\usage{
print_var_list(prey.mat, prey.ix)
}
\arguments{
\item{prey.mat}{(row-wise) dataframe of prey compositions}

\item{prey.ix}{Indexes rows of \code{prey.mat}}
}
\value{
Prints a matrix in the console with columns being 1) relative contributions of Fatty Acids to source separation, 2) matrix condition number and 3) original index used to select the fatty acids.
}
\description{
Cosine distance based separation on CAP axes and matrix condition number are printed as a function of Fatty Acids, where fatty acids are ordered by their relative contribution to Canonical axes in a CAP (Canonical Analysis of Principal Coordinates). The original indices of these fatty acids in the data are given to facilitate selection
}
\examples{
\dontrun{
# load simulated example
data('Sim')
preys <- datas$datas.FA$preys
prey.ix <- datas$prey.ix
print_var_list(preys,prey.ix)
}
}
\author{
Philipp Neubauer
}
\references{
Neubauer.P. and Jensen, O.P. (in prep)
#'
}
\seealso{
\code{\link{add_FA}},\code{\link{select_vars}},\code{\link{run_MCMC}},\code{\link{simulation}}
}

