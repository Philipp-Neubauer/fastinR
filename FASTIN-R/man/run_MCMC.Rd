\name{run_MCMC}
\alias{run_MCMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MCMC for diet proportions from FAtty acid and/or Stable Isotope data
%%  ~~function to do ... ~~
}
\description{Run an Bayesian analysis using MCMC for estimation of diet
  proportions from Fatty acid data (entered via the \code{\link{FASTIN}}
  gui or \code{\link{addFA}}) and/or Stable Isotope (entered via the \code{\link{FASTIN}}
  gui or \code{\link{addSI}}).
}
\usage{
run_MCMC(datas=NULL,nIter=10000,nBurnin=1000,nChains=1,nThin=10,Data.Type='Fatty.Acid.Profiles',Analysis.Type='Population.proportions',even=0.1,plott=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{datas}{A data object as produced by \code{\link{addFA}} and/or
  \code{\link{addSI}}. While it is possible to construct the data
  objects (see examples on github), it is advisable to use the above
  mentioned functions for convenience.}
\item{nIter}{number of MCMC iterations post- burn-in}
\item{nBurning}{number of iterations to discard as burn-in of the Markov
  Chain}
\item{nChains}{The number of MCMC chains to run simultaniously.}
\item{nThin}{Thinning interval of the Markov Chain}
\item{Data.Type}{The data to be analysed: one of
  'Fatty.Acid.Profiles', 'Stable.Isotopes' or 'Combined.Analysis' for
  estimation of diet proportions from a joint model.}
\item{Analysis.Type}{The type of analysis to perform: one of
  'Population.proportions' for estimation of population level
  proportions only, 'Individual.proportions' to estimate individual diet
  proportions, or 'Analysis.with.Covariates' for individual proportions
  and estiamtes of group contrasts and/or covariate effects.}
\item{even}{The prior eveness of diet proportions, smaller values place
  a lower prior on even proportions}
\item{plott}{If FALSE, user is not asked if MCMC runs should be
  displayed using plots from the coda package. By default, the user is
  prompted (for compatibility with the gui).}
}
\details{
This analysis can be run from the gui or from this standalone R
  function. 
}
\value{
  An object containing the MCMC runs, where the class corresponds to the
  analysis type to enable methods dispatch for plots and summaries.
}

\references{

Neubauer.P. and Jensen, O.P. (in prep)

}
\author{

Philipp Neubauer

}
\seealso{
\code{\link{addFA}},\code{\link{addSI}},\code{\link{diags}},\code{\link{addcovs}},\code{\link{plot.pop_props}}
}