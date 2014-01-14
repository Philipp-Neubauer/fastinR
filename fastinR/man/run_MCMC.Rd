\name{run_MCMC}
\alias{run_MCMC}
\title{MCMC for diet proportions from FAtty acid and/or Stable Isotope data}
\usage{
run_MCMC(datas = NULL, Covs = NULL, nIter = 10000, nBurnin = 1000,
  nChains = 1, nThin = 10, Data.Type = "Fatty.Acid.Profiles",
  Analysis.Type = "Population.proportions", even = 0.5, Rnot = 0.2,
  Rnot_SI = 1, plott = T, spawn = F)
}
\arguments{
  \item{datas}{A data object as produced by
  \code{\link{add_FA}} and/or \code{\link{add_SI}}. While
  it is possible to construct the data objects (see
  examples on github), it is advisable to use the above
  mentioned functions for convenience.}

  \item{Covs}{Covariates to use in a lienar mixed setup, as
  produced by \code{\link{add_Covs}}}

  \item{nIter}{number of MCMC iterations post- burn-in}

  \item{nBurnin}{number of iterations to discard as burn-in
  of the Markov Chain}

  \item{nChains}{The number of MCMC chains to run
  simultaniously.}

  \item{nThin}{Thinning interval of the Markov Chain}

  \item{Data.Type}{The data to be analysed: one of
  'Fatty.Acid.Profiles', 'Stable.Isotopes' or
  'Combined.Analysis' for estimation of diet proportions
  from a joint model.}

  \item{Analysis.Type}{The type of analysis to perform: one
  of'Population.proportions' for estimation of population
  level proportions only, 'Individual.proportions' to
  estimate individual diet proportions, or
  'Analysis.with.Covariates' for individual proportions and
  estiamtes of group contrasts and/or covariate effects.}

  \item{even}{The prior eveness of diet proportions,
  smaller values place a lower prior on even proportions,
  if set too low can severly impact convergence of MCMC.
  Too high a value will set a strong prior on equal
  proportions. Only needed for individual proportions or an
  analysis with covariates, else a vague prior is used and
  proportions are drawn from a Dirichlet conditional
  posterior.}

  \item{Rnot}{diagnoal of the prior for the predator stable
  isotope covariance matrix. Note that this parameter can
  significantly influence convergence, especially if there
  are few predator signatures to estiamte the covariance
  from - handle with care! If too small, the sampler will
  get stuck in local modes and extreme values, if too high
  it can produce nonsenseical estiamtes where all
  proportions are equal (posterior mean at
  \code{1/n.preys})}

  \item{Rnot_SI}{diagnoal of the prior for the predator
  fatty acid covariance matrix. Note that this parameter
  can significantly influence convergence, especially if
  there are few predator signatures to estiamte the
  covariance from - handle with care! If too small, the
  sampler will get stuck in local modes and extreme values,
  if too high it can produce nonsenseical estiamtes where
  all proportions are equal (posterior mean at
  \code{1/n.preys})}

  \item{plott}{If FALSE, user is not asked if MCMC runs
  should be displayed using plots from the coda package. By
  default, the user is prompted (for compatibility with the
  gui).}

  \item{spawn}{Logical - should separate slave R processes
  be spawned to run individual chains? This can provide
  significant speed-up for long runs using multiple chains
  on multi-core processors. A suggested use is to set
  \code{spawn=F} for short exploratory MCMC runs and once a
  satisfactory set of parameters has been found, set
  \code{spawn=T} and run multiple chains for longer.}
}
\value{
An object containing the MCMC runs, where the class
corresponds to the analysis type to enable methods dispatch
for plots and summaries.
}
\description{
MCMC for diet proportions from FAtty acid and/or Stable
Isotope data
}
\details{
This analysis can be run from the gui or from this
standalone R function. In the latter case it is recommended
to set plott = F to avoid plotting errors to prevent a
return of the fucntion (thereby loosing potential long
runs). Markov Chains can always be plotted post-hoc with
,\code{\link{MCMCplot}}.
}
\examples{
\dontrun{
data('Sim')
# a short analysis...
MCMCout <- run_MCMC(datas=datas,nIter=100,nBurnin=100,nChains=1,nThin=1,Data.Type='Fatty.Acid.Profiles',Analysis.Type='Population.proportions',even=0.5,plott=F,spawn=F)

# summaries
summary(MCMCout)

# diagnostics
diags(MCMCout)

# plot chains
MCMCplot(MCMCout)
plot(MCMCout,save=F)
}
}
\author{
Philipp Neubauer
}
\references{
Neubauer,.P. and Jensen, O.P. (in prep)
}
\seealso{
\code{\link{add_FA}},\code{\link{add_SI}},\code{\link{diags}},\code{\link{add_Covs}}
,\code{\link{plot}}
}

