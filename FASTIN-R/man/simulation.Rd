\name{simulation}
\alias{simulation}
\title{GUI: Simulate Fatty Acid profiles and Stable Isotope data}
\usage{
simulation()
}
\value{
\code{Write simulations to files} will produce a series of
files that can be used as inputs to the FASTIN gui or
individual functions. The prefix of the files is entered by
the user, the remainder of the filename suggests the
contents: \item{*_props}{Simualted diet proportions}

For Fatty Acids:

\item{*_FA_preys.csv}{Prey names (first column) and fatty
acid profiles, with fatty acids names across the first row}
\item{*_FA_preds.csv}{Predator index (first column) and
fatty acid profiles, with fatty acids names across the
first row} \item{*_fat_cont.csv}{Prey fat contents (here as
mean and variance, can be specified for each prey sample
for the main analysis, in that case the first column is the
prey id and the second column is the individual sample's
fat content)} \item{*_FA_cc_means.csv}{Prey specific
conversion coefficient means: Prey names (first column) and
an n x P matrix for n preys and P fatty acids}
\item{*_FA_cc_var.csv}{Prey specific conversion coefficient
variances, dimensions as for the means} For Stable
Isotopes: \item{*_SI_preys.csv}{Prey names (first column)
and stable isotope values, with SI names across the first
row} \item{*_SI_preds.csv}{Predator index (first column)
and stable isotope values, with SI names across the first
row} \item{*_SI_fc_means.csv}{Prey specific SI additive
fractionation coefficient means}
\item{*_SI_fc_var.csv}{Prey specific additive fractionation
coefficient variance, dimensions as for the means} For
Covariates for predator proportions and grouped predators:

\item{*_Cov_n_Grp_effects.csv}{an n*p matrix of group and
covariate influences on diet proportions for n preys and p
(groups+covariates)} \item{*_Groups.csv}{Group membership
for each predator} \item{*_Covariates.csv}{Covariate values
for each predator}
}
\description{
GUI: Simulate Fatty Acid profiles and Stable Isotope data
}
\details{
The function calls a gui to simulate data, functions are of
little use outside the gui facility.

Simulating data proceeds by selecting to simualte either
individual diet proportions for each predator as drawn from
a population level distribution for diet proportions, or by
simulating grouped diet proportions (e.g., size groups,
geographic groups) and/or by letting covariates influence
the drawn diet proportions (e.g., size based diet shifts).

Pressing \code{Simulate Fatty Acid data} or \code{Simulate
Stable Isotope data} makes sense only after proportions
were simulated - converseley, when updating proportions, FA
and/or SI data need to be re-simulated to have the new diet
proportions alter the predator fatty acid/SI makeup.

\code{Plot current simulation} will draw a NMDS (non-metric
multidimensional scaling) plot to show the simulated data.

The gui will remain open to allow for many tries at
simulating data. Variuous configurations can be written to
file to allow for exploration of the model with
\code{fastin()}.
}
\author{
Philipp Neubauer
}
\references{
Neubauer.P. and Jensen, O.P. (in prep)
}

