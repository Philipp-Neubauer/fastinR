\name{run_FASTIN}
\alias{run_FASTIN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
run_FASTIN()
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function () 
{
    add.SI <- function(predators.SI = NULL, preys.SI = NULL, 
        Frac.Coeffs.mean = NULL, Frac.Coeffs.var = NULL, FC.mean = 1, 
        FC.var = 1, R.diag.SI = 1e-10) {
    }
    add.FA <- function(predators, preys, fat.cont = NULL, Conv.Coeffs.mean = NULL, 
        Conv.Coeffs.var = NULL, FC.mean = 1, FC.var = 1, CC.mean = 1, 
        CC.var = 1, R.diag = 1e-10) {
    }
    pnorm_even <- function(eveness = 0.1) {
        p = 2 * (1 - pnorm(log(95)/2, 0, sqrt(1/eveness)))
        return(p)
    }
    guiSet("LIST_WIDTH", 50)
    guiSet("ENTRY_WIDTH", 10)
    output <- gui(FASTIN, argCommand = list(eveness = guiNestedF(pnorm_even, 
        "eveness", exec = "Set new prior", argSlider = list(eveness = c(0, 
            2, 0.01), output = c(0, 1, 0.01)), default = 0.1, 
        callback = guiExec, argText = c(output = "Probability that prop(prey x) = 0.95*prop(preys other than x)")), 
        SI.data = guiNestedF(add.SI, "SI.data", argFilter = list(predators.SI = "{{} {.csv}}", 
            preys.SI = "{{} {.csv}}", Frac.Coeffs.mean = "{{} {.csv}}", 
            Frac.Coeffs.var = "{{} {.csv}}"), title = "Stable Isotope data entry form", 
            exec = "Add Stable Isotope data", closeOnExec = TRUE, 
            output = NULL, argText = c(predators.SI = "Load predator(s) stable isotope data (csv)", 
                preys.SI = "Load prey stable isotope data (csv)", 
                Frac.Coeffs.mean = "Load prey and/or stable isotope specific fractionation means (csv)", 
                Frac.Coeffs.var = "Load prey and/or stable isotope specific fractionation variances (csv)", 
                FC.mean = "Mean fractionation coefficients (use R's c() notation)", 
                FC.var = "Variance of fractionation coefficients (use R's c() notation)", 
                R.diag.SI = "Diagonal of the prior for predator (co)-variance matrix (>0, smaller value is less informative)"), 
            cancelButton = F), FA.data = guiNestedF(add.FA, "FA.data", 
            argFilter = list(predators = "{{} {.csv}}", preys = "{{} {.csv}}", 
                fat.cont = "{{} {.csv}}", Conv.Coeffs.mean = "{{} {.csv}}", 
                Conv.Coeffs.var = "{{} {.csv}}"), title = "Fatty Acid Profile data entry form", 
            argText = c(predators = "Load predator(s) fatty acid data (csv)", 
                preys = "Load prey fatty acid data (csv)", Conv.Coeffs.mean = "Load prey and/or fatty acid specific conversion means (csv)", 
                Conv.Coeffs.var = "Load prey and/or fatty acid specific conversion variances (csv)", 
                Fat.Cont = "Prey fat content (csv)", FC.mean = "Mean fat content (use R's c() notation)", 
                FC.var = "Variance of fat content (use R's c() notation)", 
                CC.mean = "Mean conversion coefficients (use R's c() notation)", 
                CC.var = "Variance of conversion coefficients (use R's c() notation)", 
                R.diag = "Diagonal of the prior for predator (co)-variance matrix (>0, smaller value is less informative)"), 
            exec = "Add Fatty Acid data", cancelButton = F, closeOnExec = TRUE, 
            output = NULL)), argFilter = list(Groups = "{{} {.csv}}", 
        Covariates = "{{} {.csv}}"), argOption = list(Analysis.Type = c("Population.proportions", 
        "Population.and.Individual.proportions"), defaultChoice = 1), 
        argText = c(SI.data = "Add Stable Isotope data", FA.data = "Add Fatty Acid data", 
            Analysis.Type = "Choose Model Setup", Covariates = "Add Covariates (optional)", 
            Groups = "Add Groups (optional)"), helps = list(Analysis.Type = "Option 1 (Population proportions only) is faster since a Dirichlet prior is used, but Option 2 (Pop. and Individual proportions) is run by default for Covariates and Groups"), 
        closeOnExec = TRUE, output = NULL, argGridOrder = c(1, 
            1, 2, 2, 3, 3, 3, 4, 4), argGridSticky = rep("w", 
            length(formals(FASTIN))))
    detach("package:fgui", unload = TRUE)
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
