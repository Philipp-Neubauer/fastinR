\name{Pop.prop.analysis.SI}
\alias{Pop.prop.analysis.SI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Pop.prop.analysis.SI(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
}
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
function (data) 
{
    n.preys = data$n.preys
    n.preds = data$n.preds
    R.SI = data$datas.SI$R.SI
    mean_cs = data$datas.SI$mean_cs
    tau_cs = data$datas.SI$tau_cs
    Rnot.SI = data$datas.SI$Rnot.SI
    isos = data$datas.SI$isos
    ni.SI = data$datas.SI$ni.SI
    preds.SI = data$datas.SI$preds.SI
    preym.SI = data$datas.SI$preym.SI
    initials.SI = list(list(ps = rep(1/n.preys, n.preys), cs = matrix(mean_cs, 
        n.preys, isos, byrow = T), prey.means_SI = preym.SI, 
        predprec_SI = diag(0.01, isos), prey.precs_SI = array(1, 
            c(n.preys, isos, isos))))
    datas.SI = list("R_SI", "mean_cs", "tau_cs", "Rnot_SI", "n.preys", 
        "n.preds", "isos", "ni.SI", "preds.SI", "preym.SI")
    vars = c("prop")
    res <- BRugsFit("Pop.prop.analysis.SI.bugs", datas.SI, inits = initials.SI, 
        numChains = 1, vars, nBurnin = nBurnin, nIter = Iter, 
        nThin = nIter/1000, coda = T, DIC = F, working.directory = getwd(), 
        digits = 4, BRugsVerbose = T)
    res <- as.data.frame((res)[[1]])
    class(res) <- "pop.props"
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
