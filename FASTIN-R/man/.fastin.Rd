\name{fastin}
\alias{fastin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Main analysis wrapper
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fastin(SI.data = NULL, FA.data = NULL, Groups = NULL, Covariates = NULL, eveness = 0.1, nIter = 10000, nBurnin = 5000, Analysis.Type = "Pop.proportions.only")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SI.data}{
 A list with references to filenames containing predator, prey names (1st column) and stable isotope values (remaining columns). Readers are encouraged to run \code{simulation()} to simualte data and inspect the files for the correct format. If fractionation coefficients are not available, a prior mean and sd can be set instread. USers are encouraged to use the gui via \code{run_FASTIN()} rather than using this function directly.
}
  \item{FA.data}{
 A list with references to filenames containing predator, prey names (1st column) and fatty acid profiles (remaining columns). Readers are encouraged to run \code{simulation()} to simualte data and inspect the files for the correct format. If fractionation coefficients are not available, a prior mean and sd can be set instread. USers are encouraged to use the gui via \code{run_FASTIN()} rather than using this function directly.
}
  \item{Groups}{
Optional grouping variable (e.g., regions, sex, season, age classes) that may explain differences in diet proportions
}
  \item{Covariates}{
Optional covariates (e.g., age, length, latitude)
}
  \item{eveness}{
Prior eveness pof diet proportions. A low value will give high prior probablility of specialist diets, a higher value gives a higher prior probablility to a more generalist diet. See gui for numerical examples.
}
  \item{nIter}{
Number of MCMC iterations. Should be relatively low to start in order to verify that the model runs at reasonable speed. If not, fewer fatty acids and/or grouped preys may need to be used. See the gui for guidance when performing these steps.
The number of kept MCMC iterations is allways ~1000 if nIter>1000, so the chain is thinned appropriately for larger numbers of iterations. >10000 iterations are highly recommended to ensure minimal autocorrelation and independent draws for the posterior distribution.
}
  \item{nBurnin}{
Number of MCMC Iterations to discard
}
  \item{Analysis.Type}{
Choose if population proportions or individual and population proportions are to be estimated. The former is aprticularly useful when doing exploratory analysis since it tends to run much faster. The altter is automatically selected (choice is overwritten) if covariates or groups are present.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

An object of class 'post diets', grouping MCMC runs for summaries via \code{summary} and plotting via \code{plot}.

}
\references{
Neubauer.P. and Jensen, O.P. (in prep)
}
\author{
Philipp Neubauer
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{run_FASTIN}},\code{\link{simulation}}
}
\examples{

run_FASTIN()

## The function is currently defined as
function (SI.data = NULL, FA.data = NULL, Groups = NULL, Covariates = NULL, 
    eveness = 0.1, nIter = 10000, nBurnin = 5000, Analysis.Type = "Pop.proportions.only") 
{
    SI = F
    FA = F
    source.combine <- function(k, preys.ix) {
        cat("please select from source combination menu", "\n")
        combine <- menu(title = "combine sources into groups", 
            choices = c("yes", "no"), graphics = T)
        while (combine == 0) {
            cat("please select from source combination menu")
            combine <- menu(title = "combine sources into groups", 
                choices = c("yes", "no"), graphics = T)
        }
        if (combine == 1) {
            selecta <- select.list(preys.names, multiple = T, 
                graphics = T)
            preys.ix[preys.ix \%in\% selecta] <- paste("grouped prey", 
                k)
            preys.ix <- source.combine(k + 1, preys.ix)
        }
        return(preys.ix)
    }
    if (is.list(FA.data)) {
        predators = read.csv(FA.data$predators, header = F, row.names = 1)
        preys = read.csv(FA.data$preys, header = F)
        preys.ix <- as.character(preys[, 1])
        preys.names <- as.character(unique(preys.ix))
        preys = preys[, -1]
        n.fats = ncol(predators)
        m.fats = n.fats - 1
        n.preys <- length(unique(preys.ix))
        if (nchar(FA.data$Conv.Coeffs.mean) > 0 & nchar(FA.data$Conv.Coeffs.var) > 
            0) {
            mean_c = read.csv(FA.data$Conv.Coeffs.mean, header = F, 
                row.names = 1)
            tau_c = 1/read.csv(FA.data$Conv.Coeffs.var, header = F, 
                row.names = 1)^2
        }
        else if (nchar(FA.data$Conv.Coeffs.mean) == 0 & nchar(FA.data$Conv.Coeffs.var) == 
            0) {
            mean_c = matrix(FA.data$CC.mean, n.fats, n.preys)
            tau_c = matrix(1/FA.data$CC.var, n.fats, n.preys)
        }
        else {
            stop("Known conversion coefficients, or a mean AND variance for conversion coefficients need to be supplied")
        }
        if (nchar(FA.data$fat.cont) == 0) {
            fc.mean <- FA.data$FC.mean
            fc.var <- FA.data$FC.var
        }
        else {
            if (dim(fat.cont)[2] > 1) {
                fat.cont <- read.csv(FA.data$fat.cont, header = F, 
                  row.names = 1)
                fc.mean = fat.cont[, 1]
                fc.sd = fat.cont[, 2]
            }
            else {
                fat.cont <- read.csv(FA.data$fat.cont, header = F)
            }
        }
        predators <- t(apply(predators, 1, function(x) {
            x/(sum(x))
        }))
        preys <- t(apply(preys, 1, function(x) {
            x/(sum(x))
        }))
        require("robCompositions")
        dists <- matrix(, nrow(preys), nrow(preys))
        for (i in 1:nrow(preys)) {
            for (j in i:nrow(preys)) {
                dists[j, i] <- aDist(preys[i, ], preys[j, ])
            }
        }
        detach("package:robCompositions", unload = TRUE)
        dista <- as.dist(dists)
        x11()
        PR.RDA <- capscale(dista ~ as.factor(preys.ix), comm = preys)
        plot(preys \%*\% PR.RDA$CCA$v[, 1:2], pch = as.numeric(as.factor(preys.ix)), 
            col = as.numeric(as.factor(preys.ix)))
        points(predators \%*\% PR.RDA$CCA$v[, 1:2], pch = 16)
        legend(locator(2), legend = (unique(preys.ix)), xpd = T, 
            pch = 1:n.preys, col = 2:(n.preys + 1))
        prey.ix <- source.combine(1, preys.ix)
        SC <- T
        n.preys <- length(unique(prey.ix))
        n.preds <- dim(predators)[1]
        x11()
        plot(cumsum(sort(clo(rowSums(t(t(cbind(PR.RDA$CCA$v, 
            PR.RDA$CA$v)) * c(PR.RDA$CCA$eig, PR.RDA$CA$eig))^2)), 
            decreasing = T)), ylab = "cumulative proportion")
        print(cumsum(sort(clo(rowSums(t(t(cbind(PR.RDA$CCA$v, 
            PR.RDA$CA$v)) * c(PR.RDA$CCA$eig, PR.RDA$CA$eig))^2)), 
            decreasing = T)))
        answer <- menu(c("yes", "no"), "Would you like to use a subset of Fatty Acids?", 
            graphics = T)
        while (answer == 0) answer <- menu(c("yes", "no"), "Would you like to use a subset of Fatty Acids?", 
            graphics = T)
        if (answer == 1) {
            sv = sort(clo(rowSums(t(t(cbind(PR.RDA$CCA$v, PR.RDA$CA$v)) * 
                c(PR.RDA$CCA$eig, PR.RDA$CA$eig))^2)), decreasing = T, 
                index.return = T)
            par(ask = T)
            nv <- readline(prompt = "please enter number of variables for analysis \n")
            six <- sv$ix[1:as.numeric(nv)]
            n.fats <- length(six)
            m.fats = n.fats - 1
        }
        else {
            six = 1:n.fats
        }
        mprey <- matrix(, n.preys, n.fats)
        fcc.mean <- rep(NA, n.preys)
        fcc.sd <- rep(NA, n.preys)
        for (i in 1:n.preys) {
            if (is.null(dim(fat.cont))) {
                mprey[i, ] <- apply(preys[prey.ix == unique(prey.ix)[i], 
                  six] * mean_c[match(preys.ix[prey.ix == unique(prey.ix)[i]], 
                  rownames(mean_c)), six], 2, function(x) {
                  weighted.mean(x, w = fat.cont[prey.ix == unique(prey.ix)[i]])
                })
                fcc.mean[i] <- mean(fat.cont[prey.ix == unique(prey.ix)[i]])
                fcc.sd[i] <- var(fat.cont[prey.ix == unique(prey.ix)[i]])
            }
            else {
                mprey[i, ] <- apply(preys[prey.ix == unique(prey.ix)[i], 
                  six] * mean_c[match(preys.ix[prey.ix == unique(prey.ix)[i]], 
                  rownames(mean_c)), six], 2, function(x) {
                  weighted.mean(x, w = fc.mean[match(preys.ix[prey.ix == 
                    unique(prey.ix)[i]], rownames(mean_c))])
                })
                fcc.mean[i] <- mean(fc.mean[rownames(mean_c) \%in\% 
                  preys.ix[prey.ix == unique(prey.ix)[i]]])
                fcc.sd[i] <- mean(fc.sd[rownames(mean_c) \%in\% 
                  preys.ix[prey.ix == unique(prey.ix)[i]]])
            }
        }
        mean_c <- mean_c[1:n.preys, six] * 0 + 1
        preym <- unclass(alr(mprey))
        preds <- unclass(alr(predators[, six]))
        R <- array(, c(m.fats, m.fats, n.preys))
        ni <- rep(NA, n.preys)
        for (i in 1:n.preys) {
            ni[i] <- max(n.fats + 1, sum(prey.ix == unique(prey.ix)[i]) - 
                1)
            R[, , i] = cov(alr(preys[prey.ix == unique(prey.ix)[i], 
                six])) * ni[i]
        }
        Rnot = diag(R.diag, m.fats)
        FA == T
        datas.FA <- list(fc_mean = fcc.mean, fc_tau = 1/fcc.sd^2, 
            n.fats = n.fats, m.fats = m.fats, R = R, Rnot = Rnot, 
            preym = preym, preds = preds, ni = ni, mean_c = mean_c, 
            tau_c = tau_c)
    }
    if (is.list(SI.data)) {
        stopifnot(nchar(SI.data$predators.SI) > 0 & nchar(SI.data$preys.SI) > 
            0)
        predators.SI = read.csv(SI.data$predators.SI, header = F, 
            row.names = 1)
        preys.SI = read.csv(SI.data$preys.SI, header = F)
        preys.ix.SI <- as.character(preys.SI[, 1])
        preys.names.SI <- as.character(unique(preys.ix.SI))
        preys.SI = preys.SI[, -1]
        isos = ncol(predators.SI)
        Rnot_SI = diag(SI.data$R.diag.SI, isos)
        if ((nchar(SI.data$Frac.Coeffs.mean) > 0 & nchar(SI.data$Frac.Coeffs.var) == 
            0) | (nchar(SI.data$Frac.Coeffs.mean) == 0 & nchar(SI.data$Frac.Coeffs.var) > 
            0)) {
            stop("The mean AND variances of FCs for each isotope need to be supplied in the form c(FC1,FC2)")
        }
        else if (nchar(SI.data$Frac.Coeffs.mean) > 0 & nchar(SI.data$Frac.Coeffs.var) > 
            0) {
            mean_cs = read.csv(SI.data$Frac.Coeffs.mean, header = F, 
                row.names = 1)
            sd_cs = 1/read.csv(SI.data$Frac.Coeffs.var, header = F, 
                row.names = 1)^2
        }
        else if (nchar(SI.data$Frac.Coeffs.mean) == 0 & nchar(SI.data$Frac.Coeffs.var) == 
            0) {
            mean_cs = matrix(SI.data$FC.mean, isos, n.preys.SI)
            tau_cs = matrix(1/SI.data$FC.var, isos, n.preys.SI)
        }
        if (SC == F) {
            if (isos <= 2) {
                plot(preys.SI[, 1:ncol(preys.SI)], pch = as.numeric(as.factor(preys.ix.SI)), 
                  col = as.numeric(as.factor(preys.ix.SI)) + 
                    1)
                points(predators.SI[, 1:2], pch = 16)
                tkmessageBox(message = "please use the cursor to select lower right and upper left corner for legend", 
                  title = "Plot Legend")
                legend(locator(2), legend = preys.names.SI, xpd = T, 
                  pch = 1:n.preys, col = 2:(n.preys + 1))
            }
            else {
                PR.RDA <- capscale(preys.SI ~ as.factor(preys.ix.SI))
                plot(data.matrix(preys.SI) \%*\% PR.RDA$CCA$v[, 
                  1:2], pch = as.numeric(as.factor(preys.ix.SI)), 
                  col = as.numeric(as.factor(preys.ix.SI)) + 
                    1)
                points(data.matrix(predators.SI) \%*\% PR.RDA$CCA$v[, 
                  1:2], pch = 16)
                tkmessageBox(message = "please use the cursor to select lower right and upper left corner for legend", 
                  title = "Plot Legend")
                legend(locator(2), legend = preys.names.SI, xpd = T, 
                  pch = 1:n.preys, col = 2:(n.preys + 1))
            }
            prey.ix <- source.combine(1, preys.ix.SI)
            n.preys <- length(unique(prey.ix))
            n.preds <- dim(predators.SI)[1]
        }
        preym.SI <- matrix(, n.preys, isos)
        for (i in 1:n.preys) {
            preym.SI[i, ] <- apply(preys.SI[prey.ix == unique(prey.ix)[i], 
                ] + mean_cs[match(preys.ix.SI[prey.ix == unique(prey.ix)[i]], 
                rownames(mean_cs)), ], 2, mean)
        }
        mean_cs <- mean_cs[1:n.preys, ] * 0
        R.SI <- array(, c(isos, isos, n.preys))
        ni.SI <- rep(NA, n.preys)
        for (i in 1:n.preys) {
            ni.SI[i] <- max(isos + 1, sum(prey.ix == unique(prey.ix)[i]) - 
                1)
            R.SI[, , i] = cov(preys.SI[prey.ix == unique(prey.ix)[i], 
                ]) * ni.SI[i]
        }
        SI == T
        datas.SI <- list(isos = isos, R.SI = R.SI, Rnot.SI = Rnot.SI, 
            preym.SI = preym.SI, preds.SI = preds.SI, ni.SI = ni.SI, 
            mean_cs = mean_cs, tau_cs = tau_cs)
    }
    datas <- list(n.preys = n.preys, n.preds = n.preds, datas.FA = NULL, 
        datas.SI = NULL, eveness = eveness)
    if (SI == T & FA == T) {
        datas$datas.FA <- datas.FA
        datas$datas.SI <- datas.SI
        class(datas) <- "combined"
    }
    else if (SI == F & FA == T) {
        datas$datas.FA <- datas.FA
        class(datas) <- "FA"
    }
    else if (SI == T & FA == F) {
        datas$datas.SI <- datas.SI
        class(datas) <- "SI"
    }
    else {
        stop("At least one type of data (Fatty acids and/or Stable Isotopes) must be supplied")
    }
    outputs <- switch(Analysis.Type, Pop.proportions.only = Pop.prop.analysis(datas, 
        nIter, nBurnin), Pop.and.Ind.proportions = Pop.and.Ind.props(datas, 
        nIter, nBurnin), Grouped.Analysis.Linear.Model = Grouped.analysis(datas, 
        nIter, nBurnin, Groups = Groups), Analysis.with.Covariates = Analysis.with.Cov(datas, 
        nIter, nBurnin, Covariates = Covariates), Grouped.Analysis.with.Covariates = Grouped.Analysis.with.Cov(datas, 
        nIter, nBurnin, Covariates = Covariates, Groups = Groups))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{  manip}
