#' MCMC for diet proportions from FAtty acid and/or Stable Isotope data
#' 
#' @param datas A data object as produced by \code{\link{add_FA}} and/or \code{\link{add_SI}}. While it is possible to construct the data objects (see examples on github), it is advisable to use the above mentioned functions for convenience.
#' @param Covs Covariates to use in a lienar mixed setup, as produced by \code{\link{add_Covs}}
#' @param nIter number of MCMC iterations post- burn-in
#' @param nBurnin number of iterations to discard as burn-in of the Markov Chain
#' @param nChains The number of MCMC chains to run simultaniously.
#' @param nThin Thinning interval of the Markov Chain
#' @param Data.Type The data to be analysed: one of 'Fatty.Acid.Profiles', 'Stable.Isotopes' or 'Combined.Analysis' for estimation of diet proportions from a joint model.
#' @param Analysis.Type The type of analysis to perform: one of'Population.proportions' for estimation of population level proportions only, 'Individual.proportions' to estimate individual diet proportions, or 'Analysis.with.Covariates' for individual proportions and estiamtes of group contrasts and/or covariate effects.
#' @param even The prior eveness of diet proportions, smaller values place a lower prior on even proportions, if set too low can severly impact convergence of MCMC. Too high a value will set a strong prior on equal proportions. Only needed for individual proportions or an analysis with covariates, else a vague prior is used and proportions are drawn from a Dirichlet conditional posterior.
#' @param Rnot diagnoal of the prior for the predator stable isotope covariance matrix. Note that this parameter can significantly influence convergence, especially if there are few predator signatures to estiamte the covariance from - handle with care! If too small, the sampler will get stuck in local modes and extreme values, if too high it can produce nonsenseical estiamtes where all proportions are equal (posterior mean at \code{1/n.preys})
#' @param Rnot_SI diagnoal of the prior for the predator fatty acid covariance matrix. Note that this parameter can significantly influence convergence, especially if there are few predator signatures to estiamte the covariance from - handle with care! If too small, the sampler will get stuck in local modes and extreme values, if too high it can produce nonsenseical estiamtes where all proportions are equal (posterior mean at \code{1/n.preys})
#' @param plott If FALSE, user is not asked if MCMC runs should be displayed using plots from the coda package. By default, the user is prompted (for compatibility with the gui).
#' @param spawn Logical - should separate slave R processes be spawned to run individual chains? This can provide significant speed-up for long runs using multiple chains on multi-core processors. A suggested use is to set \code{spawn=F} for short exploratory MCMC runs and once a satisfactory set of parameters has been found, set \code{spawn=T} and run multiple chains for longer.
#' @return An object containing the MCMC runs, where the class corresponds to the analysis type to enable methods dispatch for plots and summaries.
#' @details This analysis can be run from the gui or from this standalone R function. In the latter case it is recommended to set plott = F to avoid plotting errors to prevent a return of the fucntion (thereby loosing potential long runs). Markov Chains can always be plotted post-hoc with ,\code{\link{MCMCplot}}.
#' @seealso \code{\link{add_FA}},\code{\link{add_SI}},\code{\link{diags}},\code{\link{add_Covs}} ,\code{\link{plot}}
#' @author Philipp Neubauer
#' @references Neubauer,.P. and Jensen, O.P. (in prep)
#' @examples \dontrun{
#' data('Sim')
#' # a short analysis...
#' MCMCout <- run_MCMC(datas=datas,nIter=100,nBurnin=100,nChains=1,nThin=1,Data.Type='Fatty.Acid.Profiles',Analysis.Type='Population.proportions',even=0.5,plott=F,spawn=F)
#' 
#' # summaries
#' summary(MCMCout)
#' 
#' # diagnostics
#' diags(MCMCout)
#' 
#' # plot chains
#' MCMCplot(MCMCout)
#' plot(MCMCout,save=F)
#' }
#' @export
run_MCMC <- function(datas=NULL,Covs=NULL,nIter=10000,nBurnin=1000,nChains=1,nThin=10,Data.Type='Fatty.Acid.Profiles',Analysis.Type='Population.proportions',even=0.5,Rnot=0.2,Rnot_SI=1,plott=T,spawn=F){
  # have three types here: FA, SI and combined, then methods dispatch based on type of arg
  
  if(missing(datas)) {datas = guiGetSafe('datas');GUI = T} else {GUI = F}
  
  if(length(datas)<=1){warning('No data processed yet')} else {   
    
    datas$even <- even
    
    # Assign data type
    if (Data.Type == 'Combined.Analysis')
    {
      datas$datas.SI$Rnot.SI <- diag(Rnot_SI,datas$datas.SI$isos)
      datas$datas.FA$Rnot <- diag(Rnot,datas$datas.FA$m.fats)
      datas$datas.FA$preys <- alr(datas$datas.FA$preys)
      class(datas) <- 'combined'
    } else if (Data.Type=='Fatty.Acid.Profiles')
    {
      datas$datas.FA$Rnot <- diag(Rnot,datas$datas.FA$m.fats)
      class(datas) <- 'FA'
      datas$datas.FA$preys <- alr(datas$datas.FA$preys)
    } else if (Data.Type=='Stable.Isotopes')
    {
      datas$datas.SI$Rnot.SI <- diag(Rnot_SI,datas$datas.SI$isos)
      class(datas) <- 'SI'
    }
    
    # check for covariates and delist data
    if(Analysis.Type == 'Analysis.with.Covariates') {
      if(GUI ==T) Covs = guiGetSafe('Covs')
      if(any(is.na(Covs)) & Analysis.Type == 'Analysis.with.Covariates') {
        stop('analysis with covariates selected, but no covariates entered.')
      } else {
        jagsdata <- .delist(datas,Covs)
      }
    } else {
      jagsdata <- .delist(datas)
    }   
    
    save(jagsdata,file='jagsdata.Rdata')
    
    sysfile <- switch(Analysis.Type,
                      Population.proportions = .Poppropanalysis(datas),
                      Individual.proportions = .PopandIndprops(datas),
                      Analysis.with.Covariates = .AnalysiswithCov(datas)
    )
    
    if(spawn==T | spawn == 1) {
      res <- .spawn(sysfile,nChains,nBurnin,nIter,nThin)
    } else {
      res <- .localrun(jagsdata,sysfile,nChains,nBurnin,nIter,nThin)
    }
    
    
    if (plott==T){
      plotta <- menu(title='plot MCMC chains?',choices = c('yes','no'),graphics=T)
      if (plotta==1) {
        if (!is.function(options()$device)){
          if (names(dev.cur())=="RStudioGD"){
            # try to open a new platform-appropriate plot window
            if (.Platform$OS.type=='windows'){
              windows()
            } else if(length(grep(R.version$platform,pattern='apple'))>0)  # is it mac?
            { 
              quartz(width=5,height=5)
            } else {  # must be unix
              x11()
            }
          }
        }
        
        try(plot(res,ask=T))
      }
    }
    
    res$res <- c(res$res,nChains = nChains)
    res$res <- c(res$res,prey.names = list(unique(datas$prey.ix)))
    if (Analysis.Type == 'Analysis.with.Covariates') res <- c(res,Covs=list(Covs))
    
    class(res$res) <- switch(Analysis.Type,
                         Population.proportions = 'pop_props',
                         Individual.proportions = 'ind_props',
                         Analysis.with.Covariates = 'cov_props'
    )
    
    ifelse(GUI==F,return(res),guiSet('MCMCout',res))
  }
}

# This is called internally in the slave processes to run jags
.jagger <- function(sysfile,nBurnin,nIter,nThin,i){
  options(warn=-1)
  load('jagsdata.Rdata')
  JM <- jags.model(file=sysfile,data=jagsdata,inits=list(.RNG.name="base::Mersenne-Twister",.RNG.seed=i),n.chains=1)
  
  cat('\n','proceeding to burn-in phase','\n')
  update(JM,n.iter=nBurnin)
  cat('\n','sampling from parameter distributions','\n')
  if(length(grep('Ind',sysfile))>0){
    res<- coda.samples(model=JM,variable.names=c('prop','pop.prop'),n.iter=nIter,thin=nThin)
  } else if(length(grep('Cov',sysfile))>0) {
    res<- coda.samples(model=JM,variable.names=c('prop','pop.prop','beta'),n.iter=nIter,thin=nThin)    
  } else {
    res<- coda.samples(model=JM,variable.names='prop',n.iter=nIter,thin=nThin)
  }
  return(res)
  
}

#spawn slave processes to run separate chains - ultimately this should integrate with runjags
.spawn <- function(sysfile,nChains,nBurnin,nIter,nThin){
  
  orgtime <- Sys.time()
  options(useFancyQuotes=F)
  # spawn seperate R slave processes for each chain and run jagger in them
  for (i in 1:nChains){
    outfile <- paste('MCout',i,'.Rdata',sep='')
    
    cmd <- paste("Rscript -e 'options(warn=-1);library(fastinR,quietly=T,verbose=F);res<- .jagger(",dQuote(sysfile),",",eval(nBurnin),",",eval(nIter),",",eval(nThin),",",i,");save(res,file=",dQuote(outfile),");print(",dQuote("all done"),")'",sep='')
    
    system(cmd,wait=F)
    
  }
  
  while(!any(grep('MCout',dir()))){Sys.sleep(3);cat('+')}
  while(length(grep('MCout',dir()))<nChains | sum(file.info(dir()[grep('MCout',dir())])$mtime >  orgtime)<nChains) {Sys.sleep(3);cat('+')}
  
  
  results <- list()
  
  if (nChains>=1) {
    for(i in 1:nChains){
      load(paste('MCout',i,'.Rdata',sep=''))
      results[[i]] <- res[[1]]
    }
  }
  
  return(results)
  
}

# run MCMC locally
#' @import rstan
.localrun <- function(jagsdata,sysfile,nChains,nBurnin,nIter,nThin){
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  if(length(grep('Ind',sysfile))>0){
    variable.names=c('prop','pop.prop')
  } else if(length(grep('Cov',sysfile))>0) {
    variable.names=c('prop','pop.prop','beta')   
  } else {
    variable.names='prop'
  }
  #browser()
  #stanmodels[[sysfile]]
  JM <- sampling(stanmodels[[sysfile]],
                 init = function() list(props = rep(1/jagsdata$n_preys, jagsdata$n_preys),
                                        prey_means = jagsdata$preym,
                                        cons_prey  = jagsdata$preym,
                                        prey_means_SI = jagsdata$preym_SI,
                                        cons_prey_SI  = jagsdata$preym_SI,
                                        cc = jagsdata$mean_cs,
                                        cs = jagsdata$mean_c,
                                        fc = exp(jagsdata$fc_mean)
                                        ),
                 chains = nChains, 
                 warmup = nBurnin,
                 thin =nThin,
                 iter =nIter,
                 data=jagsdata)
  res <- As.mcmc.list(JM, pars=variable.names)
  
  
  list(res=res, mod=JM)
}

#' @title Raftery-Lewis and Gelman-Rubin diagnostics for MCMC chains from \code{\link{run_MCMC}}. Multiple chains are needed for the latter type of diagnostic.
#' 
#' @param MCMCout An object produced by \code{\link{run_MCMC}}
#' @param accuracy Accuracy with which parameters are to be estiamted
#' @param proba probability with which estiamtes are within the itnerval [quant,1-quant]
#' @param quant interval within which to estiamte the parameter
#' @details The function allows to diagnose if a given MCMC run was likely long enough to have produced reliable estiamtes. The number of iterations and a thinning interval are suggested - if these are significantly larger than the parameters used to produce MCMCout, then the chains should be re-run with the suggested parameters.
#' @seealso \code{\link{run_MCMC}}
#' @author Philipp Neubauer
#' @references Neubauer,.P. and Jensen, O.P. (in prep)
#' @import coda
#' @export
diags <- function(MCMCout=NULL,accuracy=0.01,proba=0.95,quant=0.025){
  
  # revert to mcmc.list compatible format
  y  <- vector("list", MCMCout$nChains)
  for (i in 1:MCMCout$nChains) y[[i]] <- MCMCout[[i]]
  MCMCout <- y
  
  # do RL diag
  
  cat('\n','\n')
  
  cat('#################################','\n')
  
  cat('Raftery-Lewis diagnostics','\n')
  
  cat('#################################','\n','\n')
  
  class(MCMCout) <- 'mcmc.list'
  
  rd <- raftery.diag(MCMCout, q=quant, r=accuracy, s=proba)
  
  print(rd)
  
  if(rd[[1]]$resmatrix[1]=='Error'){
    
    cat('\n','Based on these diagnostics you should repeat the pilot MCMC with at least',rd[[1]]$resmatrix[2],'iterations to calculate diagnostics',' \n',' \n')
    
  } else{
    rit <- max(unlist(lapply(rd,function(x){x$resmatrix[,2]})))
    
    thin <- max(unlist(lapply(rd,function(x){x$resmatrix[,4]})))
    
    cat('\n','Based on these diagnostics you should repeat the MCMC with ',rit,' iterations','\n','and a thinning interval of ',round(thin),' ,if these values are higher than the values','\n','used to produce these diagnostics','\n','\n')
  }
  
  # do GR diag
  
  if(length(rd)>1){
    
    cat('\n','\n')
    
    cat('#################################','\n')
    
    cat('Gelman-Rubin diagnostics','\n')
    
    cat('#################################','\n','\n')
    
    print(gelman.diag(MCMCout,transform=T))
    
    cat('\n','Both univariate upper C.I. and multivariate psrf','\n','should be close to 1 if the chains converged','\n','\n','\n')
    
  }
}

.delist <- function(datas,Covs=NULL){
  
  options(warn=-1)
  
  n_preys =datas$n.preys
  m_preys =datas$n.preys-1
  n_preds =datas$n.preds
  eveness =datas$even
  if (!is.null(Covs)) {
    n_covs = ncol(Covs)
    ind = c(0,rep(1,n_covs-1))
    Covs = Covs[,1:n_covs]
  }
  
  jagsdata <- list(
    
    prey_ix = as.numeric(factor(datas$prey.ix,levels = unique(datas$prey.ix))),
    n_prey_samps = length(datas$prey.ix),
    n_preys =datas$n.preys,
    m_preys =datas$n.preys-1,
    n_preds =datas$n.preds,
    n_covs = ncol(Covs),
    
    S = diag(eveness,m_preys),
    SS = diag(1,m_preys),
    zeros = rep(0,m_preys),
    alpha = rep(1,n_preys),
    Covs = if(!is.null(Covs)) Covs,
    ind = if(!is.null(Covs)) ind
  )
  
  if(!is.null(datas$datas.FA))    jagsdata <- c(jagsdata,list(
    n_fats =datas$datas.FA$n.fats,
    R =datas$datas.FA$R,
    fc_mean =datas$datas.FA$fc_mean,
    fc_tau =datas$datas.FA$fc_tau,
    mean_c = data.frame(datas$datas.FA$mean_c),
    tau_coeffs =datas$datas.FA$tau_c,
    Rnot =datas$datas.FA$Rnot,
    m_fats =datas$datas.FA$m.fats,
    ni =datas$datas.FA$ni,
    preds = data.frame(datas$datas.FA$preds),
    preys = data.frame(datas$datas.FA$preys),
    preym =datas$datas.FA$preym)
  )
  
  
  if(!is.null(datas$datas.SI))    jagsdata <- c(jagsdata,list(
    prey_ix_SI = as.numeric(factor(datas$prey.ix.SI,levels = unique(datas$prey.ix.SI))),
    n_prey_samps_SI = length(datas$prey.ix.SI),
    isos =datas$datas.SI$isos,
    R_SI =datas$datas.SI$R.SI,
    mean_cs = as.matrix(datas$datas.SI$mean_cs),
    sigma_cs =as.matrix(1/datas$datas.SI$tau_cs),
    Rnot_SI =datas$datas.SI$Rnot.SI,
    ni_SI =datas$datas.SI$ni.SI,
    preds_SI = as.matrix(datas$datas.SI$preds.SI),
    preys_SI = as.matrix(datas$datas.SI$preys.SI),
    preym_SI =datas$datas.SI$preym.SI)
  )
  
  return(jagsdata)
  
}
