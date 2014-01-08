delist <- function(datas,Covs=NULL){
  
  n.preys =datas$n.preys
  m.preys =datas$n.preys-1
  n.preds =datas$n.preds
  eveness =datas$even
  if (!is.null(Covs)) {
    n.covs = ncol(Covs)
    ind = c(0,rep(1,n.covs-1))
    Covs = Covs[,1:n.covs]
  }
  
  jagsdata <- list(
    n.fats =datas$datas.FA$n.fats,
    R =datas$datas.FA$R,
    fc_mean =datas$datas.FA$fc_mean,
    fc_tau =datas$datas.FA$fc_tau,
    mean_c = data.frame(datas$datas.FA$mean_c),
    tau_coeffs =datas$datas.FA$tau_c,
    Rnot =datas$datas.FA$Rnot,
    m.fats =datas$datas.FA$m.fats,
    ni =datas$datas.FA$ni,
    preds = data.frame(datas$datas.FA$preds),
    preym =datas$datas.FA$preym,   
    n.preys =datas$n.preys,
    m.preys =datas$n.preys-1,
    n.preds =datas$n.preds,
    n.covs = ncol(Covs),
    isos =datas$datas.SI$isos,
    R_SI =datas$datas.SI$R.SI,
    mean_cs = data.frame(datas$datas.SI$mean_cs),
    tau_cs =datas$datas.SI$tau_cs,
    Rnot_SI =datas$datas.SI$Rnot.SI,
    ni.SI =datas$datas.SI$ni.SI,
    preds.SI = data.frame(datas$datas.SI$preds.SI),
    preym.SI =datas$datas.SI$preym.SI,
    S = diag(eveness,m.preys),
    SS = diag(1,m.preys),
    zeros = rep(0,m.preys),
    Covs = ifelse(!is.null(Covs),Covs,NA),
    ind = ifelse(!is.null(Covs),ind,NA)
  )
  
  return(jagsdata)
  
}

run_MCMC <- function(datas=NULL,nIter=10000,nBurnin=1000,nChains=1,nThin=10,Data.Type='Fatty.Acid.Profiles',Analysis.Type='Population.proportions',even=0.1,plott=T){
  # have three types here: FA, SI and combined, then methods dispatch based on type of arg
  
  if(missing(datas)) {datas = guiGetSafe('datas');GUI = T} else {GUI = F}
    
  if(length(datas)<=1){warning('No data processed yet')} else {   
    
    # Assign data type
    if (Data.Type == 'Combined.Analysis')
    {
      class(datas) <- 'combined'
    } else if (Data.Type=='Fatty.Acid.Profiles')
    {
      class(datas) <- 'FA'
    } else if (Data.Type=='Stable.Isotopes')
    {
      class(datas) <- 'SI'
    }
    
    # check for covariates
    if(GUI==T & Analysis.Type == 'Analysis.with.Covariates') Covs = guiGetSafe('Covs')
    
    if(any(is.na(Covs)) & Analysis.Type == 'Analysis.with.Covariates'){stop('analysis with covariates selected, but no covariates entered.')}
    
    jagsdata <- delist(datas)
    
    save(jagsdata,file='jagsdata.Rdata')
    
    sysfile <- switch(Analysis.Type,
                      Population.proportions = .Poppropanalysis(datas),
                      Individual.proportions = .PopandIndprops(datas),
                      Analysis.with.Covariates = .AnalysiswithCov(datas)
    )
    
    res <- spawn(sysfile,nBurnin,nIter,nThin)
    
    if (plott){
      plotta <- menu(title='plot MCMC chains?',choices = c('yes','no'),graphics=T)
      if (plotta==1) plot(res,ask=T)}
    
    class(res) <- switch(Analysis.Type,
                      Population.proportions = 'pop_props',
                      Individual.proportions = 'ind_props',
                      Analysis.with.Covariates = 'cov_props'
    )
    
    ifelse(GUI==F,return(outputs),guiSet('MCMCout',outputs))
  }
}

# this is called in the slave processes to run jags
jagger <- function(sysfile,nBurnin,nIter,nThin,i){
  
  load('jagsdata.Rdata')
  JM <- jags.model(file=sysfile,data=jagsdata,inits=list(.RNG.name="base::Mersenne-Twister",.RNG.seed=i),n.chains=1)
  
  cat('\n','proceeding to burn-in phase','\n')
  update(JM,n.iter=nBurnin)
  cat('\n','sampling from parameter distributions','\n')
  res<- coda.samples(model=JM,variable.names='prop',n.iter=nIter,thin=nThin)
  
  return(res)
  
}

#spawn slave processes to run separate chains - ultimately this should integrate with runjags0
spawn <- function(sysfile,nBurnin,nIter,nThin){
  
  options(useFancyQuotes=F)
  # spawn seperate R slave processes for each chain and run jagger in them
  for (i in 1:nChains){
    outfile <- paste('MCout',i,'.Rdata',sep='')
    
    cmd <- paste("Rscript -e 'library(FASTIN,quietly=T,verbose=F);options(warn=-1);res",i,"<- jagger(",dQuote(sysfile),",",eval(nBurnin),",",eval(nIter),",",eval(nThin),",",i,");save.image(file=",dQuote(outfile),");print(",dQuote("all done"),")'",sep='')
    
    if(i<nChains)
    {
      system(cmd,wait=F)
    }
    else {
      system(cmd,wait=T)
    }
    
  }
  
  res <- eval(parse(text=eval(load('MCout1.Rdata'))))
  
  if (nChains>1) {
    for(i in 2:nChains)
      res[[i]] <- eval(parse(text=eval(load(paste('MCout',i,'.Rdata',sep='')))))[[1]]
  }
  
  return(res)
  
}