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

run_MCMC <- function(datas=NULL,Covs=NULL,nIter=10000,nBurnin=1000,nChains=1,nThin=10,Data.Type='Fatty.Acid.Profiles',Analysis.Type='Population.proportions',even=5,plott=T){
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
    if(GUI==T & Analysis.Type == 'Analysis.with.Covariates') {
      Covs = guiGetSafe('Covs')
      if(any(is.na(Covs)) & Analysis.Type == 'Analysis.with.Covariates')
      {
        stop('analysis with covariates selected, but no covariates entered.')
      }
    }
    
    datas$even <- even
    jagsdata <- delist(datas)
    
    save(jagsdata,file='jagsdata.Rdata')
    
    sysfile <- switch(Analysis.Type,
                      Population.proportions = .Poppropanalysis(datas),
                      Individual.proportions = .PopandIndprops(datas),
                      Analysis.with.Covariates = .AnalysiswithCov(datas)
    )
    
    res <- spawn(sysfile,nChains,nBurnin,nIter,nThin)
    
    if (plott){
      plotta <- menu(title='plot MCMC chains?',choices = c('yes','no'),graphics=T)
      if (plotta==1) {
        try(plot(res,ask=T))
        }
      }
    
    class(res) <- switch(Analysis.Type,
                      Population.proportions = 'pop_props',
                      Individual.proportions = 'ind_props',
                      Analysis.with.Covariates = 'cov_props'
    )
    
    ifelse(GUI==F,return(res),guiSet('MCMCout',res))
  }
}

# this is called in the slave processes to run jags
jagger <- function(sysfile,nBurnin,nIter,nThin,i){
  
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

#spawn slave processes to run separate chains - ultimately this should integrate with runjags0
spawn <- function(sysfile,nChains,nBurnin,nIter,nThin){
  
  orgtime <- Sys.time()
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
      system(cmd,wait=F)
    }
    
  }
  
  grep('MCout',dir())
  while(!any(grep('MCout',dir()))){Sys.sleep(3);cat('+')}
  while(length(grep('MCout',dir()))<nChains | any(file.info(dir()[grep('MCout',dir())])$mtime <  orgtime)) {Sys.sleep(3);cat('+')}
  
  
  res <- eval(parse(text=eval(load('MCout1.Rdata'))))
  
  if (nChains>1) {
    for(i in 2:nChains)
      res[[i]] <- eval(parse(text=eval(load(paste('MCout',i,'.Rdata',sep='')))))[[1]]
  }
  
  return(res)
  
}

diags <- function(MCMCout=NULL,accuracy=0.01,proba=0.95,quant=0.025){
  
  # do RL diag
 
  cat('\n','\n')
  
  cat('#################################','\n')
 
  cat('Raftery-Lewis diagnostics','\n')
 
  cat('#################################','\n','\n')

  class(MCMCout) <- 'mcmc.list'
  
  rd <- raftery.diag(MCMCout, q=quant, r=accuracy, s=proba)
  
  print(rd)
  
  if(rd[[1]]$resmatrix[1]=='Error'){
    
    cat('\n','Based on these diagnostics you should repeat the pilot MCMC with at least ',rd[[1]]$resmatrix[2],' \n','iterations to calculate diagnostics',' \n',' \n')
    
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