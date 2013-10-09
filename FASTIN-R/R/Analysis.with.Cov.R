.AnalysiswithCov <- function(datas,Covs,nIter=10000,nBurnin=1000,nChains=1,nThin=10) UseMethod("AnalysiswithCov",datas)

.AnalysiswithCov.FA <- function(datas,Covs,nIter=10000,nBurnin=1000,nChains=1,nThin=10)
{

    n.preys =datas$n.preys
  m.preys =datas$n.preys-1
  n.preds =datas$n.preds
   eveness =datas$even
n.covs = ncol(Covs)

     jagsdata <- list(
  n.preys =datas$n.preys,
  m.preys =datas$n.preys-1,
  n.fats =datas$datas.FA$n.fats,
  n.covs = ncol(Covs),
  R =datas$datas.FA$R,
  fc_mean =datas$datas.FA$fc_mean,
  fc_tau =datas$datas.FA$fc_tau,
  mean_c = matrix(unlist(datas$datas.FA$mean_c),n.preys,datas$datas.FA$n.fats),
  tau_coeffs =datas$datas.FA$tau_c,
  Rnot =datas$datas.FA$Rnot,
  n.preds =datas$n.preds,
  m.fats =datas$datas.FA$m.fats,
  ni =datas$datas.FA$ni,
  preds =datas$datas.FA$preds,
  preym =datas$datas.FA$preym,   
  S = diag(eveness,m.preys),
  SS = diag(1,m.preys),
  zeros = rep(0,m.preys),
  Covs = Covs[,1:n.covs],
  ind = c(0,rep(1,n.covs-1)))
    
  JM <- jags.model(file=paste(system.file("exec",package = "FASTIN"),"/Analysis.with.Cov.FA.bugs",sep=''),n.chains=nChains,data=jagsdata)
  
  cat('\n','proceeding to burn-in phase','\n')
  update(JM,n.iter=nBurnin)
  cat('\n','sampling from parameters','\n')
  res<- coda.samples(model=JM,variable.names=c('prop','pop.prop','beta'),n.iter=nIter,thin=nThin)
    
  res <- as.data.frame((res)[[1]])
  output <- list(MCMC=res,Covnames=names(Covs))
  class(output) <- 'cov_props'
    return(output)
  
}

.AnalysiswithCov.SI <- function(datas,Covs,nIter=10000,nBurnin=1000,nChains=1,nThin=10)
{
  n.preys =datas$n.preys
  m.preys =datas$n.preys-1
  n.preds =datas$n.preds
   eveness =datas$even
n.covs = ncol(Covs)
  
  jagsdata <- list(
  n.preys =datas$n.preys,
  m.preys =datas$n.preys-1,
  n.preds =datas$n.preds,
  n.covs = ncol(Covs),
  isos =datas$datas.SI$isos,
  R_SI =datas$datas.SI$R.SI,
  mean_cs = matrix(unlist(datas$datas.SI$mean_cs),n.preys,datas$datas.SI$isos),
  tau_cs =datas$datas.SI$tau_cs,
  Rnot_SI =datas$datas.SI$Rnot.SI,
  ni.SI =datas$datas.SI$ni.SI,
  preds.SI =datas$datas.SI$preds.SI,
  preym.SI =datas$datas.SI$preym.SI,
  S = diag(eveness,m.preys),
  SS = diag(1,m.preys),
  zeros = rep(0,m.preys),
  Covs = Covs[,1:n.covs],
  ind = c(0,rep(1,n.covs-1))
      )
  
  JM <- jags.model(file=paste(system.file("exec",package = "FASTIN"),"/Analysis.with.Cov.SI.bugs",sep=''),inits=list(beta=matrix(0,m.preys,n.covs)),n.chains=nChains,data=jagsdata)
  
  cat('\n','proceeding to burn-in phase','\n')
  update(JM,n.iter=nBurnin)
  cat('\n','sampling from parameters','\n')
  res<- coda.samples(model=JM,variable.names=c('prop','pop.prop','beta'),n.iter=nIter,thin=nThin)

  res <- as.data.frame((res)[[1]])
  
  output <- list(MCMC=res,Covnames=names(Covs))
  class(output) <- 'cov_props'
  return(output)
}

.AnalysiswithCov.combined <- function(datas,Covs,nIter=10000,nBurnin=1000,nChains=1,nThin=10)
{

    n.preys =datas$n.preys
  m.preys =datas$n.preys-1
  n.preds =datas$n.preds
   eveness =datas$even
n.covs = ncol(Covs)
  
  jagsdata <- list(
  n.fats =datas$datas.FA$n.fats,
  R =datas$datas.FA$R,
  fc_mean =datas$datas.FA$fc_mean,
  fc_tau =datas$datas.FA$fc_tau,
  mean_c = matrix(unlist(datas$datas.FA$mean_c),n.preys,datas$datas.FA$n.fats),
  tau_coeffs =datas$datas.FA$tau_c,
  Rnot =datas$datas.FA$Rnot,
  m.fats =datas$datas.FA$m.fats,
  ni =datas$datas.FA$ni,
  preds =datas$datas.FA$preds,
  preym =datas$datas.FA$preym,   
  n.preys =datas$n.preys,
  m.preys =datas$n.preys-1,
  n.preds =datas$n.preds,
  n.covs = ncol(Covs),
  isos =datas$datas.SI$isos,
  R_SI =datas$datas.SI$R.SI,
  mean_cs = matrix(unlist(datas$datas.SI$mean_cs),n.preys,datas$datas.SI$isos),
  tau_cs =datas$datas.SI$tau_cs,
  Rnot_SI =datas$datas.SI$Rnot.SI,
  ni.SI =datas$datas.SI$ni.SI,
  preds.SI =datas$datas.SI$preds.SI,
  preym.SI =datas$datas.SI$preym.SI,
  S = diag(eveness,m.preys),
  SS = diag(1,m.preys),
  zeros = rep(0,m.preys),
  Covs = Covs[,1:n.covs],
  ind = c(0,rep(1,n.covs-1))
      )
 
  JM <- jags.model(file=paste(system.file("exec",package = "FASTIN"),"/Analysis.with.Cov.combined.bugs",sep=''),n.chains=nChains,data=jagsdata)
  
  cat('\n','proceeding to burn-in phase','\n')
  update(JM,n.iter=nBurnin)
  cat('\n','sampling from parameters','\n')
  res<- coda.samples(model=JM,variable.names=c('prop','pop.prop','beta'),n.iter=nIter,thin=nThin)
  
  res <- as.data.frame((res)[[1]])
  output <- list(MCMC=res,Covnames=names(Covs))
  class(output) <- 'cov_props'
  return(output)
  
}
