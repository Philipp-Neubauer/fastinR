PopandIndprops <- function(datas,nIter=10000,nBurnin=1000,nChains=1,nThin=10) UseMethod("PopandIndprops", datas)

PopandIndprops.FA <- function(datas,nIter=10000,nBurnin=1000,nChains=1,nThin=10)
{
  n.preys = datas$n.preys
  n.fats = datas$datas.FA$n.fats
  R = datas$datas.FA$R
  fc_mean = datas$datas.FA$fc_mean
  fc_tau = datas$datas.FA$fc_tau
  mean_c = matrix(unlist(datas$datas.FA$mean_c),n.preys,n.fats)
  tau_coeffs = datas$datas.FA$tau_c
  Rnot = datas$datas.FA$Rnot
  n.preds = datas$n.preds
  m.fats = datas$datas.FA$m.fats
  ni = datas$datas.FA$ni
  preds = datas$datas.FA$preds
  preym = datas$datas.FA$preym
  
  S = diag(eveness,n.preys)
  SS = diag(1,n.preys)
  zeros = rep(0,n.preys)
  
  JM <- jags.model(file=paste(system.file("exec",package = "FASTIN"),"/Pop.and.Ind.props.FA.bugs",sep=''),n.chains=nChains)
  
  cat('\n','proceeding to burn-in phase','\n')
  update(JM,n.iter=nBurnin)
  cat('\n','sampling from parameters','\n')
  res<- coda.samples(model=JM,variable.names=c('prop','pop.prop'),n.iter=nIter,thin=nThin)
  
  res <- as.data.frame((res)[[1]])
  output <- list(MCMC=res)
  class(output) <- 'ind_props'
  return(output)
  
}

PopandIndprops.SI <- function(datas,nIter=10000,nBurnin=1000,nChains=1,nThin=10)
{
  n.preys = datas$n.preys
  n.preds = datas$n.preds
  
  isos = datas$datas.SI$isos
  R_SI = datas$datas.SI$R.SI
  mean_cs = matrix(unlist(datas$datas.SI$mean_cs),n.preys,isos)
  tau_cs = datas$datas.SI$tau_cs
  Rnot_SI = datas$datas.SI$Rnot.SI
  ni.SI = datas$datas.SI$ni.SI
  preds.SI = datas$datas.SI$preds.SI
  preym.SI = datas$datas.SI$preym.SI
  eveness = datas$even
    
  S = diag(1,n.preys)
  SS = diag(1,n.preys)
  zeros = rep(0,n.preys)
  
  JM <- jags.model(file=paste(system.file("exec",package = "FASTIN"),"/Pop.and.Ind.props.SI.bugs",sep=''),n.chains=nChains)
  
  cat('\n','proceeding to burn-in phase','\n')
  update(JM,n.iter=nBurnin)
  cat('\n','sampling from parameters','\n')
  res<- coda.samples(model=JM,variable.names=c('prop','pop.prop'),n.iter=nIter,thin=nThin)

summary(res)
plot(res[,1])
  
  xres <- as.data.frame((res)[[1]])
  output <- list(MCMC=res)
  class(output) <- 'ind_props'
  return(output)
}

PopandIndprops.combined <- function(datas,nIter=10000,nBurnin=1000,nChains=1,nThin=10)
{
  
  n.preys = datas$n.preys
  n.preds = datas$n.preds
  
  isos = datas$datas.SI$isos
  R_SI = datas$datas.SI$R.SI
  mean_cs = matrix(unlist(datas$datas.SI$mean_cs),n.preys,isos)
  tau_cs = datas$datas.SI$tau_cs
  Rnot_SI = datas$datas.SI$Rnot.SI
  ni.SI = datas$datas.SI$ni.SI
  preds.SI = datas$datas.SI$preds.SI
  preym.SI = datas$datas.SI$preym.SI
  
  n.fats = datas$datas.FA$n.fats
  R = datas$datas.FA$R
  fc_mean = datas$datas.FA$fc_mean
  fc_tau = datas$datas.FA$fc_tau
  mean_c =  matrix(unlist(datas$datas.FA$mean_c),n.preys,n.fats)
  tau_coeffs = datas$datas.FA$tau_c
  Rnot = datas$datas.FA$Rnot
  m.fats = datas$datas.FA$m.fats
  ni = datas$datas.FA$ni
  preds = datas$datas.FA$preds
  preym = datas$datas.FA$preym
  eveness = datas$even
  
  S = diag(eveness,n.preys)
  SS = diag(1,n.preys)
  zeros=rep(0,n.preys)
  
  JM <- jags.model(file=paste(system.file("exec",package = "FASTIN"),"/Pop.and.Ind.props.combined.bugs",sep=''),n.chains=nChains)
  
  cat('\n','proceeding to burn-in phase','\n')
  update(JM,n.iter=nBurnin)
  cat('\n','sampling from parameters','\n')
  res<- coda.samples(model=JM,variable.names=c('prop','pop.prop'),n.iter=nIter,thin=nThin)
  
  res <- as.data.frame((res)[[1]])
  output <- list(MCMC=res)
  class(output) <- 'ind_props'
  return(output)
  
}
