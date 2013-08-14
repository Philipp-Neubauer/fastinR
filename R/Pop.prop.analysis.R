Poppropanalysis <- function(datas,nIter=1000,nBurnin=1000,nChains=1,nThin=10) UseMethod("Poppropanalysis", datas)

Poppropanalysis.FA <- function(datas,nIter=10000,nBurnin=1000,nChains=1,nThin=10)
{
  n.preys = datas$n.preys
  R = datas$datas.FA$R
  fc_mean = datas$datas.FA$fc_mean
  fc_tau = datas$datas.FA$fc_tau
  n.fats = datas$datas.FA$n.fats
  mean_c = matrix(unlist(datas$datas.FA$mean_c),n.preys,n.fats)
  tau_coeffs = datas$datas.FA$tau_c
  Rnot = datas$datas.FA$Rnot
  n.preds = datas$n.preds
  m.fats = datas$datas.FA$m.fats
  ni = datas$datas.FA$ni
  preds = matrix(datas$datas.FA$preds,n.preds,m.fats)
  preym = datas$datas.FA$preym
  
  JM <- jags.model(file=paste(system.file("exec",package = "FASTIN"),"/Pop.prop.analysis.FA.bugs",sep=''),n.chains=nChains)
  
  cat('\n','proceeding to burn-in phase','\n')
  update(JM,n.iter=nBurnin)
  cat('\n','sampling from parameters','\n')
  res<- coda.samples(model=JM,variable.names='prop',n.iter=nIter,thin=nThin)
  
  res <- as.data.frame((res)[[1]])
  output <- list(MCMC=res)
  class(output) <- 'pop_props'
  return(output)
  
}

Poppropanalysis.SI <- function(datas,nIter=10000,nBurnin=1000,nChains=1,nThin=10)
{
  n.preys = datas$n.preys
  n.preds = datas$n.preds
  
  isos = datas$datas.SI$isos
  R_SI = datas$datas.SI$R.SI
  mean_cs = matrix(unlist(datas$datas.SI$mean_cs),n.preys,isos)
  tau_cs = datas$datas.SI$tau_cs
  Rnot_SI = datas$datas.SI$Rnot.SI
  ni.SI = datas$datas.SI$ni.SI
  preds.SI = matrix(unlist(datas$datas.SI$preds.SI),n.preds,isos)
  preym.SI = datas$datas.SI$preym.SI
  
  JM <- jags.model(file=paste(system.file("exec",package = "FASTIN"),"/Pop.prop.analysis.SI.bugs",sep=''),n.chains=nChains)
  
  cat('\n','proceeding to burn-in phase','\n')
  update(JM,n.iter=nBurnin)
  cat('\n','sampling from parameters','\n')
  res<- coda.samples(model=JM,variable.names='prop',n.iter=nIter,thin=nThin)
  res <- as.data.frame((res)[[1]])
  output <- list(MCMC=res)
  class(output) <- 'pop_props'
  return(output)
}

Poppropanalysis.combined <- function(datas,nIter=10000,nBurnin=1000,nChains=1,nThin=10)
{
  
  n.preys = datas$n.preys
  n.preds = datas$n.preds
  
  R_SI = datas$datas.SI$R.SI
  isos = datas$datas.SI$isos
  mean_cs = matrix(unlist(datas$datas.SI$mean_cs),n.preys,isos)
  tau_cs = datas$datas.SI$tau_cs
  Rnot_SI = datas$datas.SI$Rnot.SI
  ni.SI = datas$datas.SI$ni.SI
  preds.SI = datas$datas.SI$preds.SI
  preym.SI = datas$datas.SI$preym.SI
  
  R = datas$datas.FA$R
  fc_mean = datas$datas.FA$fc_mean
  fc_tau = datas$datas.FA$fc_tau
  n.fats = datas$datas.FA$n.fats
  mean_c = matrix(unlist(datas$datas.FA$mean_c),n.preys,n.fats)
  tau_coeffs = data.matrix(datas$datas.FA$tau_c)
  Rnot = datas$datas.FA$Rnot
  m.fats = datas$datas.FA$m.fats
  ni = datas$datas.FA$ni
  preds = datas$datas.FA$preds[,]
  preym = datas$datas.FA$preym[,]
  
  JM <- jags.model(file=paste(system.file("exec",package = "FASTIN"),"/Pop.prop.analysis.combined.bugs",sep=''),n.chains=nChains)
  
  cat('\n','proceeding to burn-in phase','\n')
  update(JM,n.iter=nBurnin)
  cat('\n','sampling from parameters','\n')
  res<- coda.samples(model=JM,variable.names='prop',n.iter=nIter,thin=nThin)
  
  res <- as.data.frame((res)[[1]])
  output <- list(MCMC=res)
  class(output) <- 'pop_props'
  return(output)
  
}