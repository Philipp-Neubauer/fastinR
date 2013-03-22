AnalysiswithCov <- function(data,Covs,nIter=1000,nBurnin=1000) UseMethod("PopandIndprops", data)

AnalysiswithCov.FA <- function(data,Covs,nIter=1000,nBurnin=1000)
{
  n.covs = ncol(Covs)
  R = data$datas.FA$R
  fc_mean = data$datas.FA$fc_mean
  fc_tau = data$datas.FA$fc_tau
  mean_c = data$datas.FA$mean_c
  tau_coeffs = data$datas.FA$tau_c
  Rnot = data$datas.FA$Rnot
  n.preys = data$n.preys
  n.preds = data$n.preds
  n.fats = data$datas.FA$n.fats
  m.fats = data$datas.FA$m.fats
  ni = data$datas.FA$ni
  preds = data$datas.FA$preds
  preym = data$datas.FA$preym
  eveness = data$eveness
  
  S = diag(eveness,n.preys)
  SS = diag(1,n.preys)
  
  initials=list(list(
    beta=matrix(0,n.preys,n.covs),
    fc = fc_mean,
    fracs=matrix(mean_c,n.preys,n.fats),     
    pprec = diag(1,n.preys),
    prey.means=preym,
    predprec = diag(1,m.fats),
    prey.precs = array(1,c(n.preys,m.fats,m.fats))
  ))
    
  datas=list('Covs','n.covs','S','SS','R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','ni','preds','preym')
  
  vars = c('prop','pop.prop','beta')
  
  # compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...
  
  res <- BRugsFit('Pop.prop.analysis.FA.bugs', datas, inits=initials, numChains = 1, vars,
                                 nBurnin = nBurnin, nIter = nIter, nThin = round(nIter/1000), coda = T,
                                 DIC = F, working.directory = getwd(), digits = 4,
                                 BRugsVerbose = T)
  
  res <- as.data.frame((res)[[1]])
  output <- list(MCMC=res,Covnames=names(Covs),Preynames=rownames(preys))
  class(output) <- 'cov_props'
  return(output)
  
}

AnalysiswithCov.SI <- function(data,Covs,nIter=1000,nBurnin=1000)
{
  n.preys = data$n.preys
  n.preds = data$n.preds
  
  n.covs = ncol(Covs)
  
  R.SI = data$datas.SI$R.SI
  mean_cs = data$datas.SI$mean_cs
  tau_cs = data$datas.SI$tau_cs
  Rnot.SI = data$datas.SI$Rnot.SI
  isos = data$datas.SI$isos
  ni.SI = data$datas.SI$ni.SI
  preds.SI = data$datas.SI$preds.SI
  preym.SI = data$datas.SI$preym.SI
  
  S = diag(eveness,n.preys)
  SS = diag(1,n.preys)
  
  initials.SI=list(list(
    beta=matrix(0,n.preys,n.covs),
    pmean=rep(0,n.preys),
    pprec = diag(1,n.preys), 
    cs=matrix(mean_cs,n.preys,isos,byrow=T),
    prey.means_SI=preym.SI,
    predprec_SI = diag(0.01,isos),
    prey.precs_SI = array(1,c(n.preys,isos,isos))
    
  ))
  
  
  
  datas.SI=list('Covs','n.covs','S','SS','R_SI','mean_cs','tau_cs','Rnot_SI','n.preys','n.preds','isos','ni.SI','preds.SI','preym.SI')
  
  vars = c('prop','pop.prop','beta')
  
  # compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...
  
  res <- BRugsFit('Pop.prop.analysis.SI.bugs', datas.SI, inits=initials.SI, numChains = 1, vars,
                  nBurnin = nBurnin, nIter = Iter, nThin = round(nIter/1000), coda = T,
                  DIC = F, working.directory = getwd(), digits = 4, 
                  BRugsVerbose = T)
  
  output <- list(MCMC=res,Covnames=names(Covs),Preynames=rownames(preys.SI))
  class(output) <- 'cov_props'
  return(output)
}

AnalysiswithCov.combined <- function(data,Covs,nIter=1000,nBurnin=1000)
{
  
  n.preys = data$n.preys
  n.preds = data$n.preds
  n.covs = ncol(Covs)
  R.SI = data$datas.SI$R.SI
  mean_cs = data$datas.SI$mean_cs
  tau_cs = data$datas.SI$tau_cs
  Rnot.SI = data$datas.SI$Rnot.SI
  isos = data$datas.SI$isos
  ni.SI = data$datas.SI$ni.SI
  preds.SI = data$datas.SI$preds.SI
  preym.SI = data$datas.SI$preym.SI
  
  R = data$datas.FA$R
  fc_mean = data$datas.FA$fc_mean
  fc_tau = data$datas.FA$fc_tau
  mean_c = data$datas.FA$mean_c
  tau_coeffs = data$datas.FA$tau_c
  Rnot = data$datas.FA$Rnot
  n.fats = data$datas.FA$n.fats
  m.fats = data$datas.FA$m.fats
  ni = data$datas.FA$ni
  preds = data$datas.FA$preds
  preym = data$datas.FA$preym
  
  S = diag(eveness,n.preys)
  SS = diag(1,n.preys)
 
  initials.comb=list(list(
    beta=matrix(0,n.preys,n.covs),
    pprec = diag(1,n.preys), 
    fc = fc_mean,
    fracs=matrix(mean_c,n.preys,n.fats),             
    prey.means=preym,
    predprec = diag(1,m.fats),
    prey.precs = array(0.1,c(n.preys,m.fats,m.fats)),
    cs=matrix(mean_cs,n.preys,isos,byrow=T),
    prey.means_SI=preym.SI,
    predprec_SI = diag(0.01,isos),
    prey.precs_SI = array(1,c(n.preys,isos,isos))
    
  ))
  
  datas.comb=list('Covs','n.covs','S','SS','R','R_SI','fc_mean','fc_tau','mean_c','tau_coeffs','mean_cs','tau_cs','Rnot','Rnot_SI','n.preys','n.preds','isos','n.fats','m.fats','ni.SI','ni','preds','preds.SI','preym.SI','preym')
  
  vars = c('prop','pop.prop','beta')
  
  # compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...
  
  res <- BRugsFit('Pop.prop.analysis.combined.bugs', datas.comb, inits=initials.comb, numChains = 1, vars,
                                   nBurnin = nBurnin, nIter = Iter, nThin = round(nIter/1000), coda = T,
                                   DIC = F, working.directory = getwd(), digits = 4, 
                                   BRugsVerbose = T)
  res <- as.data.frame((res)[[1]])
  output <- list(MCMC=res,Covnames=names(Covs),Preynames=rownames(preys))
  class(output) <- 'cov_props'
  return(output)
  
}