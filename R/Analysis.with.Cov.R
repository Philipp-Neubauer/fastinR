AnalysiswithCov <- function(datas,Covs,nIter=1000,nBurnin=1000) UseMethod("AnalysiswithCov",datas)

AnalysiswithCov.FA <- function(datas,Covs,nIter=1000,nBurnin=1000)
{
  n.covs = ncol(Covs)
  R =datas$datas.FA$R
  fc_mean =datas$datas.FA$fc_mean
  fc_tau =datas$datas.FA$fc_tau
  mean_c = matrix(unlist(data$datas.SI$mean_cs),n.preys,isos)
  tau_coeffs =datas$datas.FA$tau_c
  Rnot =datas$datas.FA$Rnot
  n.preys =datas$n.preys
  n.preds =datas$n.preds
  n.fats =datas$datas.FA$n.fats
  m.fats =datas$datas.FA$m.fats
  ni =datas$datas.FA$ni
  preds =datas$datas.FA$preds
  preym =datas$datas.FA$preym
  eveness =datas$eveness
  
  S = diag(eveness,n.preys)
  SS = diag(1,n.preys)
  
  initials=list(list(
    beta=matrix(0,n.preys,n.covs),
    fc = fc_mean,
    fracs=mean_c,     
    pprec = diag(1,n.preys),
    prey.means=preym,
    predprec = diag(1,m.fats),
    prey.precs = array(1,c(n.preys,m.fats,m.fats))
  ))
    
 datass=list('Covs','n.covs','S','SS','R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','ni','preds','preym')
  
  vars = c('prop','pop.prop','beta')
  
  # compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...
  
  res <- BRugsFit(system.file("exec","Analysis.with.Cov.FA.bugs",package = 'FASTIN'),datass, inits=initials, numChains = 1, vars,
                                 nBurnin = nBurnin, nIter = nIter, nThin = round(nIter/1000), coda = T,
                                 DIC = F, working.directory = getwd(), digits = 4,
                                 BRugsVerbose = T)
  
  res <- as.data.frame((res)[[1]])
  output <- list(MCMC=res,Covnames=names(Covs))
  class(output) <- 'cov_props'
  return(output)
  
}

AnalysiswithCov.SI <- function(datas,Covs,nIter=1000,nBurnin=1000)
{
  n.preys =datas$n.preys
  n.preds =datas$n.preds
  
  n.covs = ncol(Covs)
  
  R_SI =datas$datas.SI$R.SI
  mean_cs = matrix(unlist(data$datas.SI$mean_cs),n.preys,isos)
  tau_cs =datas$datas.SI$tau_cs
  Rnot_SI =datas$datas.SI$Rnot.SI
  isos =datas$datas.SI$isos
  ni.SI =datas$datas.SI$ni.SI
  preds.SI =datas$datas.SI$preds.SI
  preym.SI =datas$datas.SI$preym.SI
  eveness =datas$eveness
  
  S = diag(eveness,n.preys)
  SS = diag(1,n.preys)
  
  initials.SI=list(list(
    beta=matrix(0,n.preys,n.covs),
    pmean=rep(0,n.preys),
    pprec = diag(1,n.preys), 
    cs=mean_cs,
    prey.means_SI=preym.SI,
    predprec_SI = diag(0.01,isos),
    prey.precs_SI = array(1,c(n.preys,isos,isos))
    
  ))
  
  
  
 datass.SI=list('Covs','n.covs','S','SS','R_SI','mean_cs','tau_cs','Rnot_SI','n.preys','n.preds','isos','ni.SI','preds.SI','preym.SI')
  
  vars = c('prop','pop.prop','beta')
  
  # compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...
  
  res <- BRugsFit(system.file("exec","Analysis.with.Cov.SI.bugs",package = 'FASTIN'),datass.SI, inits=initials.SI, numChains = 1, vars,
                  nBurnin = nBurnin, nIter = Iter, nThin = round(nIter/1000), coda = T,
                  DIC = F, working.directory = getwd(), digits = 4, 
                  BRugsVerbose = T)
  
  output <- list(MCMC=res,Covnames=names(Covs))
  class(output) <- 'cov_props'
  return(output)
}

AnalysiswithCov.combined <- function(datas,Covs,nIter=1000,nBurnin=1000)
{
  
  n.preys =datas$n.preys
  n.preds =datas$n.preds
  n.covs = ncol(Covs)
  R_SI =datas$datas.SI$R.SI
  mean_cs = matrix(unlist(datas$datas.SI$mean_cs),n.preys,isos)
  tau_cs =datas$datas.SI$tau_cs
  Rnot_SI =datas$datas.SI$Rnot.SI
  isos =datas$datas.SI$isos
  ni.SI =datas$datas.SI$ni.SI
  preds.SI =datas$datas.SI$preds.SI
  preym.SI =datas$datas.SI$preym.SI
  
  R =datas$datas.FA$R
  fc_mean =datas$datas.FA$fc_mean
  fc_tau =datas$datas.FA$fc_tau
  mean_c =matrix(unlist(datas$datas.FA$mean_c),n.preys,n.fats)
  tau_coeffs =datas$datas.FA$tau_c
  Rnot =datas$datas.FA$Rnot
  n.fats =datas$datas.FA$n.fats
  m.fats =datas$datas.FA$m.fats
  ni =datas$datas.FA$ni
  preds =datas$datas.FA$preds
  preym =datas$datas.FA$preym
  eveness =datas$eveness
  
  S = diag(eveness,n.preys)
  SS = diag(1,n.preys)
 
  initials.comb=list(list(
    beta=matrix(0,n.preys,n.covs),
    pprec = diag(1,n.preys), 
    fc = fc_mean,
    fracs=mean_c,             
    prey.means=preym,
    predprec = diag(1,m.fats),
    prey.precs = array(0.1,c(n.preys,m.fats,m.fats)),
    cs=mean_cs,
    prey.means_SI=preym.SI,
    predprec_SI = diag(0.01,isos),
    prey.precs_SI = array(1,c(n.preys,isos,isos))
    
  ))
  
 datass.comb=list('Covs','n.covs','S','SS','R','R_SI','fc_mean','fc_tau','mean_c','tau_coeffs','mean_cs','tau_cs','Rnot','Rnot_SI','n.preys','n.preds','isos','n.fats','m.fats','ni.SI','ni','preds','preds.SI','preym.SI','preym')
  
  vars = c('prop','pop.prop','beta')
  
  # compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...
  system.file("exec","Analysis.with.Cov.combined.bugs",package = 'FASTIN')
  res <- BRugsFit(system.file("exec","Analysis.with.Cov.combined.bugs",package = 'FASTIN'),datass.comb, inits=initials.comb, numChains = 1, vars,
                                   nBurnin = nBurnin, nIter = Iter, nThin = round(nIter/1000), coda = T,
                                   DIC = F, working.directory = getwd(), digits = 4, 
                                   BRugsVerbose = T)
  res <- as.data.frame((res)[[1]])
  output <- list(MCMC=res,Covnames=names(Covs))
  class(output) <- 'cov_props'
  return(output)
  
}