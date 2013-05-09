PopandIndprops <- function(data,nIter=1000,nBurnin=1000) UseMethod("PopandIndprops", data)

PopandIndprops.FA <- function(data,nIter=1000,nBurnin=1000)
{
  
  R = data$datas.FA$R
  fc_mean = data$datas.FA$fc_mean
  fc_tau = data$datas.FA$fc_tau
  mean_c = matrix(unlist(data$datas.FA$mean_c),n.preys,n.fats)
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
  zeros = rep(0,n.preys)
  
  initials=list(list(
    fc = fc_mean,
    fracs=mean_c,     
    pmean=rep(0,n.preys),
    pnorm=matrix(0,n.preds,n.preys),
    pprec = diag(1,n.preys),
    prey.means=preym,
    predprec = diag(1,m.fats),
    prey.precs = array(1,c(n.preys,m.fats,m.fats))
  ))
    
  datas=list('zeros','S','SS','R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','ni','preds','preym')
  
  vars = c('prop','pop.prop')
  
  # compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...
  
  res <- BRugsFit(system.file("exec","Pop.and.Ind.props.FA.bugs",package = 'FASTIN'), datas, inits=initials, numChains = 1, vars,
                                 nBurnin = nBurnin, nIter = nIter, nThin = round(nIter/1000), coda = T,
                                 DIC = F, working.directory = getwd(), digits = 4,
                                 BRugsVerbose = T)
  
  res <- as.data.frame((res)[[1]])
  output <- list(MCMC=res,Preynames=rownames(preys))
  class(output) <- 'ind_props'
  return(output)
  
}

PopandIndprops.SI <- function(data,nIter=1000,nBurnin=1000)
{
  n.preys = data$n.preys
  n.preds = data$n.preds
  
  R_SI = data$datas.SI$R.SI
  mean_cs = matrix(unlist(data$datas.SI$mean_cs),n.preys,isos)
  tau_cs = data$datas.SI$tau_cs
  Rnot_SI = data$datas.SI$Rnot.SI
  isos = data$datas.SI$isos
  ni.SI = data$datas.SI$ni.SI
  preds.SI = data$datas.SI$preds.SI
  preym.SI = data$datas.SI$preym.SI
  
  S = diag(eveness,n.preys)
  SS = diag(1,n.preys)
  zeros = rep(0,n.preys)
  
  initials.SI=list(list(
    pnorm=matrix(0,n.preds,n.preys),
    pmean=rep(0,n.preys),
    pprec = diag(1,n.preys), 
    cs=mean_cs,
    prey.means_SI=preym.SI,
    predprec_SI = diag(0.01,isos),
    prey.precs_SI = array(1,c(n.preys,isos,isos))
    
  ))
  
  
  
  datas.SI=list('zeros','S','SS','R_SI','mean_cs','tau_cs','Rnot_SI','n.preys','n.preds','isos','ni.SI','preds.SI','preym.SI')
  
  vars = c('prop','pop.prop')
  
  # compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...

  res <- BRugsFit(system.file("exec","Pop.and.Ind.props.SI.bugs",package = 'FASTIN'), datas.SI, inits=initials.SI, numChains = 1, vars,
                  nBurnin = nBurnin, nIter = Iter, nThin = round(nIter/1000), coda = T,
                  DIC = F, working.directory = getwd(), digits = 4, 
                  BRugsVerbose = T)
  
  res <- as.data.frame((res)[[1]])
  output <- list(MCMC=res)
  class(output) <- 'ind_props'
  return(output)
}

PopandIndprops.combined <- function(data,nIter=1000,nBurnin=1000)
{
  
  n.preys = data$n.preys
  n.preds = data$n.preds
  
  R_SI = data$datas.SI$R.SI
  mean_cs = matrix(unlist(data$datas.SI$mean_cs),n.preys,isos)
  tau_cs = data$datas.SI$tau_cs
  Rnot_SI = data$datas.SI$Rnot.SI
  isos = data$datas.SI$isos
  ni.SI = data$datas.SI$ni.SI
  preds.SI = data$datas.SI$preds.SI
  preym.SI = data$datas.SI$preym.SI
  
  R = data$datas.FA$R
  fc_mean = data$datas.FA$fc_mean
  fc_tau = data$datas.FA$fc_tau
  mean_c =  matrix(unlist(data$datas.FA$mean_c),n.preys,n.fats)
  tau_coeffs = data$datas.FA$tau_c
  Rnot = data$datas.FA$Rnot
  n.fats = data$datas.FA$n.fats
  m.fats = data$datas.FA$m.fats
  ni = data$datas.FA$ni
  preds = data$datas.FA$preds
  preym = data$datas.FA$preym
  
  S = diag(eveness,n.preys)
  SS = diag(1,n.preys)
  zeros=rep(0,n.preys)
  
  initials.comb=list(list(
    pnorm=matrix(0,n.preds,n.preys),
    pmean=rep(0,n.preys),
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
  
  datas.comb=list('zeros','S','SS','R','R_SI','fc_mean','fc_tau','mean_c','tau_coeffs','mean_cs','tau_cs','Rnot','Rnot_SI','n.preys','n.preds','isos','n.fats','m.fats','ni.SI','ni','preds','preds.SI','preym.SI','preym')
  
  vars = c('prop','pop.prop')
  
  # compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...
  res <- BRugsFit(system.file("exec","Pop.and.Ind.props.combined.bugs",package = 'FASTIN'), datas.comb, inits=initials.comb, numChains = 1, vars,
                                   nBurnin = nBurnin, nIter = Iter, nThin = round(nIter/1000), coda = T,
                                   DIC = F, working.directory = getwd(), digits = 4, 
                                   BRugsVerbose = T)
  res <- as.data.frame((res)[[1]])
  output <- list(MCMC=res)
  class(output) <- 'ind_props'
  return(output)
  
}