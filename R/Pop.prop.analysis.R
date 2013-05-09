Poppropanalysis <- function(data,nIter=1000,nBurnin=1000) UseMethod("Poppropanalysis", data)

Poppropanalysis.FA <- function(data,nIter=1000,nBurnin=1000)
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
  preds = matrix(data$datas.FA$preds,n.preds,m.fats)
  preym = data$datas.FA$preym
  
  
  initials=list(list(
    fc = fc_mean,
    fracs=mean_c,             
    ps = rep(1/n.preys,n.preys),
    prey.means=preym,
    predprec = diag(1,m.fats),
    prey.precs = array(1,c(n.preys,m.fats,m.fats))
  ))
  
  datas=list('R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','ni','preds','preym')
  
  vars = c('prop')
  
  # compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...
  res <- BRugsFit(system.file("exec","Pop.prop.analysis.FA.bugs",package = 'FASTIN'), datas, inits=initials, numChains = 1, vars,
                  nBurnin = nBurnin, nIter = nIter, nThin = round(nIter/1000), coda = T,
                  DIC = F, working.directory = getwd(), digits = 4,
                  BRugsVerbose = T)
  
  res <- as.data.frame((res)[[1]])
  output <- list(MCMC=res)
  class(output) <- 'pop_props'
  return(output)
  
}

Poppropanalysis.SI <- function(data,nIter=1000,nBurnin=1000)
{
  n.preys = data$n.preys
  n.preds = data$n.preds
  
  R_SI = data$datas.SI$R.SI
  mean_cs = matrix(unlist(data$datas.SI$mean_cs),n.preys,isos)
  tau_cs = data$datas.SI$tau_cs
  Rnot_SI = data$datas.SI$Rnot.SI
  isos = data$datas.SI$isos
  ni.SI = data$datas.SI$ni.SI
  preds.SI = matrix(unlist(data$datas.SI$preds.SI),n.preds,isos)
  preym.SI = data$datas.SI$preym.SI
  
  
  initials.SI=list(list(
    
    ps = rep(1/n.preys,n.preys),#matrix(1/n.preys,n.preds,n.preys),
    cs=mean_cs,
    prey.means_SI=preym.SI,
    predprec_SI = diag(0.01,isos),
    prey.precs_SI = array(1,c(n.preys,isos,isos))
    
  ))
  
  
  
  datas.SI=list('R_SI','mean_cs','tau_cs','Rnot_SI','n.preys','n.preds','isos','ni.SI','preds.SI','preym.SI')
  
  vars = c('prop')
  
  # compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...
  
  res <- BRugsFit(system.file("exec","Pop.prop.analysis.SI.bugs",package = 'FASTIN'), datas.SI, inits=initials.SI, numChains = 1, vars,
                  nBurnin = nBurnin, nIter = nIter, nThin = round(nIter/1000), coda = T,
                  DIC = F, working.directory = getwd(), digits = 4, 
                  BRugsVerbose = T)
  
  res <- as.data.frame((res)[[1]])
  output <- list(MCMC=res)
  class(output) <- 'pop_props'
  return(output)
}

Poppropanalysis.combined <- function(data,nIter=1000,nBurnin=1000)
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
  mean_c = matrix(unlist(data$datas.FA$mean_c),n.preys,n.fats)
  tau_coeffs = data.matrix(data$datas.FA$tau_c)
  Rnot = data$datas.FA$Rnot
  n.fats = data$datas.FA$n.fats
  m.fats = data$datas.FA$m.fats
  ni = data$datas.FA$ni
  preds = data$datas.FA$preds[,]
  preym = data$datas.FA$preym[,]
  
  initials.comb=list(list(
    fc = fc_mean,
    fracs=mean_c,             
    ps = rep(1/n.preys,n.preys),#matrix(1/n.preys,n.preds,n.preys),
    prey.means=preym,
    predprec = diag(1,m.fats),
    prey.precs = array(0.1,c(n.preys,m.fats,m.fats)),
    cs=mean_cs,
    prey.means_SI=preym.SI,
    predprec_SI = diag(0.01,isos),
    prey.precs_SI = array(1,c(n.preys,isos,isos))
    
  ))
  
  
  datas.comb=list('R','R_SI','fc_mean','fc_tau','mean_c','tau_coeffs','mean_cs','tau_cs','Rnot','Rnot_SI','n.preys','n.preds','isos','n.fats','m.fats','ni.SI','ni','preds','preds.SI','preym.SI','preym')
  
  vars = c('prop')
  
  # compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...
  
  res <- BRugsFit(system.file("exec","Pop.prop.analysis.combined.bugs",package = 'FASTIN'), datas.comb, inits=initials.comb, numChains = 1, vars,
                  nBurnin = nBurnin, nIter = nIter, nThin = round(nIter/1000), coda = T,
                  DIC = F, working.directory = getwd(), digits = 4, 
                  BRugsVerbose = T)
  res <- as.data.frame((res)[[1]])
  output <- list(MCMC=res,Preynames=rownames(preys))
  class(output) <- 'pop_props'
  return(output)
  
}