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