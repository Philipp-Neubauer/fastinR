addcovs <- function(Groups='',Covariates=''){
  
  if (nchar(Covariates)>0 & nchar(Groups)==0)
  {
    Covs <- read.csv(Covariates,header=T)
    Covs <- cbind(rep(1,nrow(Covs)),Covs)
    n.covs <- ncol(Covs)
    guiSet('Covs',Covs)
  } else if (nchar(Covariates)==0 & nchar(Groups)>0) 
  {
    Grps <- read.csv(Groups,header=T)
    Grp.names <- unlist(unique(Grps)) 
    Grps <- read.csv(Groups,header=T)
    for (i in 1:ncol(Grps)){
      vg <- as.vector(Grps[,i])
      Grps[,i] <- as.factor(vg)
    }
    
    Covs <- model.matrix(attr(model.frame(1:nrow(Grps)~.,data=Grps),'terms'),data=Grps)[,]
    colnames(Covs) <- Grp.names[length(Grp.names):1]
    guiSet('Covs',Covs)
    
  } else if (nchar(Covariates)>0 & nchar(Groups)>0) 
  {
    Covs <- read.csv(Covariates,header=T)
    Grps <- read.csv(Groups,header=T)
    Covnames <- names(Covs)
    Grp.names <- unlist(unique(Grps)) 
    
    for (i in 1:ncol(Grps)){
      vg <- as.vector(Grps[,i])
      Grps[,i] <- as.factor(vg)
    }
    
    Covs <- cbind(model.matrix(attr(model.frame(1:nrow(Grps)~.,data=Grps),'terms'),data=Grps)[,],Covs)
    colnames(Covs) <- c(Grp.names[length(Grp.names):1],Covnames)
    guiSet('Covs',Covs)
  }
}

#### SI -----

addSI <- function(SI.predators=NULL,SI.preys=NULL,Frac.Coeffs.mean='',Frac.Coeffs.var='',FC.mean=1,FC.var=1,R.diag.SI=0.2,datas=NULL){
  
  # check if GUI is being used
  if(exists('GUI',envir=.GlobalEnv)){
    GUI <- get('GUI',envir=.GlobalEnv)
    if (GUI) datas <- guiGetSafe('datas')
  } else {
    GUI=F
  }
  
  ## first check for potential conflicts
  
  stopifnot(nchar(SI.predators)>0 & nchar(SI.preys)>0)
  
  # import predator and prey data - note that the first column is names, or an index
  predators.SI = read.csv(SI.predators,header=T,row.names=1)
  preys.SI = read.csv(SI.preys,header=T)
  
  n.preds <- dim(predators.SI)[1]
  
  preys.ix  <- as.character(preys.SI[,1])
  preys.SI  <- preys.SI[,-1]
  
  # check that samples are in the same order
  if(length(datas)>1) stopifnot(preys.ix==datas$prey.ix)
  
  preys.names  <- unique(preys.ix)
  if(GUI) guiSet('prey.names',preys.names )
  
  # set number of isotopes
  isos=ncol(predators.SI)
  #number of preys species/groups
  n.preys <- length(preys.names)
  
  # improt prior predator variance
  Rnot_SI = diag(R.diag.SI,isos)
  
  # deal with fractionation coeffs
  if ((nchar(Frac.Coeffs.mean)>0 & nchar(Frac.Coeffs.var)==0) | (nchar(Frac.Coeffs.mean)==0 & nchar(Frac.Coeffs.var)>0))
  {
    stop('The mean AND variances of FCs for each isotope need to be supplied')
  } else if (nchar(Frac.Coeffs.mean)>0 & nchar(Frac.Coeffs.var)>0)
  {     
    mean_cs = read.csv(Frac.Coeffs.mean,header=T,row.names=1)
    var_cs  = read.csv(Frac.Coeffs.var,header=T,row.names=1)
    stopifnot(dim(mean_cs)[1]==n.preys & dim(mean_cs)[2]==isos)
    stopifnot(dim(var_cs)[1]==n.preys & dim(var_cs)[2]==isos)
  } else if (nchar(Frac.Coeffs.mean)==0 & nchar(Frac.Coeffs.var)==0)
  {
    mean_cs = matrix(FC.mean,isos,n.preys,byrow=T)
    var_cs = matrix(FC.var,isos,n.preys,byrow=T)
  }
  
  # prey means and vars
  preym.SI <- matrix(,n.preys,isos) 
  # calc prey emans and FC means
  
  for (i in 1:n.preys){
    
    preym.SI[i,] <- apply(preys.SI[preys.ix==unique(preys.ix)[i],],2,mean)
    
  }
  
  # now prepare data for analysis ---
  
  R.SI <- array(,c(isos,isos,n.preys))
  ni.SI<-rep(NA,n.preys)
  for (i in 1:n.preys){
    ni.SI[i] <- max(isos+1,sum(preys.ix==unique(preys.ix)[i])-1)
    R.SI[,,i]=cov(preys.SI[preys.ix==unique(preys.ix)[i],])*ni.SI[i]
  }
  
  datas.SI <- list(isos=isos,R.SI=R.SI,Rnot.SI=Rnot_SI,preys.SI=preys.SI,preym.SI=preym.SI,preds.SI=predators.SI,ni.SI=ni.SI,mean_cs=mean_cs,tau_cs=1/var_cs)
  
  if(length(datas)<=1){
    datas <- list(n.preys = n.preys,n.preds=n.preds,prey.ix=preys.ix,datas.FA=NULL,datas.SI=datas.SI,even=NULL)
  } else {
    datas$datas.SI = datas.SI
    datas$n.preys = n.preys
    datas$n.preds = n.preds
    datas$prey.ix = preys.ix
  }
  
  ifelse(GUI,guiSet('datas',datas),return(datas))
}

####### FAs ----

addFA <- function(FA.predators=NULL,FA.preys=NULL,fat.conts = '',Conv.Coeffs.mean='',Conv.Coeffs.var='',FC.mean=1,FC.var=1,CC.mean=1,CC.var=1,R.diag=0.2,datas=NULL){
  
  # check if GUI is being used
  if(exists('GUI',envir=.GlobalEnv)){
    GUI <- get('GUI',envir=.GlobalEnv)
    if (GUI) datas <- guiGetSafe('datas')
  } else {
    GUI=F
  }
  
  # import predator and prey FA profiles
  predators = read.csv(FA.predators,header=T,row.names=1)
  preys = read.csv(FA.preys,header=T)
  n.preds <- dim(predators)[1]
  preys.ix <- as.character(preys[,1])
  
  if(length(datas)>1) stopifnot(preys.ix==datas$prey.ix)
  
  preys.names <- as.character(unique(preys.ix))
  
  if (GUI) guiSet('prey.names',preys.names )
  
  preys = preys[,-1]
  
  n.fats = ncol(predators)
  m.fats = n.fats-1
  
  #number of preys species/groups
  n.preys <- length(unique(preys.ix))
  
  # treat conversion coeffs
  if (nchar(Conv.Coeffs.mean)>0 & nchar(Conv.Coeffs.var)>0)
  {     
    mean_c = read.csv(Conv.Coeffs.mean,header=T,row.names=1)
    var_c  = read.csv(Conv.Coeffs.var,header=T,row.names=1)
    stopifnot(dim(mean_c)[1]==n.preys & dim(mean_c)[2]==n.fats)
    stopifnot(dim(var_c)[1]==n.preys & dim(var_c)[2]==n.fats)
  } else if (nchar(Conv.Coeffs.mean)==0 & nchar(Conv.Coeffs.var)==0)
  {
    mean_c = matrix(CC.mean,n.preys,n.fats,byrow=T)
    var_c =matrix(CC.var,n.preys,n.fats,byrow=T)
  } else
  {
    print('Known conversion coefficients, or a mean AND variance for conversion coefficients need to be supplied')
  }
  
  # deal with fat content
  if(nchar(fat.conts)==0) 
  {
    if(length(FC.mean) == n.preys & length(FC.var) == n.preys)
    {
      fc.mean <- FC.mean; fc.var <- FC.var
    } else if(length(FC.mean) == 1 & length(FC.var) == 1){
      fc.mean <- rep(FC.mean,n.preys); fc.var <- rep(FC.var,n.preys)
    } else {stop('Fat content mean and variance need to be either a single number, or supplied as a vector of length equal to the number of prey items - use R c() notation in that case. In the latter case, or for individual sample fat content please supply a file')}
    
  } else
  {
    fat.cont <- read.csv(fat.conts,header=F)
    if (dim(fat.cont)[2]>2){
      fat.cont <- read.csv(fat.conts,header=F,row.names=1) 
      fc.mean <- fat.cont[,1];fc.var <- fat.cont[,2]
      if (any(fc.mean>1)){
        fc.mean <- fc.mean/100
        fc.var <- fc.var/100
      }
    } else if (any(fat.cont>1)){
      fat.cont <- fat.cont/100
      fc.mean <- tapply(fat.cont,preys.ix,mean)
      fc.var <- tapply(fat.cont,preys.ix,var)
    }
  }
  
  # make sure everything sums to 1
  
  predators <- clo(predators)
  preys <- clo(preys)
  
  # get prey means
  mprey <- aggregate(preys,list(preys.ix),gmean)[,2:(n.fats+1)]
  
  preym <- unclass(alr(mprey))
  preds <- unclass(alr(predators))    
  
  # now prepare data for analysis
  
  R <- array(,c(m.fats,m.fats,n.preys))
  ni<-rep(NA,n.preys)
  for (i in 1:n.preys){
    ni[i] <- max(n.fats+1,sum(preys.ix==unique(preys.ix)[i])-1)
    R[,,i]=cov(alr(preys[preys.ix==unique(preys.ix)[i],]))*ni[i]
  }
  
  ## first some data and inits ----
  
  # set uninformative prior SS matrix for wishart prior alr transformed predator data
  
  Rnot =diag(R.diag,m.fats)
  
  datas.FA <- list(fc_mean=fc.mean,fc_tau=1/fc.var,n.fats=n.fats,m.fats=m.fats,R=R,Rnot=Rnot,preys=preys,preds.FA=predators,preym=preym,preds=preds,ni=ni,mean_c=mean_c,tau_c=1/var_c)
  
  if(length(datas)<=1){
    datas <- list(n.preys = n.preys,n.preds=n.preds,prey.ix=preys.ix,datas.FA=datas.FA,datas.SI=NULL,even=NULL)
  } else {
    datas$datas.FA = datas.FA
    datas$n.preys = n.preys
    datas$n.preds=n.preds
    datas$prey.ix=preys.ix
  }
  ifelse(GUI,guiSet('datas',datas),return(datas))
}

### dataplot ---

dataplot <- function(datas=NULL){
  
  # check if GUI is being used
  if(exists('GUI',envir=.GlobalEnv)){
    GUI <- get('GUI',envir=.GlobalEnv)
    if(GUI) datas <- guiGetSafe('datas')
  } else {
    GUI=F
  }
  
  preya <- {}
  preda <- {}
  if(all(!is.na(datas$datas.FA$preys)) & !is.null(datas$datas.FA$preys)) {
      preya=cbind(preya,clr(datas$datas.FA$preys))
      preda=cbind(preda,clr(datas$datas.FA$preds.FA))     
  }
  if(all(!is.na(datas$datas.SI$preys.SI)) & !is.null(datas$datas.SI$preys.SI)) {
      preya=cbind(preya,as.matrix(datas$datas.SI$preys.SI))
      preda=cbind(preda,as.matrix(t(t(datas$datas.SI$preds.SI)-colMeans(datas$datas.SI$mean_cs))))
  }
  
  names(preya)  <- names(preda)
  
  dista <- dist(rbind(preya,preda))
  mds <- metaMDS(dista)
  
  X11()
  pl <- plot(mds,type='n')
  points(pl,'sites',pch=cbind(as.numeric(as.factor(datas$prey.ix)),rep(16,datas$n.preds)),col=cbind(1+as.numeric(as.factor(datas$prey.ix)),rep(1,datas$n.preds)))
  legend('bottomright',c('Predators',unique(datas$prey.ix)),xpd=T,pch=c(16,1:datas$n.preys),col=c(1,2:(datas$n.preys+1)))
  
}