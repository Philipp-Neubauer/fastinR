#' Add Covariates and Group indices for predators
#'     
#' Files must be in .csv format
#'             
#' @param Groups Index of group membership for each predator, one (named) column per grouping variable
#' @param Covariates Covariate values for each predator, one (named) column per covariate
#' @details Use \code{\link{simulation}} to simulate and write these files to inspect the file structure.
#' @seealso \code{\link{addFA}},\code{\link{addSI}},\code{\link{run_MCMC}},\code{\link{simulation}}
#' @author Philipp Neubauer
#' @references Neubauer,.P. and Jensen, O.P. (in prep)
#' @export
#' @examples 
#' Cov_path <- system.file("extdata", "Simdata_Covariates.csv", package="FASTIN")
#' Group_path<- system.file("extdata", "Simdata_Groups.csv", package="FASTIN")
#' addCovs(Groups=Group_path,Covariates=Cov_path)

addCovs <- function(Groups='',Covariates=''){
  
  # check if GUI is being used
  if(exists('GUI',envir=.GlobalEnv)){
    GUI <- get('GUI',envir=.GlobalEnv)
    } else {
    GUI=F
  }
  
  if (nchar(Covariates)>0 & nchar(Groups)==0)
  {
    Covs <- read.csv(Covariates,header=T)
    Covs <- cbind(rep(1,nrow(Covs)),Covs)
    n.covs <- ncol(Covs)
    if(GUI) guiSet('Covs',Covs)
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
    if(GUI) guiSet('Covs',Covs)
    
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
    if(GUI) guiSet('Covs',Covs)
  }
  if(!GUI) return(Covs)
}

#' Add Stable Isotope data for predators and prey items
#'     
#' Files must be in .csv format.
#'             
#' @param SI.predators Predator index/names (first column) and Stable Isotopes (1 row pey predator), with Stable Isotope named across the first row 
#' @param SI.preys Prey names/sample id (first column) and fatty acid profiles (1 row pey prey item), ith Stable Isotope named across the first row 
#' @param Frac.Coeffs.mean Prey specific additive fractionation coefficient means: Prey names (first column) and an n x P matrix for n preys and P Stable Isotopes
#' @param Frac.Coeffs.var Prey specific Fractionation coefficient variances, dimensions as for the means
#' @param FC.mean optional - if no prey specific fractionation coefficiants are supplied via Frac.Coeffs.mean, FC mean can provide either a global (single) mean coefficient or fatty acid specific mean coefficients using R's c(FA_1,FA_2,...) notation for ALL fatty acids.
#' @param FC.var optional - if no prey specific fractionation coefficiants are supplied via Frac.Coeffs.mean, FC var can provide either a global (single) coefficient variance or fatty acid specific coefficient variances using R's c(FA_1,FA_2,...) notation for ALL fatty acids.
#' @param datas a data structure as produced by \code{\link{addSI}}, needed if fatty acids and stable isotopes are added sequentially.
#' @details Use \code{\link{simulation}} to simulate and write these files to inspect the file structure.
#' @seealso \code{\link{addFA}},\code{\link{addCovs}},\code{\link{run_MCMC}},\code{\link{simulation}}
#' @author Philipp Neubauer
#' @references Neubauer,.P. and Jensen, O.P. (in prep)
#' @examples 
#' SI.predators <- system.file("extdata", "Simdata_SI_preds.csv", package="FASTIN")
#' SI.preys <- system.file("extdata", "Simdata_SI_preys.csv", package="FASTIN")
#' Frac.Coeffs.mean <- system.file("extdata", "Simdata_SI_fc_means.csv", package="FASTIN")
#' Frac.Coeffs.var <- system.file("extdata", "Simdata_SI_fc_var.csv", package="FASTIN")
#' dats <- addSI(SI.predators=SI.predators,SI.preys=SI.preys,Frac.Coeffs.mean=Frac.Coeffs.mean,Frac.Coeffs.var=Frac.Coeffs.var)
#' @export
addSI <- function(SI.predators=NULL,SI.preys=NULL,Frac.Coeffs.mean='',Frac.Coeffs.var='',FC.mean=1,FC.var=1,datas=NULL){
  
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
  
  datas.SI <- list(isos=isos,R.SI=R.SI,Rnot.SI=NULL,preys.SI=preys.SI,preym.SI=preym.SI,preds.SI=predators.SI,ni.SI=ni.SI,mean_cs=mean_cs,tau_cs=1/var_cs)
  
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

#' Add Fatty Acid profile data for predators and prey items
#'     
#' Files must be in .csv format.
#'             
#' @param FA.predators Predator index/names (first column) and fatty acid profiles (1 row pey predator), with fatty acids named across the first row 
#' @param FA.preys Prey names/sample id (first column) and fatty acid profiles (1 row pey prey item), with fatty acids names across the first row
#' @param fat.conts Prey fat contents, as (columnwise) mean and variance per prey species or specified for each prey sample for the main analysis, in that case the first column is the prey sample id id and the second column is the individual sample's fat content
#' @param Conv.Coeffs.mean Prey specific conversion coefficient means: Prey names (first column) and an n x P matrix for n preys and P fatty acids
#' @param Conv.Coeffs.var Prey specific conversion coefficient variances, dimensions as for the means
#' @param FC.mean optional - if no prey or sample specific fat content means are supplied in a fat.conts file, prey specific coefficients can be entered here using R's c(FC_1,FC_2,...) notation.
#' @param FC.var optional - if no prey or sample specific fat content variances are supplied in a fat.conts file, prey specific coefficients can be entered here using R's c(FC_1,FC_2,...) notation.
#' @param CC.mean optional - if no prey specific fractionation coefficiants are supplied via Conv.Coeffs.mean, CC.mean can provide either a global (single) mean coefficient or fatty acid specific mean coefficients using R's c(FA_1,FA_2,...) notation for ALL fatty acids.
#' @param CC.var optional - if no prey specific fractionation coefficiants are supplied via Conv.Coeffs.mean, CC.var can provide either a global (single) coefficient variance or fatty acid specific coefficient variances using R's c(FA_1,FA_2,...) notation for ALL fatty acids.
#' @param datas a data structure as produced by \code{\link{addSI}}, needed if fatty acids and stable isotopes are added sequentially.
#' @details Use \code{\link{simulation}} to simulate and write these files to inspect the file structure.
#' @seealso \code{\link{addSI}},\code{\link{addCovs}},\code{\link{selectvars}},\code{\link{run_MCMC}},\code{\link{simulation}}
#' @author Philipp Neubauer
#' @references Neubauer,.P. and Jensen, O.P. (in prep)
#' @examples 
#' FA.predators <- system.file("extdata", "Simdata_FA_preds.csv", package="FASTIN")
#' FA.preys <- system.file("extdata", "Simdata_FA_preys.csv", package="FASTIN")
#' Conv.Coeffs.mean <- system.file("extdata", "Simdata_FA_cc_means.csv", package="FASTIN")
#' Conv.Coeffs.var <- system.file("extdata", "Simdata_FA_cc_var.csv", package="FASTIN")
#' fat.conts <- system.file("extdata", "Simdata_fat_cont.csv", package="FASTIN")
#' dats <- addFA(FA.predators=FA.predators,FA.preys=FA.preys,fat.conts=fat.conts,Conv.Coeffs.mean=Conv.Coeffs.mean,Conv.Coeffs.var=Conv.Coeffs.var)
#' @export
addFA <- function(FA.predators=NULL,FA.preys=NULL,fat.conts = '',Conv.Coeffs.mean='',Conv.Coeffs.var='',FC.mean=1,FC.var=1,CC.mean=1,CC.var=1,datas=NULL){
  
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
  
 datas.FA <- list(fc_mean=fc.mean,fc_tau=1/fc.var,n.fats=n.fats,m.fats=m.fats,R=R,Rnot=NULL,preys=preys,preds.FA=predators,preym=preym,preds=preds,ni=ni,mean_c=mean_c,tau_c=1/var_c)
  
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