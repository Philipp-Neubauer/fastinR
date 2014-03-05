#' Add Covariates and Group indices for predators
#'     
#' Files must be in .csv format
#'             
#' @param Groups Index of group membership for each predator, one (named) column per grouping variable
#' @param Covariates Covariate values for each predator, one (named) column per covariate
#' @details Use \code{\link{simulation}} to simulate and write these files to inspect the file structure.
#' @seealso \code{\link{add_FA}},\code{\link{add_SI}},\code{\link{run_MCMC}},\code{\link{simulation}}
#' @author Philipp Neubauer
#' @references Neubauer,.P. and Jensen, O.P. (in prep)
#' @export
#' @examples 
#' Cov_path <- system.file("extdata", "Simdata_Covariates.csv", package="fastinR")
#' Group_path<- system.file("extdata", "Simdata_Groups.csv", package="fastinR")
#' add_Covs(Groups=Group_path,Covariates=Cov_path)

add_Covs <- function(Groups='',Covariates=''){
  
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
#' @param datas a data structure as produced by \code{\link{add_SI}}, needed if fatty acids and stable isotopes are added sequentially.
#' @details Use \code{\link{simulation}} to simulate and write these files to inspect the file structure.
#' @seealso \code{\link{add_FA}},\code{\link{add_Covs}},\code{\link{run_MCMC}},\code{\link{simulation}}
#' @author Philipp Neubauer
#' @references Neubauer,.P. and Jensen, O.P. (in prep)
#' @examples 
#' SI.predators <- system.file("extdata", "Simdata_SI_preds.csv", package="fastinR")
#' SI.preys <- system.file("extdata", "Simdata_SI_preys.csv", package="fastinR")
#' Frac.Coeffs.mean <- system.file("extdata", "Simdata_SI_fc_means.csv", package="fastinR")
#' Frac.Coeffs.var <- system.file("extdata", "Simdata_SI_fc_var.csv", package="fastinR")
#' dats <- add_SI(SI.predators=SI.predators,SI.preys=SI.preys,Frac.Coeffs.mean=Frac.Coeffs.mean,Frac.Coeffs.var=Frac.Coeffs.var)
#' @export
add_SI <- function(SI.predators=NULL,SI.preys=NULL,Frac.Coeffs.mean='',Frac.Coeffs.var='',FC.mean=1,FC.var=1,datas=NULL){
  
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
  if(is.character(SI.predators)) {
    predators.SI = read.csv(SI.predators,header=T,row.names=1)
  } else {
    predators.SI = SI.predators
  }
  if(is.character(SI.preys)) {
    preys.SI = read.csv(SI.preys,header=T)
  } else {
    preys.SI = SI.preys
  }  
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
#' @param datas a data structure as produced by \code{\link{add_SI}}, needed if fatty acids and stable isotopes are added sequentially.
#' @param LN.par - are fat content means and variances given as log-normal parameters or sample mean and variance? 
#' @details Use \code{\link{simulation}} to simulate and write these files to inspect the file structure.
#' @seealso \code{\link{add_SI}},\code{\link{add_Covs}},\code{\link{select_vars}},\code{\link{run_MCMC}},\code{\link{simulation}}
#' @author Philipp Neubauer
#' @references Neubauer,.P. and Jensen, O.P. (in prep)
#' @examples 
#' FA.predators <- system.file("extdata", "Simdata_FA_preds.csv", package="fastinR")
#' FA.preys <- system.file("extdata", "Simdata_FA_preys.csv", package="fastinR")
#' Conv.Coeffs.mean <- system.file("extdata", "Simdata_FA_cc_means.csv", package="fastinR")
#' Conv.Coeffs.var <- system.file("extdata", "Simdata_FA_cc_var.csv", package="fastinR")
#' fat.conts <- system.file("extdata", "Simdata_fat_cont.csv", package="fastinR")
#' dats <- add_FA(FA.predators=FA.predators,FA.preys=FA.preys,fat.conts=fat.conts,Conv.Coeffs.mean=Conv.Coeffs.mean,Conv.Coeffs.var=Conv.Coeffs.var)
#' @export
add_FA <- function(FA.predators=NULL,FA.preys=NULL,fat.conts = '',Conv.Coeffs.mean='',Conv.Coeffs.var='',FC.mean=1,FC.var=1,CC.mean=1,CC.var=1,datas=NULL,LN.par=F){
  
  # check if GUI is being used
  if(exists('GUI',envir=.GlobalEnv)){
    GUI <- get('GUI',envir=.GlobalEnv)
    if (GUI) datas <- guiGetSafe('datas')
  } else {
    GUI=F
  }
  
  # import predator and prey FA profiles
  if (is.character(FA.predators)) {
    predators = read.csv(FA.predators,header=T,row.names=1)
  } else {predators = FA.predators}
  
  if (is.character(FA.preys)){
    preys = read.csv(FA.preys,header=T)
  } else {
    preys = FA.preys
  }
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
    if(length(mean_c)==n.fats){
      mean_c = matrix(CC.mean,n.preys,n.fats,byrow=T)
      var_c =matrix(CC.var,n.preys,n.fats,byrow=T)
    } else {
      mean_c = matrix(CC.mean,n.preys,n.fats)
      var_c =matrix(CC.var,n.preys,n.fats)
    }
    
  } else
  {
    print('Known conversion coefficients, or a mean AND variance for conversion coefficients need to be supplied')
  }
  
#   # make sure that cs sum to one
#   if(any(rowSums(mean_c))!=1){
#     mean_c <- t(apply(mean_c,1,function(x) x/sum(x)))
#     sums = (apply(mean_c,1,function(x) sum(x)))
#     var_c = var_c/(sums^2)
#   }
  
  # covert to gamma parameters
  rate <- mean_c/var_c
  shape <- mean_c^2/var_c
  
  mean_c <- shape
  var_c <- rate
  
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
      
    } else if (dim(fat.cont)[2]==2){
      fc.mean <- fat.cont[,1];fc.var <- fat.cont[,2]
    } else {
      fc.mean <- tapply(fat.cont,preys.ix,mean)
      fc.var <- tapply(fat.cont,preys.ix,var)
    }
  }
  
  if(LN.par == F){
    fc.var = log(fc.var + fc.mean^2) - 2*log(fc.mean)
    fc.mean = log(fc.mean)-fc.var/2
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