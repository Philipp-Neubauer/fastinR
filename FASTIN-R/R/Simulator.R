#' @name simulation
#' @title GUI: Simulate Fatty Acid profiles and Stable Isotope data
#' @usage simulation()
#' @details The function calls a gui to simulate data, functions are of little use outside the gui facility.
#'
#' Simulating data proceeds by selecting to simualte either individual diet proportions for each predator as drawn from a population level distribution for diet proportions, or by simulating grouped diet proportions (e.g., size groups, geographic groups) and/or by letting covariates influence the drawn diet proportions (e.g., size based diet shifts).
#' 
#' Pressing \code{Simulate Fatty Acid data} or \code{Simulate Stable Isotope data} makes sense only after proportions were simulated - converseley, when updating proportions, FA and/or SI data need to be re-simulated to have the new diet proportions alter the predator fatty acid/SI makeup.
#' 
#' \code{Plot current simulation} will draw a NMDS (non-metric multidimensional scaling) plot to show the simulated data.
#' 
#' The gui will remain open to allow for many tries at simulating data. Variuous configurations can be written to file to allow for exploration of the model with \code{fastin()}.
#' @return \code{Write simulations to files} will produce a series of files that can be used as inputs to the FASTIN gui or individual functions. The prefix of the files is entered by the user, the remainder of the filename suggests the contents:
#' \item{*_props}{Simualted diet proportions}
#' 
#'     For Fatty Acids:
#'     
#'             \item{*_FA_preys.csv}{Prey names (first column) and fatty acid profiles, with fatty acids names across the first row}
#'             \item{*_FA_preds.csv}{Predator index (first column) and fatty acid profiles, with fatty acids names across the first row}
#'             \item{*_fat_cont.csv}{Prey fat contents (here as mean and variance, can be      specified for each prey sample for the main analysis, in that case the                                          first column is the prey id and the second column is the individual                                           sample's fat content)}
#'             \item{*_FA_cc_means.csv}{Prey specific conversion coefficient means:                                         Prey names (first column) and an n x P matrix for n preys and P fatty acids} 
#'             \item{*_FA_cc_var.csv}{Prey specific conversion coefficient variances, dimensions as for the means}
#'             For Stable Isotopes:
#'             \item{*_SI_preys.csv}{Prey names (first column) and stable isotope values, with SI names across the first row}
#'             \item{*_SI_preds.csv}{Predator index (first column) and stable isotope values, with SI names across the first row}
#'             \item{*_SI_fc_means.csv}{Prey specific SI additive fractionation coefficient means}
#'             \item{*_SI_fc_var.csv}{Prey specific additive fractionation coefficient variance, dimensions as for the means}
#'             For Covariates for predator proportions and grouped predators:
#'             
#'             \item{*_Cov_n_Grp_effects.csv}{an n*p matrix of group and covariate                                   influences on diet proportions for n preys and p (groups+covariates)}
#'             \item{*_Groups.csv}{Group membership for each predator}
#'             \item{*_Covariates.csv}{Covariate values for each predator}
#' @seealso \code{\link{fastin}},\code{\link{addFA}},\code{\link{addSI}},\code{\link{run_MCMC}}
#'             
#' @references  Neubauer.P. and Jensen, O.P. (in prep)
#' @author Philipp Neubauer
#' @examples
#'\dontrun{simulation()}
#' @export                 
simulation <- function(){
  
  ######### sim functions
#  require(fgui)
  options(warn = -2)
  
  simulator <- function(n.preys=3,n.preds=10,nsamples=30,simProps=NULL,simGroups=NULL,simCovs=NULL,FAsim.data=NULL,SIsim.data=NULL,simplot=NULL,simwrite=NULL){
    
    if( is.na(guiGetSafe("n.preys")) |guiGetSafe("n.preys") != n.preys | guiGetSafe("n.preds") != n.preds | guiGetSafe("nsamples") != nsamples){
      
      preys.ix={};
      for (i in 1:n.preys){
        preys.ix <- c(preys.ix,rep(paste('Prey',i,sep='_'),nsamples))
        
      }
      
      preys.ix.num = as.numeric(as.factor(preys.ix))
      
      guiSet("n.preys",n.preys)
      guiSet("n.preds",n.preds)
      guiSet("nsamples",nsamples)
      guiSet("preys.ix",preys.ix)
      
      switcheroo <- guiGetSafe("gswitch")
      
      if (!is.na(switcheroo)){
        
        if (switcheroo == 'gprops')
        {
          simProps(guiGetSafe("eveness"))       
        } else if (switcheroo == 'Covas')
        {
          simCovs(guiGetSafe("n.covs"),guiGetSafe("covariate.effect.size"))       
        } else if (switcheroo == 'groups')
        {
          simGroups(guiGetSafe("n.groups"),guiGetSafe("group.effect.size"))
        }  else if (switcheroo == 'combined')
        { guiSet("gswitch",NA) # rest switch to make sure no wrong matrices are concatenated
          simGroups(guiGetSafe("n.groups"),guiGetSafe("group.effect.size"))
          simCovs(guiGetSafe("n.covs"),guiGetSafe("covariate.effect.size"))       
        }      
        
      }
    }
   }
  
  simProps <- function(eveness=1){
    
    n.preys <- guiGetSafe("n.preys")
    n.preds <- guiGetSafe("n.preds")
    
    #require(MCMCpack)
    
    v1 = 0.5; v2 = 0.7
    Q = rlnorm(n.preys,c(rep(0,n.preys-1),eveness),v2)
    #rlnorm(n.preys,c(rep(0,n.preys-1),eveness),v2)# individual proportions
    props=MCMCpack::rdirichlet(n.preds,Q)
    
    guiSet("props",props)
    guiSet("gswitch",'gprops')
    guiSet("eveness",eveness)
    
    if(all(!is.na(guiGetSafe("preds"))) | all(!is.na(guiGetSafe("preds.SI"))))
    {
      tkmessageBox( message="simulating new diet signatures", title="New FA/SI signatures" )
    }
    if(all(!is.na(guiGetSafe("preds")))){
      sim.FA(guiGetSafe("sep.FA"),guiGetSafe("n.fats"),guiGetSafe("cvar"))
    }
    if(all(!is.na(guiGetSafe("preds.SI")))){
      sim.SI(guiGetSafe("sep.SI"),guiGetSafe("isos"))
    }
    
  }
  
  simGroups <- function(n.groups=2,group.effect.size=1){
    
    n.preys <- guiGetSafe("n.preys")
    n.preds <- guiGetSafe("n.preds")
    
    switchin <- guiGetSafe("gswitch")
   
    if (!is.na(switchin)){
      # if I only have Covs so far, get them and concatenate
      if(switchin == 'Covas'){
        Covs    <- guiGetSafe("Covs")
        beta    <- guiGetSafe("beta")
      } else if (switchin == 'combined') {
        # if I already have combined data, only take the covs and discard prior groups
        n.covs.old <- guiGetSafe("n.covs")        
        Covs    <- guiGetSafe("Covs")
        beta    <- guiGetSafe("beta")[,c(1,(ncol(Covs)-(n.covs.old-1)):ncol(Covs))]
        Covs <- as.data.frame(Covs[,c(1,(ncol(Covs)-(n.covs.old-1)):ncol(Covs))])
                
      } else {(Covs=NA)}
    } else {(Covs=NA)}
    
    Grps <- as.data.frame(as.factor(sample(1:n.groups,n.preds,replace=T)))
    colnames(Grps) <- 'Grouping_Variable'
    
    Covs.new<-as.data.frame(model.matrix(attr(model.frame(1:nrow(Grps)~.,data=Grps),'terms'),data=Grps)[,])
      
    n.covs <- ncol(Covs.new)
    
    names(Covs.new)[1] <- 'mean'
    for (i in 1:n.covs){
      names(Covs.new)[i] <- paste('Group',i)
    }
    
    beta.new <-matrix(,n.preys,n.covs)
       
    for (j in 1:n.covs)
    {
      beta.new[1:n.preys,j] <- mvrnorm(1,rep(0,n.preys),diag(group.effect.size,n.preys))
    }
    
    colnames(beta.new) <- names(Covs.new)
        
    if(all(is.na(Covs))){
      guiSet("gswitch",'groups')
      Covs <- Covs.new
      beta <- beta.new
    } else {
      guiSet("gswitch",'combined')
      Covs <- cbind(Covs.new,Covs[,2:ncol(Covs)])
      beta <- cbind(beta.new,beta[,2:ncol(beta)])
    }
    
    n.covs <- ncol(beta)
    
    mus <- matrix(,n.preds,n.preys)
    pnorm <- mus
    props <- mus
    
    for(i in 1:n.preds) {
      for(j in 1:n.preys) {
        mus[i,j] <-beta[j,1:n.covs] %*% as.numeric(Covs[i,1:n.covs])
      }
      pnorm[i,1:n.preys] <- mvrnorm(1,mus[i,],diag(0.05,n.preys))
      p_unn <- exp(pnorm)
      for(j in 1:n.preys) {    
        
        #closure
        props[i,j] <- (p_unn[i,j])/sum(p_unn[i,1:n.preys])
        
      }  
    }
    
    guiSet("beta",beta)
    guiSet("props",props)
    guiSet("Covs",Covs)
    guiSet("Grps",Grps)
    guiSet("n.groups",n.groups)
    guiSet("group.effect.size",group.effect.size)
    
    if(all(!is.na(guiGetSafe("preds"))) | all(!is.na(guiGetSafe("preds.SI"))))
    {
      tkmessageBox( message="simulating new diet signatures", title="New FA/SI signatures" )
    }
    if(all(!is.na(guiGetSafe("preds")))){
      sim.FA(guiGetSafe("sep.FA"),guiGetSafe("n.fats"),guiGetSafe("cvar"))
    }
    if(all(!is.na(guiGetSafe("preds.SI")))){
      sim.SI(guiGetSafe("sep.SI"),guiGetSafe("isos"))
    }
    
  }
  
  simCovs <- function(n.covs=1,covariate.effect.size=1){
       
       
    n.preys <- guiGetSafe("n.preys")
    n.preds <- guiGetSafe("n.preds")
    
    switchin <- guiGetSafe("gswitch")
    
    if (!is.na(switchin)){
      # if I only have groups so far, get them and concatenate
      if(switchin == 'groups'){
        Covs    <- guiGetSafe("Covs")
        beta    <- guiGetSafe("beta")
      } else if (switchin == 'combined') {
        # if I already have combined data, only take the groups and discard prior Covs
        n.covs.old    <- guiGetSafe("n.covs")
        Covs    <- guiGetSafe("Covs")
        beta    <- guiGetSafe("beta")[,1:(ncol(Covs)-n.covs.old)]
        Covs <- as.data.frame(Covs[,1:(ncol(Covs)-n.covs.old)])       
      } else {(Covs=NA)}
    } else {(Covs=NA)}
    
    
    guiSet("n.covs",n.covs)
    n.covs=n.covs+1
 
    beta.new <-matrix(,n.preys,n.covs)
    
    Covs.new <- as.data.frame(matrix(,n.preds,n.covs))
    
    Covs.new[,1] <- 1
    names(Covs.new)[1] <- 'mean'
    
    for (i in 2:(n.covs))
    {
      Covs.new[,i]<-rnorm(n.preds)
      names(Covs.new)[i] <- paste('Covariate',i-1)
    }
    
    for (j in 1:(n.covs))
    {
      beta.new[1:n.preys,j] <- mvrnorm(1,rep(0,n.preys),diag(covariate.effect.size,n.preys))
    }
    
    colnames(beta.new) <- names(Covs.new)
    
    if(all(is.na(Covs))){
      guiSet("gswitch",'Covas')
      Covs <- Covs.new
      beta <- beta.new
    } else {
      guiSet("gswitch",'combined')
      Covs <- cbind(Covs,Covs.new[,2:ncol(Covs.new)])
      beta <- cbind(beta,beta.new[,2:ncol(Covs.new)])
    }
    
    n.covs <- ncol(beta)
    
    mus <- matrix(,n.preds,n.preys)
    pnorm <- mus
    props <- mus
    
    for(i in 1:n.preds) {
      for(j in 1:n.preys) {
        mus[i,j] <-beta[j,1:n.covs] %*% as.numeric(Covs[i,1:n.covs])
      }
      pnorm[i,1:n.preys] <- mvrnorm(1,mus[i,],diag(0.05,n.preys))
      p_unn <- exp(pnorm)
      for(j in 1:n.preys) {    
        
        #closure
        props[i,j] <- (p_unn[i,j])/sum(p_unn[i,1:n.preys])
        
      }  
    }
    
    guiSet("beta",beta)
    guiSet("props",props)
    guiSet("Covs",Covs)
    guiSet("Covariates",Covs.new[,2:ncol(Covs.new)])
    
    guiSet("covariate.effect.size",covariate.effect.size)
    
    if(all(!is.na(guiGetSafe("preds"))) | all(!is.na(guiGetSafe("preds.SI"))))
    {
      tkmessageBox( message="simulating new diet signatures", title="New FA/SI signatures" )
    }
    if(all(!is.na(guiGetSafe("preds")))){
      sim.FA(guiGetSafe("sep.FA"),guiGetSafe("n.fats"),guiGetSafe("cvar"))
    }
    if(all(!is.na(guiGetSafe("preds.SI")))){
      sim.SI(guiGetSafe("sep.SI"),guiGetSafe("isos"))
    }
  }
  
  sim.FA <- function(sep=1,n.fats=10,cvar=0.5){
    
   # require(MCMCpack)
    
    props <- guiGetSafe("props")
    
    if (any(is.na(props)))
    {
      stop("Please simulate proportions first to initialize parameters")
    }
    
    n.preys <- guiGetSafe("n.preys")
    n.preds <- guiGetSafe("n.preds")
    nsamples <- guiGetSafe("nsamples")
    preys.ix <- guiGetSafe("preys.ix")
    
    As = matrix(rlnorm(n.preys*n.fats,rnorm(n.preys*n.fats,sep)),n.preys,n.fats,byrow=T)
    
    preys <- array(,c(nsamples[1],n.fats,n.preys))
    for (i in 1:n.preys){
      preys[,,i] <- MCMCpack::rdirichlet(nsamples,As[i,])
    }
    
    preya={}
    for (i in 1:n.preys){
      preya<- rbind(preya,preys[,,i])
    }

    colnames(preya) <- sprintf(paste0("FA_%d"),1:n.fats)
    
    #randomly draw fat content
    fc_mean <- as.matrix(sample(40,n.preys))/10+1
    rownames(fc_mean) <- unique(preys.ix)
    # set variance of fat content to 0.4 - this is loosely based on literature
    fc_sd <- matrix(0.5,n.preys,1)
    
    # draw mean conversion coeffs - cvar is the variance of mean conversion coeffs around 1; 0.2 seems reasonable  (c.f., Rosen & Tollit histogram)
    mean_css=data.frame(matrix(rlnorm(n.fats*n.preys,0,cvar),n.preys,n.fats))
    rownames(mean_css) <- unique(preys.ix)
    
    sd_css=data.frame(matrix(0.05,n.preys,n.fats))
    rownames(sd_css) <- unique(preys.ix)
    
    mprey <-  clo(t(apply(preys,c(3),function(x){exp(colMeans(log(x)))})))
    preds <-  clo(props%*%as.matrix(as.vector(fc_mean)*mprey*mean_css))
    
    
    guiSet("preds",preds)
    guiSet("preys",preya)
    
    
    guiSet("mean_c",mean_css)
    guiSet("sd_c",sd_css)
    
    guiSet("fc_mean",fc_mean)
    guiSet("fc_sd",fc_sd)
    
    # set for re-sim
    guiSet("sep.FA",sep)
    guiSet("n.fats",n.fats)
    guiSet("cvar",cvar)
  }
  
  sim.SI <- function(sep=1,isos=2){
    
    props <- guiGetSafe("props")
    
    if (any(is.na(props)))
    {
      stop("Please click 'Update Parameters' first to initialize parameters")
    }
    
    n.preys <- guiGetSafe("n.preys")
    n.preds <- guiGetSafe("n.preds")
    nsamples <- guiGetSafe("nsamples")
    preys.ix <- guiGetSafe("preys.ix")
    
    
    prey.m <- mvrnorm(n.preys,c(-16,18,rnorm(isos-2,0,sep)),diag(sep,isos))
    prey.cov <- diag(0.2/sep,isos)
    
    preys.SI <- array(NA,c(nsamples,isos,n.preys))
    for (i in 1:n.preys){
      preys.SI[,,i] <-mvrnorm(nsamples,prey.m[i,],prey.cov)
    }
    
    preya.SI={}
    for (i in 1:n.preys){
      preya.SI<- rbind(preya.SI,preys.SI[,,i])
    }
    
    mean_cs = data.frame(matrix(rnorm(isos*n.preys,c(3,1,rnorm(isos-2))),n.preys,isos,byrow=T))
    rownames(mean_cs) <- unique(preys.ix)
    
    sd_cs = data.frame(matrix(rep(0.2,isos),n.preys,isos,byrow=T)) # 1.61
    rownames(sd_cs) <- unique(preys.ix)
    
    # the mixture in matrix algebra (just a matrix product)
    preym.SI <-  t(apply(preys.SI,3,colMeans))
    preds.SI <-  props %*% as.matrix(preym.SI+mean_cs)
    
    
    guiSet("preds.SI",preds.SI)
    guiSet("preys.SI",preya.SI)
    
    guiSet("mean_cs",mean_cs)
    guiSet("sd_cs",sd_cs)
    
    # set for re-sim
    guiSet("sep.SI",sep)
    guiSet("isos",isos)

    
  }    
  
  simwrite <- function(filename='Simdata'){
    preys <- guiGetSafe("preys")
    preds <- guiGetSafe("preds")
    preys.SI <- guiGetSafe("preys.SI")
    preds.SI <- guiGetSafe("preds.SI")
    preys.ix <- guiGetSafe("preys.ix")
    fc_mean <- guiGetSafe("fc_mean")
    fc_sd <- guiGetSafe("fc_sd")
    mean_c <- guiGetSafe("mean_c")
    sd_c <- guiGetSafe("sd_c")
    mean_cs <- guiGetSafe("mean_cs")
    sd_cs <- guiGetSafe("sd_cs")
    props <- guiGetSafe("props")
    Covs <- guiGetSafe("Covariates")
    beta <- guiGetSafe("beta")
    Grps <- guiGetSafe("Grps")
    switchin <- guiGetSafe("gswitch")
    
    if(all(!is.na(preys))){
      write.table(cbind(preys.ix,preys),file=paste(filename,'_FA_preys.csv',sep=''),sep=',',quote=F,col.names=T,row.names=F)
      write.table(preds,file=paste(filename,'_FA_preds.csv',sep=''),sep=',',quote=F,col.names=T,row.names=T)
      write.table(cbind(fc_mean,fc_sd),file=paste(filename,'_fat_cont.csv',sep=''),sep=',',quote=F,col.names=F,row.names=T)
      write.table(mean_c,file=paste(filename,'_FA_cc_means.csv',sep=''),sep=',',quote=F,col.names=T,row.names=T)
      write.table(sd_c,file=paste(filename,'_FA_cc_var.csv',sep=''),sep=',',quote=F,col.names=T,row.names=T)
    }
    if(all(!is.na(preys.SI))){
      write.table(cbind(preys.ix,preys.SI),file=paste(filename,'_SI_preys.csv',sep=''),sep=',',quote=F,col.names=T,row.names=F)
      write.table(preds.SI,file=paste(filename,'_SI_preds.csv',sep=''),sep=',',quote=F,col.names=T,row.names=T)
      write.table(mean_cs,file=paste(filename,'_SI_fc_means.csv',sep=''),sep=',',quote=F,col.names=T,row.names=T)
      write.table(sd_cs,file=paste(filename,'_SI_fc_var.csv',sep=''),sep=',',quote=F,col.names=T,row.names=T)
    }
    if(switchin == 'Covas' | switchin == 'combined'){
      write.table(Covs,file=paste(filename,'_Covariates.csv',sep=''),sep=',',quote=F,col.names=T,row.names=F)
      write.table(beta,file=paste(filename,'_Cov_n_Grp_effects.csv',sep=''),sep=',',quote=F,col.names=T,row.names=F)}
     if(switchin == 'groups' | switchin == 'combined'){
      write.table(Grps,file=paste(filename,'_Groups.csv',sep=''),sep=',',quote=F,col.names=T,row.names=F)
      write.table(beta,file=paste(filename,'_Cov_n_Grp_effects.csv',sep=''),sep=',',quote=F,col.names=T,row.names=F)
    }
          
    write.table(props,file=paste(filename,'_props.csv',sep=''),sep=',',quote=F,col.names=F,row.names=T)
    
  }
  
  simplott <- function(){
    preya <- {}
    preda <- {}

    n.preds <- guiGetSafe("n.preds")
    n.preys <- guiGetSafe("n.preys")
    preys   <- guiGetSafe("preys")
    preds   <- guiGetSafe("preds")
    preys.SI   <- guiGetSafe("preys.SI")
    preds.SI   <- guiGetSafe("preds.SI")
    mean_cs   <- guiGetSafe("mean_cs")
    preys.ix <- guiGetSafe("preys.ix")
    
    plott <-F
    if(all(!is.na(preys))){preya=cbind(preya,clr(preys));preda=cbind(preda,clr(preds));plott <-T}
    if(all(!is.na(preys.SI))){preya=cbind(preya,preys.SI);preda=cbind(preda,t(t(preds.SI)-colMeans(mean_cs)));plott <-T}
    
    if(plott==T){
           
      dista <- dist(rbind(preya,preda))
      mds <- metaMDS(dista)
      
      if (!is.function(options()$device)){
        if (names(dev.cur())=="RStudioGD"){
          # try to open a new platform-appropriate plot window
          if (.Platform$OS.type=='windows'){
            windows()
          } else if(length(grep(R.version$platform,pattern='apple'))>0)  # is it mac?
          { 
            quartz(width=5,height=5)
          } else {  # must be unix
            x11()
          }
          externalDevice<-TRUE
        }
      }
      pl <- plot(mds,type='n')
      points(pl,'sites',pch=cbind(as.numeric(as.factor(preys.ix)),rep(16,n.preds)),col=cbind(1+as.numeric(as.factor(preys.ix)),rep(1,n.preds)))
      legend('bottomright',c('Predators',unique(preys.ix)),xpd=T,pch=c(16,1:n.preys),col=c(1,2:(n.preys+1)))
      
    }
}
  
  #require(fgui)
  
  output <- gui(simulator,callback=guiExec,
                argCommand=list(simplot=simplott,
                                simProps = guiNestedF(simProps,'Props',exec='Simulate simple proportions',
                                                       argText=list(eveness = 'Eveness of simulated diet proportions'),
                                                       argSlider=list(eveness=c(0.05,1,0.01)),helps=''),
                                simGroups = guiNestedF(simGroups,'Groups',exec='Simulate grouped diet proportions',
                                                       argText=list(n.groups = 'Number of groups'),
                                                       argSlider=list(group.effect.size=c(0,5,0.1)),helps=''),
                                simCovs = guiNestedF(simCovs,'Covars',exec='Simulate diet proportions from covariates',
                                                     argText=list(n.covs = 'Number of covariates'),
                                                     argSlider=list(covariate.effect.size=c(0,5,0.1)),helps=''),
                                simwrite = guiNestedF(simwrite,'simwrite',exec='Write simulation to files',helps=''),
                                FAsim.data = guiNestedF(sim.FA,'FAsim.data',
                                                        exec='Simulate Fatty Acid data',
                                                        cancelButton=F,closeOnExec = TRUE,output = NULL,
                                                        argSlider=list(
                                                          sep=c(0.1,10,0.1),
                                                          n.fats=c(3,40,1),
                                                          cvar=c(0.05,0.5,0.01)),
                                                        argText=list(sep = 'Prey separation in FA space',
                                                                     n.fats= 'Number of fatty acids',
                                                                     cvar='Variability of conversion coefficients'
                                                        ),helps=''
                                ),
                                SIsim.data = guiNestedF(sim.SI,'SIsim.data',
                                                        exec='Simulate Stable Isotope data',
                                                        cancelButton=F,closeOnExec = TRUE,output = NULL,
                                                        argSlider=list(
                                                          sep=c(0.1,5,0.1),
                                                          isos=c(2,10,1)
                                                        ),
                                                        argText=list(sep = 'Prey separation in SI space',
                                                                     isos= 'Number of Stable Isotopes'                                                                  
                                                        )
                                                        ,helps='')
                ),                
                argText=list(FAsim.data = 'Simulate Fatty Acid data',
                             SIsim.data = 'Simulate Stable Isotope data',                      
                             simwrite = 'Write simulation to files',
                             simplot = 'Plot current simulation on NMDS axes',
                             simCovs = 'Simulate covariates',
                             simGroups = 'Simulate grouped predator data',
                             simProps = 'Simulate simple proportions',
                             n.preys = 'Number of prey species',
                             n.preds = 'Number of individual predators',
                             nsamples = 'Number of samples per prey species'),
                argSlider=list(n.preys=c(3,10,1),
                               n.preds=c(2,40,1),
                               nsamples=c(5,100,5)
                              ),exec='',closeOnExec = F,argGridOrder=c(1,2,3,4,4,4,5,5,6,6),output=NULL,helps=list(simwrite=simulation)
  ) 
  #detach("package:fgui", unload=TRUE)
  #detach("package:tcltk", unload=TRUE)
}    