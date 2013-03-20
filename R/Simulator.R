# SIMULATOR
simulation <- function(){
  
  ######### sim functions
  
  simulator <- function(n.preys=3,n.preds=10,nsamples=30,eveness=0.5,FAsim.data=NULL,SIsim.data=NULL,simwrite=NULL){
    
    preys.ix={};
    for (i in 1:n.preys){
      preys.ix <- c(preys.ix,rep(paste('Prey',i,sep='_'),nsamples))
      
    }
    
    preys.ix.num = as.numeric(as.factor(preys.ix))
       
    v1 = 0.5; v2 = 0.7
    Q = rlnorm(n.preys,c(rep(0,n.preys-1),eveness),v2)
    #rlnorm(n.preys,c(rep(0,n.preys-1),eveness),v2)# individual proportions
    props=rdirichlet(n.preds,Q)
    
    guiSet("props",props)
    guiSet("n.preys",n.preys)
    guiSet("n.preds",n.preds)
    guiSet("nsamples",nsamples)
    guiSet("preys.ix",preys.ix)
    
    preya <- {}
    preda <- {}
    
    #n.preys <- guiGetSafe("n.preys")
    preys   <- guiGetSafe("preys")
    preds   <- guiGetSafe("preds")
    preys.SI   <- guiGetSafe("preys.SI")
    preds.SI   <- guiGetSafe("preds.SI")
    #preys.ix <- guiGetSafe("preys.ix")
    plott <-F
    if(all(!is.na(preys))){preya=cbind(preya,alr(preys));preda=cbind(preda,alr(preds));plott <-T}
    if(all(!is.na(preys.SI))){preya=cbind(preya,preys.SI);preda=cbind(preda,preds.SI);plott <-T}
    
    if(plott==T){
     x11() 
    PR.RDA <- capscale(preya~as.factor(preys.ix))
    plot(preya%*%PR.RDA$CCA$v[,1:2],pch=as.numeric(as.factor(preys.ix)),col=as.numeric(as.factor(preys.ix))+1)
    points(preda%*%PR.RDA$CCA$v[,1:2],pch=16)
    tkmessageBox( message="please use the cursor to select lower right and upper left corner for legend", title="Plot Legend" )
     legend(locator(2),(unique(preys.ix)),xpd=T,pch=1:n.preys,col=2:(n.preys+1))
    
    #list(preys = preya,preys.ix = preys.ix,preds = preds,props = props,mean_c = mean_css,sd_c = sd_css, fc_mean = fc_mean,fc_sd = fc_sd)
    }
  }
  
  sim.FA <- function(sep=1,n.fats=10,cvar=0.5){
    
    props <- guiGetSafe("props")
    n.preys <- guiGetSafe("n.preys")
    n.preds <- guiGetSafe("n.preds")
    nsamples <- guiGetSafe("nsamples")
    preys.ix <- guiGetSafe("preys.ix")
    
    As = matrix(rlnorm(n.preys*n.fats,rnorm(n.preys*n.fats,sep),sep),n.preys,n.fats,byrow=T)
    
    preys <- array(,c(nsamples[1],n.fats,n.preys))
    for (i in 1:n.preys){
      preys[,,i] <- rdirichlet(nsamples,As[i,])
    }
    
    preya={}
    for (i in 1:n.preys){
      preya<- rbind(preya,preys[,,i])
    }
    
    preys.ix.num = as.numeric(as.factor(preys.ix))
    
    
    #randomly draw fat content
    fc_mean <- as.matrix(sample(40,n.preys))/10+1
    rownames(fc_mean) <- unique(preys.ix)
    # set variance of fat content to 0.4 - this is loosely based on literature
    fc_sd <- matrix(0.5,n.preys,1)
    
    # draw mean conversion coeffs - cvar is the variance of mean conversion coeffs around 1; 0.2 seems reasonable  (c.f., Rosen & Tollit histogram)
    mean_css=matrix(rlnorm(n.fats*n.preys,0,cvar),n.preys,n.fats) # rtnorm is a truncated normal, with lower truncation set to 0
    rownames(mean_css) <- unique(preys.ix)
    
    sd_css=matrix(0.05,n.preys,n.fats)
    rownames(sd_css) <- unique(preys.ix)
    
    mprey <-  clo(t(apply(preys,3,colMeans)))
    preds <-  clo(props%*%(as.vector(fc_mean)*mprey*mean_css))
    
    
    guiSet("preds",preds)
    guiSet("preys",preya)
    
    
    guiSet("mean_c",mean_css)
    guiSet("sd_c",sd_css)
    
    guiSet("fc_mean",fc_mean)
    guiSet("fc_sd",fc_sd)
  }
  
  sim.SI <- function(sep=1,isos=2){
    
    props <- guiGetSafe("props")
    n.preys <- guiGetSafe("n.preys")
    n.preds <- guiGetSafe("n.preds")
    n.samples <- guiGetSafe("n.samples")
    preys.ix <- guiGetSafe("preys.ix")
    
    
    prey.m <- mvrnorm(isos*n.preys,c(-16,18,rnorm(isos-2,0,sep)),diag(sep,isos))
    prey.cov <- diag(1/sep,isos)
    
    preys.SI <- array(NA,c(nsamples,isos,n.preys))
    for (i in 1:n.preys){
      preys.SI[,,i] <-mvrnorm(nsamples,prey.m[i,],prey.cov)
    }
    
    preya.SI={}
    for (i in 1:n.preys){
      preya.SI<- rbind(preya.SI,preys.SI[,,i])
    }
    
    mean_cs = matrix(rnorm(isos*n.preys,c(3,1,rnorm(isos-2))),n.preys,isos,byrow=T)
    rownames(mean_cs) <- unique(preys.ix)
    
    sd_cs =matrix(c(0.2,0.2),n.preys,isos) # 1.61
    rownames(sd_cs) <- unique(preys.ix)
    
    # the mixture in matrix algebra (just a matrix product)
    preym.SI <-  t(apply(preys.SI,3,colMeans))
    preds.SI <-  props %*% (preym.SI+mean_cs)
    
    
    guiSet("preds.SI",preds.SI)
    guiSet("preys.SI",preya.SI)
    
    guiSet("mean_cs",mean_cs)
    guiSet("sd_cs",sd_cs)
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
    
    if(all(!is.na(preys))){
      write.table(cbind(preys.ix,preys),file=paste(filename,'_preys.csv',sep=''),sep=',',quote=F,col.names=F,row.names=F)
      write.table(preds,file=paste(filename,'_preds.csv',sep=''),sep=',',quote=F,col.names=F,row.names=T)
      write.table(cbind(fc_mean,fc_sd),file=paste(filename,'_fat_cont.csv',sep=''),sep=',',quote=F,col.names=F,row.names=T)
      write.table(mean_c,file=paste(filename,'_FA_cc_means.csv',sep=''),sep=',',quote=F,col.names=F,row.names=T)
      write.table(sd_c,file=paste(filename,'_FA_cc_sd.csv',sep=''),sep=',',quote=F,col.names=F,row.names=T)
    }
    if(all(!is.na(preys.SI))){
      write.table(cbind(preys.ix,preys.SI),file=paste(filename,'_SI_preys.csv',sep=''),sep=',',quote=F,col.names=F,row.names=F)
      write.table(preds.SI,file=paste(filename,'_SI_preds.csv',sep=''),sep=',',quote=F,col.names=F,row.names=T)
      write.table(mean_cs,file=paste(filename,'_SI_fc_means.csv',sep=''),sep=',',quote=F,col.names=F,row.names=T)
      write.table(sd_cs,file=paste(filename,'_SI_fc_sd.csv',sep=''),sep=',',quote=F,col.names=F,row.names=T)
    }
    write.table(props,file=paste(filename,'_props.csv',sep=''),sep=',',quote=F,col.names=F,row.names=T)
    
  }
 
  
  output <- gui(simulator,
                argCommand=list(
                  simwrite = guiNestedF(simwrite,'simwrite',exec='Write simulation to files'),
                                FAsim.data = guiNestedF(sim.FA,'FA.data',
                                                     exec='Simulate Fatty Acid data',
                                                     cancelButton=F,closeOnExec = TRUE,output = NULL,
                                                     argSlider=list(
                                                       sep=c(0,2,0.1),
                                                       n.fats=c(2,40,1),
                                                       cvar=c(0.05,0.5,0.01)),
                                                     argText=list(sep = 'Prey separation in FA space',
                                                                  n.fats= 'Number of fatty acids',
                                                                  cvar='variability of conversion coefficients'
                                                     )
                                ),
                                SIsim.data = guiNestedF(sim.SI,'SI.data',
                                                     exec='Simulate Stable Isotope data',
                                                     cancelButton=F,closeOnExec = TRUE,output = NULL,
                                                     argSlider=list(
                                                       sep=c(0.1,5,0.1),
                                                       isos=c(1,10,1)
                                                     ),
                                                     argText=list(sep = 'Prey separation in FA space',
                                                                  isos= 'Number of Stable Isotopes'                                                                  
                                                     )
                                )
                ),                
                argText=list(FAsim.data = 'Simulate Fatty Acid data',
                             SIsim.data = 'Simulate Stable Isotope data',                      
                             simwrite = 'Write simulation to files',
                             n.preys = 'Number of prey species',
                             n.preds = 'Number of individual predators',
                             nsamples = 'number of samples per prey species',
                             eveness = 'Eveness of simulated diet proportions'),
                  argSlider=list(n.preys=c(2,10,1),
                  n.preds=c(2,40,1),
                  nsamples=c(5,100,5),
                  eveness=c(0.05,1,0.01)),exec='Update simulation',closeOnExec = F,argGridOrder=c(1,2,3,4,5,5,5),output=NULL
  )
  
}    