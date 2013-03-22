fastin <- function(SI.data=NULL,FA.data=NULL,Groups=NULL,Covariates=NULL,eveness=0.1,nIter=10000,nBurnin=5000,Data.Type='Fatty.Acid.Profiles',Analysis.Type='Population.proportions'){
  
  # combine sources function
  source.combine <- function(k,preys.ix){
    # combination choice
    cat('please select from source combination menu','\n')
    combine <- menu(title='combine sources into groups',choices = c('yes','no'),graphics=T)
    
    #in case the menyu is clicked away without choice
    while (combine==0) {cat('please select from source combination menu')
                        combine <- menu(title='combine sources into groups',choices = c('yes','no'),graphics=T)
    }
    
    if (combine==1) {
      
      selecta <- select.list(preys.names,multiple = T,graphics=T)
      #       cat(selecta,"\n")
      #       cat(k,"\n")
      #       unique(preys.ix)[selecta]
      preys.ix[preys.ix %in% selecta] <- paste('grouped prey',k)
      
      preys.ix <- source.combine(k+1,preys.ix)
      
    }
    
    return(preys.ix)
    
  }
  
  ## FA data prep ---------------
  if (Data.Type=='Fatty.Acid.Profiles' | Data.Type == 'Combined.Analysis'){
    if(is.list(FA.data)){
      
      # import predator and prey FA profiles
      predators = read.csv(FA.data$predators,header=F,row.names=1)
      preys = read.csv(FA.data$preys,header=F)
      
      preys.ix <- as.character(preys[,1])
      
      preys.names <- as.character(unique(preys.ix))
      
      
      preys = preys[,-1]
      
      n.fats = ncol(predators)
      m.fats = n.fats-1
      
      #number of preys species/groups
      n.preys <- length(unique(preys.ix))
      
      # treat conversion coeffs
      if (nchar(FA.data$Conv.Coeffs.mean)>0 & nchar(FA.data$Conv.Coeffs.var)>0)
      {     
        mean_c = read.csv(FA.data$Conv.Coeffs.mean,header=F,row.names=1)
        tau_c  = 1/read.csv(FA.data$Conv.Coeffs.var,header=F,row.names=1)^2
      }
      else if (nchar(FA.data$Conv.Coeffs.mean)==0 & nchar(FA.data$Conv.Coeffs.var)==0)
      {
        mean_c = matrix(FA.data$CC.mean,n.fats,n.preys)
        tau_c =matrix(1/FA.data$CC.var,n.fats,n.preys)
      } else
      {
        stop('Known conversion coefficients, or a mean AND variance for conversion coefficients need to be supplied')
      }
      
      # deal with fat content
      if(nchar(FA.data$fat.cont)==0) 
      {
        fc.mean <- FA.data$FC.mean; fc.var <- FA.data$FC.var
      } else
      {
        
        if (dim(fat.cont)[2]>1){
          fat.cont <- read.csv(FA.data$fat.cont,header=F,row.names=1) 
          fc.mean = fat.cont[,1];fc.sd = fat.cont[,2]
        } else { fat.cont <- read.csv(FA.data$fat.cont,header=F) }
      }
      
      # make sure everything sums to 1
      
      predators <- t(apply(predators,1,function(x){x/(sum(x))}))
      preys <- t(apply(preys,1,function(x){x/(sum(x))}))
      
      ## first, calculate distances  ---------------
      
      require('robCompositions')
      
      dists <- matrix(,nrow(preys),nrow(preys))
      for (i in 1:nrow(preys)){
        for (j in i:nrow(preys)){
          dists[j,i] <- aDist(preys[i,],preys[j,])
        }}
      detach("package:robCompositions", unload=TRUE)
      dista <- as.dist(dists)
      
      x11()
      PR.RDA <- capscale(dista~as.factor(preys.ix),comm=preys)
      #plot(PR.RDA,t='n',xlim=c(-0.5,0.5),ylim=c(-1,1))
      
      plot(preys%*%PR.RDA$CCA$v[,1:2],pch=as.numeric(as.factor(preys.ix)),col=as.numeric(as.factor(preys.ix)))
      points(predators%*%PR.RDA$CCA$v[,1:2],pch=16)
      #cat('please select lower right and upper left corner for legend','\n','(can be outside of plot region)')
      
      legend(locator(2),legend=(unique(preys.ix)),xpd=T,pch=1:n.preys,col=2:(n.preys+1))
      
      #recursive call to combine sources
      prey.ix <- source.combine(1,preys.ix)
      
      SC <- T # any(prey.ix != preys.ix)         
      
      n.preys <- length(unique(prey.ix))
      n.preds <- dim(predators)[1]
      
      ## NOW DO variable selection -----
      x11()  
      plot(cumsum(sort(clo(rowSums(t(t(cbind(PR.RDA$CCA$v,PR.RDA$CA$v))*c(PR.RDA$CCA$eig,PR.RDA$CA$eig))^2)),decreasing =T)),ylab='cumulative proportion')
      
      print(cumsum(sort(clo(rowSums(t(t(cbind(PR.RDA$CCA$v,PR.RDA$CA$v))*c(PR.RDA$CCA$eig,PR.RDA$CA$eig))^2)),decreasing =T)))
      
      # new CAP for vs
      answer <- menu(c('yes','no'),'Would you like to use a subset of Fatty Acids?',graphics=T)
      
      while(answer==0) answer <- menu(c('yes','no'),'Would you like to use a subset of Fatty Acids?',graphics=T)
      
      if (answer==1){
        #   PR.RDA <- capscale(dista~as.factor(prey.ix),comm=preys)
        #   plot(preys%*%PR.RDA$CCA$v[,1:2],pch=1:n.preys,col=2:(n.preys+1))
        #   points(predators%*%PR.RDA$CCA$v[,1:2],pch=16)
        #cat('please select lower right and upper left corner for legend','\n','(can be outside of plot region)')
        #legend(locator(2),(unique(preys.ix)),xpd=T,pch=1:n.preys,col=2:(n.preys+1))
        
        sv = sort(clo(rowSums(t(t(cbind(PR.RDA$CCA$v,PR.RDA$CA$v))*c(PR.RDA$CCA$eig,PR.RDA$CA$eig))^2)),decreasing =T,index.return=T)
        par(ask=T)
        nv <- readline(prompt = "please enter number of variables for analysis \n")
        six <- sv$ix[1:as.numeric(nv)]
        
        n.fats <- length(six)
        m.fats=n.fats-1
      } else {six = 1:n.fats}
      
      # get prey means - loop
      mprey <- matrix(,n.preys,n.fats)
      fcc.mean <- rep(NA,n.preys)
      fcc.sd <- rep(NA,n.preys)
      
      for (i in 1:n.preys){
        
        if (is.null(dim(fat.cont))){
          mprey[i,] <- apply(preys[prey.ix==unique(prey.ix)[i],six]*mean_c[match(preys.ix[prey.ix==unique(prey.ix)[i]],rownames(mean_c)),six],
                             2,function(x){weighted.mean(x,w=fat.cont[prey.ix==unique(prey.ix)[i]])})
          
          fcc.mean[i] <- mean(fat.cont[prey.ix==unique(prey.ix)[i]])
          fcc.sd[i] <- var(fat.cont[prey.ix==unique(prey.ix)[i]])
        } else # combine means and variance
        {
          mprey[i,] <- apply(preys[prey.ix==unique(prey.ix)[i],six]*mean_c[match(preys.ix[prey.ix==unique(prey.ix)[i]],rownames(mean_c)),six],
                             2,function(x){weighted.mean(x,w=fc.mean[match(preys.ix[prey.ix==unique(prey.ix)[i]],rownames(mean_c))])})
          
          fcc.mean[i] <- mean(fc.mean[rownames(mean_c) %in% preys.ix[prey.ix==unique(prey.ix)[i]]])
          fcc.sd[i] <- mean(fc.sd[rownames(mean_c) %in% preys.ix[prey.ix==unique(prey.ix)[i]]])
        }
      }
      
      mean_c <- mean_c[1:n.preys,six]*0+1         
      
      preym <- unclass(alr(mprey))
      preds <- unclass(alr(predators[,six]))    
      
      # now prepare data for analysis
      
      R <- array(,c(m.fats,m.fats,n.preys))
      ni<-rep(NA,n.preys)
      for (i in 1:n.preys){
        ni[i] <- max(n.fats+1,sum(prey.ix==unique(prey.ix)[i])-1)
        R[,,i]=cov(alr(preys[prey.ix==unique(prey.ix)[i],six]))*ni[i]
      }
      
      ## first some data and inits ----
      
      # set uninformative prior SS matrix for wishart prior alr transformed predator data
      
      Rnot =diag(R.diag,m.fats)
      
      FA==T
      datas.FA <- list(fc_mean=fcc.mean,fc_tau=1/fcc.sd^2,n.fats=n.fats,m.fats=m.fats,R=R,Rnot=Rnot,preym=preym,preds=preds,ni=ni,mean_c=mean_c,tau_c=tau_c)
      
      
    } else {
      stop("Must supply Fatty Acid Profile data when choosing option 'Fatty Acid Profiles' or 'Combined analysis' as data type")
    }
  }
  # switch between cases for analyses
  
  if (Data.Type=='Stable.Isotopes' | Data.Type == 'Combined.Analysis'){
    ## SI data prep -------------
  if(is.list(SI.data)){
    
    ## first check for potential conflicts
    
    stopifnot(nchar(SI.data$predators.SI)>0 & nchar(SI.data$preys.SI)>0)
    
    # import predator and prey data - note that the first column is names, or an index
    predators.SI = read.csv(SI.data$predators.SI,header=F,row.names=1)
    preys.SI = read.csv(SI.data$preys.SI,header=F)
    
    preys.ix.SI  <- as.character(preys.SI[,1])
    
    preys.names.SI  <- as.character(unique(preys.ix.SI ))
    
    
    preys.SI  = preys.SI [,-1]
    
    # set number of isotopes
    isos=ncol(predators.SI)
    
    # improt prior predator variance
    Rnot_SI = diag(SI.data$R.diag.SI,isos)
    
    # deal with fractionation coeffs
    if ((nchar(SI.data$Frac.Coeffs.mean)>0 & nchar(SI.data$Frac.Coeffs.var)==0) | (nchar(SI.data$Frac.Coeffs.mean)==0 & nchar(SI.data$Frac.Coeffs.var)>0))
    {
      stop('The mean AND variances of FCs for each isotope need to be supplied in the form c(FC1,FC2)')
    } else if (nchar(SI.data$Frac.Coeffs.mean)>0 & nchar(SI.data$Frac.Coeffs.var)>0)
    {     
      mean_cs = read.csv(SI.data$Frac.Coeffs.mean,header=F,row.names=1)
      sd_cs  = 1/read.csv(SI.data$Frac.Coeffs.var,header=F,row.names=1)^2
    } else if (nchar(SI.data$Frac.Coeffs.mean)==0 & nchar(SI.data$Frac.Coeffs.var)==0)
    {
      mean_cs = matrix(SI.data$FC.mean,isos,n.preys.SI)
      tau_cs =matrix(1/SI.data$FC.var,isos,n.preys.SI)
    }
    
    # if we've allready combine the preys based on FAs, then jsut use the index and combine
    # else combine preys here if desired
    
    if (SC==F) # query for prey combination
    {
      if (isos<=2){
        
        plot(preys.SI[,1:ncol(preys.SI)],pch=as.numeric(as.factor(preys.ix.SI)),col=as.numeric(as.factor(preys.ix.SI))+1)
        points(predators.SI[,1:2],pch=16)
        tkmessageBox( message="please use the cursor to select lower right and upper left corner for legend", title="Plot Legend" )
        legend(locator(2),legend=preys.names.SI,xpd=T,pch=1:n.preys,col=2:(n.preys+1))
        #cat('please select lower right and upper left corner for legend','\n','(can be outside of plot region)')
        #legend(locator(2),(unique(preys.ix.SI)),xpd=T,pch=1:n.preys.SI,col=2:(n.preys.SI+1))
        
      } else {
        
        PR.RDA <- capscale(preys.SI~as.factor(preys.ix.SI))
        plot(data.matrix(preys.SI)%*%PR.RDA$CCA$v[,1:2],pch=as.numeric(as.factor(preys.ix.SI)),col=as.numeric(as.factor(preys.ix.SI))+1)
        points(data.matrix(predators.SI)%*%PR.RDA$CCA$v[,1:2],pch=16)
        tkmessageBox(message="please use the cursor to select lower right and upper left corner for legend", title="Plot Legend" )
        legend(locator(2),legend=preys.names.SI,xpd=T,pch=1:n.preys,col=2:(n.preys+1))
        #cat('please select lower right and upper left corner for legend','\n','(can be outside of plot region)')
        #legend(locator(2),(unique(preys.ix.SI)),xpd=T,pch=1:n.preys.SI,col=2:(n.preys.SI+1))
        
      }
      
      #recursive call to combine sources
      prey.ix <- source.combine(1,preys.ix.SI)
      
      n.preys <- length(unique(prey.ix))
      n.preds <- dim(predators.SI)[1]
    }
    
    
    # combine preys
    preym.SI <- matrix(,n.preys,isos) 
    # in either case, combine preys....
    for (i in 1:n.preys){
      
      preym.SI[i,] <- apply(preys.SI[prey.ix==unique(prey.ix)[i],]+mean_cs[match(preys.ix.SI[prey.ix==unique(prey.ix)[i]],rownames(mean_cs)),],2,mean)
      
    }
    
    #set fractionation to 0 since it's allready applied
    mean_cs <- mean_cs[1:n.preys,]*0         
    
    # now prepare data for analysis
    
    R.SI <- array(,c(isos,isos,n.preys))
    ni.SI<-rep(NA,n.preys)
    for (i in 1:n.preys){
      ni.SI[i] <- max(isos+1,sum(prey.ix==unique(prey.ix)[i])-1)
      R.SI[,,i]=cov(preys.SI[prey.ix==unique(prey.ix)[i],])*ni.SI[i]
    }
    
    ## first some data and inits ----
    
    # set uninformative prior SS matrix for wishart prior alr transformed predator data
    
    SI==T
    datas.SI <- list(isos=isos,R.SI=R.SI,Rnot.SI=Rnot.SI,preym.SI=preym.SI,preds.SI=preds.SI,ni.SI=ni.SI,mean_cs=mean_cs,tau_cs=tau_cs)
    
  } else {
    stop("Must supply Stable Isotope data when choosing option 'Stable Isotopes' or 'Combined analysis' as data type")
  }
}

datas <- list(n.preys = n.preys,n.preds=n.preds,datas.FA=NULL,datas.SI=NULL,eveness=eveness)

if (Data.Type == 'Combined.Analysis')
{
  datas$datas.FA <- datas.FA
  datas$datas.SI <- datas.SI
  class(datas) <- 'combined'
} else if (Data.Type=='Fatty.Acid.Profiles')
{
  datas$datas.FA <- datas.FA
  class(datas) <- 'FA'
} else if (Data.Type=='Stable.Isotopes')
{
  datas$datas.SI <- datas.SI
  class(datas) <- 'SI'
}


if (nchar(Covariates)>0 & nchar(Groups)==0)
{
  Analysis.Type = 'Analysis.with.Covariates'    
  Covs <- read.csv(Covariates,header=T)
  Covs <- cbind(rep(1,nrow(Covs)),Covs)
  n.covs <- ncol(Covs)
} else if (nchar(Covariates)==0 & nchar(Groups)>0) 
  {
  Analysis.Type = 'Analysis.with.Covariates'
  Grps <- read.csv(Groups,header=F)
  Grp.names <- Grps[,1]    
  Grps <- Grps[,2]  
  
  for (i in 1:ncol(Grps)){
    vg <- as.vector(Grps[,i])
    Grps[,i] <- as.factor(vg)
  }
  
  Covs <- model.matrix(attr(model.frame(1:nrow(Grps)~.,data=Grps),'terms'),data=Grps)[,]
  names(Covs) <- Grp.names
  
  
} else if (nchar(Covariates)>0 & nchar(Groups)>0) 
  {
  Analysis.Type = 'Analysis.with.Covariates'
  Covs <- read.csv(Covariates,header=T)
  Covnames <- names(Covs)
  Grps <- read.csv(Groups,header=F)
  Grp.names <- Grps[,1]    
  Grps <- Grps[,2] 
  
  for (i in 1:ncol(Grps)){
    vg <- as.vector(Grps[,i])
    Grps[,i] <- as.factor(vg)
  }
  
  Covs <- cbind(model.matrix(attr(model.frame(1:nrow(Grps)~.,data=Grps),'terms'),data=Grps)[,],Covs)
  names(Covs) <- c(Grp.names,Covnames)
  
}

# have three types here: FA, SI and combined, then methods dispatch based on type of arg
outputs <- switch(Analysis.Type,
                  Population.proportions = Poppropanalysis(datas,nIter,nBurnin),
                  Individual.proportions = PopandIndprops(datas,nIter,nBurnin),
                  Analysis.with.Covariates = AnalysiswithCov(datas,Covs,nIter,nBurnin)
)

return(output)

# plot output in here? - methods dispatch here too
#   * ggdensity of individual, and all psots
#   * shaded density plot
#   * boxplot
#   * trace of MCMC

# sumamry of output - methdos dispatch

}