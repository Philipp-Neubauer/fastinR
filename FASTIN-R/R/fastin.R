.fastin <- function(SI.data=NULL,FA.data=NULL,groupings,Load.Data=NULL,Save.Data=NULL,add.covs,MCMC,Save.Outputs=NULL,Display.Summaries=NULL,Plot.Outputs=NULL,disp.diags=NULL){}

addSI <- function(predators.SI=NULL,preys.SI=NULL,Frac.Coeffs.mean=NULL,Frac.Coeffs.var=NULL,FC.mean=1,FC.var=1,R.diag.SI=0.01){
        
    # combine sources function
        source.combine <- function(k,preys.ix,preys.names){

            n.preys <- length(unique(preys.ix))
            preys.names <- as.character(unique(preys.ix))
            
            if (isos<=2){
                x11()
                plot(preys.SI[,1:ncol(preys.SI)],pch=as.numeric(as.factor(preys.ix)),col=as.numeric(as.factor(preys.ix))+1)
                points(predators.SI[,1:2],pch=16)
                legend("topright",legend=preys.names,xpd=T,pch=1:n.preys,col=2:(n.preys+1))
                                  
        
            } else {
                x11()
                dista <- dist(rbind(preys.SI,predators.SI))
                mds <- metaMDS(dista)
                pl <- plot(mds,type='n')
                points(pl,'sites',pch=cbind(as.numeric(as.factor(preys.ix)),rep(16,n.preds)),col=cbind(1+as.numeric(as.factor(preys.ix)),rep(1,n.preds)))
                legend('bottomright',c('Predators',unique(preys.ix)),xpd=T,pch=c(16,1:n.preys),col=c(1,2:(n.preys+1)))
          
      }
      # combination choice
      cat('please select from source combination menu','\n')
      combine <- menu(title='combine sources into groups',choices = c('yes','no'),graphics=T)
      
      #in case the menyu is clicked away without choice
      while (combine==0) {cat('please select from source combination menu')
                          combine <- menu(title='combine sources into groups',choices = c('yes','no'),graphics=T)
      }
      
        if (combine==1) {
        dev.off()
        selecta <- select.list(preys.names,multiple = T,graphics=T)
        #       cat(selecta,"\n")
        #       cat(k,"\n")
        #       unique(preys.ix)[selecta]
        preys.ix[preys.ix %in% selecta] <- paste('grouped prey',k)

        preys.ix <- source.combine(k+1,preys.ix)
        
      }
      
      return(preys.ix)
            return(datas)
    }
    
    datas <- guiGetSafe('datas')
    
    ## first check for potential conflicts
    
    stopifnot(nchar(predators.SI)>0 & nchar(preys.SI)>0)
    
    # import predator and prey data - note that the first column is names, or an index
    predators.SI = read.csv(predators.SI,header=T,row.names=1)
    preys.SI = read.csv(preys.SI,header=T)
    
    n.preds <- dim(predators.SI)[1]
    
    preys.ix.SI  <- as.character(preys.SI[,1])
    
    preys.names.SI  <- as.character(unique(preys.ix.SI))
    guiSet('prey.names',preys.names.SI )
    
    preys.SI  = preys.SI [,-1]
    
    # set number of isotopes
    isos=ncol(predators.SI)
    #number of preys species/groups
    n.preys <- length(unique(preys.ix.SI))
    
    # improt prior predator variance
    Rnot_SI = diag(R.diag.SI,isos)
    
    # deal with fractionation coeffs
    if ((nchar(Frac.Coeffs.mean)>0 & nchar(Frac.Coeffs.var)==0) | (nchar(Frac.Coeffs.mean)==0 & nchar(Frac.Coeffs.var)>0))
    {
      stop('The mean AND variances of FCs for each isotope need to be supplied in the form c(FC1,FC2)')
    } else if (nchar(Frac.Coeffs.mean)>0 & nchar(Frac.Coeffs.var)>0)
    {     
      mean_cs = read.csv(Frac.Coeffs.mean,header=F,row.names=1)
      sd_cs  = read.csv(Frac.Coeffs.var,header=F,row.names=1)
    } else if (nchar(Frac.Coeffs.mean)==0 & nchar(Frac.Coeffs.var)==0)
    {
      mean_cs = matrix(FC.mean,isos,n.preys.SI)
      sd_cs =matrix(FC.var,isos,n.preys.SI)
    }
    
    # if we've allready combine the preys based on FAs, then jsut use the index and combine
    if(length(datas)<=1) {
      SC=F
    } else if (datas$SC==T) {
      prey.ix=datas$prey.ix ; n.preys=datas$n.preys;SC=T
    } else {SC=F}
    
    # else combine preys here if desired
    
    if (SC==F) # query for prey combination
    {
          
      #recursive call to combine sources
      prey.ix <- source.combine(1,preys.ix.SI,preys.names.SI)

      n.preys <- length(unique(prey.ix))
      preys.names <- as.character(unique(prey.ix))
        
      guiSet('prey.names',preys.names )
      
    }else{warning('using previously combined sources')}
    SC=T
    if (dev.cur()!=1)
      dev.off()
    # combine preys
    preym.SI <- matrix(,n.preys,isos) 
    var_cs<- matrix(,n.preys,isos) 
    # in either case, combine preys....
    for (i in 1:n.preys){
      
      preym.SI[i,] <- apply(preys.SI[prey.ix==unique(prey.ix)[i],]+mean_cs[match(preys.ix.SI[prey.ix==unique(prey.ix)[i]],rownames(mean_cs)),],2,mean)
      var_cs[i,] <- colMeans(sd_cs[rownames(mean_cs) %in% preys.ix.SI[prey.ix==unique(prey.ix)[i]],])^2
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
    
    datas.SI <- list(isos=isos,R.SI=R.SI,Rnot.SI=Rnot_SI,preym.SI=preym.SI,preds.SI=predators.SI,ni.SI=ni.SI,mean_cs=mean_cs,tau_cs=1/var_cs)
    
    if(length(datas)<=1){
      datas <- list(n.preys = n.preys,n.preds=n.preds,prey.ix=prey.ix,SC=SC,datas.FA=NULL,datas.SI=datas.SI,even=NULL)
    } else {
      datas$datas.SI = datas.SI
      datas$SC = SC
      datas$n.preys = n.preys
      datas$n.preds=n.preds
      datas$prey.ix=prey.ix
    }
    
    guiSet('datas',datas)
    return(datas)
  }

addFA <- function(predators.FA=NULL,preys.FA=NULL,fat.conts = NULL,Conv.Coeffs.mean=NULL,Conv.Coeffs.var=NULL,FC.mean=1,FC.var=1,CC.mean=1,CC.var=1,R.diag=0.01){


    clo <- function(x){
        if (is.null(dim(x))){xc <- x/sum(x)} else {xc <- t(apply(x,1,function(y){y/sum(y)}))}
        return(xc)
    }

    alr <- function(x){
        x<-clo(x)
        if (is.null(dim(x))){log(x[1:(length(x)-1)]/x[length(x)])} else { t(apply(x,1,function(y){log(y[1:(length(y)-1)]/y[length(y)])}))}
    }
    
    kappas <- function(mat,vars){
        ks <- vector(,length(vars))
        for (i in 1:length(vars)){
            ks[i] <- kappa(mat[,vars[1:i]],exact=T)
        }
        return(ks)
    }

    adist <- function(mat){

        dims <- dim(mat)
        dists <- matrix(,dims[1],dims[1])
        for (i in 1:(dims[1]-1)){
            for (j in (i+1):dims[1]){
                dists[j,i] <- robCompositions::aDist(mat[i,],mat[j,])
            }}
        dista <- as.dist(dists)
              return(dista)
    }
        
    # combine sources function
    source.combine <- function(k,preys.ix){
          n.preys <- length(unique(preys.ix))
          preys.names <- as.character(unique(preys.ix))
          x11()
          dista <- dist(rbind(preys,predators))
          mds <- metaMDS(dista)
          pl <- plot(mds,type='n')
          points(pl,'sites',pch=cbind(as.numeric(as.factor(preys.ix)),rep(16,n.preds)),col=cbind(1+as.numeric(as.factor(preys.ix)),rep(1,n.preds)))
          legend('bottomright',c('Predators',unique(preys.ix)),xpd=T,pch=c(16,1:n.preys),col=c(1,2:(n.preys+1)))
         
                                        # combination choice
          cat('please select from source combination menu','\n')
          combine <- menu(title='combine sources into groups?',choices = c('yes','no'),graphics=T)
      
      #in case the menyu is clicked away without choice
          while (combine==0) {cat('please select from source combination menu')
                              combine <- menu(title='combine sources into groups',choices = c('yes','no'),graphics=T)
                          }
      
          if (combine==1) {
              dev.off()
              selecta <- select.list(preys.names,multiple = T,graphics=T)
                                        #       cat(selecta,"\n")
                                        #       cat(k,"\n")
                                        #       unique(preys.ix)[selecta]
              preys.ix[preys.ix %in% selecta] <- paste('grouped prey',k)
        
              preys.ix <- source.combine(k+1,preys.ix)
        
          }
      
          return(preys.ix)
      
      }
    
    
    # import predator and prey FA profiles
    predators = read.csv(predators.FA,header=T,row.names=1)
    preys = read.csv(preys.FA,header=T)
    n.preds <- dim(predators)[1]
    preys.ix <- as.character(preys[,1])
    
    preys.names <- as.character(unique(preys.ix))
    
    guiSet('prey.names',preys.names )
    
    preys = preys[,-1]
    
    n.fats = ncol(predators)
    m.fats = n.fats-1
    
    #number of preys species/groups
    n.preys <- length(unique(preys.ix))
    
    # treat conversion coeffs
    if (nchar(Conv.Coeffs.mean)>0 & nchar(Conv.Coeffs.var)>0)
    {     
      mean_c = read.csv(Conv.Coeffs.mean,header=F,colClasses=c('character',rep('numeric',n.fats)),row.names=1)
      sd_c  = read.csv(Conv.Coeffs.var,header=F,row.names=1)
    } else if (nchar(Conv.Coeffs.mean)==0 & nchar(Conv.Coeffs.var)==0)
    {
      mean_c = matrix(CC.mean,n.preys,n.fats)
      sd_c =matrix(CC.var,n.preys,n.fats)
    } else
    {
      print('Known conversion coefficients, or a mean AND variance for conversion coefficients need to be supplied')
    }
    
    # deal with fat content
      if(nchar(fat.conts)==0) 
          {
              fc.mean <- FC.mean; fc.var <- FC.var
          } else
              {
                  fat.cont <- read.csv(fat.conts,header=F)
                  if (dim(fat.cont)[2]>2){
                      fat.cont <- read.csv(fat.conts,header=F,row.names=1) 
                      fc.mean <- fat.cont[,1];fc.sd <- fat.cont[,2]
                      if (any(fc.mean>1)){
                          fc.mean <- fc.mean/100
                          fc.sd <- fc.sd/100
                      }
                  } else if (any(fat.cont>1)){
                      fat.cont <- fat.cont/100
                  }
              }
  
    
    # make sure everything sums to 1
    
    predators <- t(apply(predators,1,function(x){x/(sum(x))}))
    preys <- t(apply(preys,1,function(x){x/(sum(x))}))
    
    datas <- guiGetSafe('datas')
    if(length(datas)<=1) {
      SC=F
    } else if (datas$SC==T) {
      prey.ix=datas$prey.ix ; n.preys=datas$n.preys; SC=T
    } else {SC=F}
    
    ## first, calculate distances for preys
    dista <- adist(preys)

    # combine sources? SC refers to potential prior combination based on other sources.
    if(SC==F){
            
      #recursive call to combine sources
      prey.ix <- source.combine(1,preys.ix) 
      
      n.preys <- length(unique(prey.ix))
      preys.names <- as.character(unique(prey.ix))
      
      guiSet('prey.names',preys.names )
    } else{warning('using previously combined sources')}
    SC=T
    if (dev.cur()!=1)dev.off()

    PR.RDA <- capscale(dista~as.factor(preys.ix),comm=preys)
    ## NOW DO variable selection -----
    x11()
    par(mar = c(7, 4, 4, 2) + 0.1)
    par(mfcol=c(2,1))
    sv = sort(clo(rowSums(t(t(cbind(PR.RDA$CCA$v,PR.RDA$CA$v))*c(PR.RDA$CCA$eig,PR.RDA$CA$eig))^2)),decreasing =T,index.return=T)
    plot(cumsum(sort(clo(rowSums(t(t(cbind(PR.RDA$CCA$v,PR.RDA$CA$v))*c(PR.RDA$CCA$eig,PR.RDA$CA$eig))^2)),decreasing =T)),axes=F,xlab='',ylab='Cumulative source separation',ylim=c(0,1))
    axis(2)
    axis(1,at=1:n.fats,labels=F)
    text(1:n.fats, par("usr")[3] - 0.2, srt = 45, adj = 1,
     labels = colnames(preys)[sv$ix], xpd = TRUE)

    ks <- kappas(preys,sv$ix)
    plot(ks,axes=F,ylab='Prey matrix condition number',xlab='',ylim=c(0,max(ks)))
    axis(2)
    axis(1,at=1:n.fats,labels=F)
    text(1:n.fats, par("usr")[3] - max(ks)/5 , srt = 45, adj = 1,
     labels = colnames(preys)[sv$ix], xpd = TRUE)

    
    cumsums <- as.data.frame(matrix(,n.fats,1))
    cumsums[,1] <- cumsum(sort(clo(rowSums(t(t(cbind(PR.RDA$CCA$v,PR.RDA$CA$v))*c(PR.RDA$CCA$eig,PR.RDA$CA$eig))^2)),decreasing =T))
    names(cumsums) <- 'Cumulative Proportion'
    rownames(cumsums) <-  colnames(preys)[sv$ix]
    print(cumsums)
    
    # new CAP for vs
    answer <- menu(c('yes','no'),'Would you like to use a subset of Fatty Acids?',graphics=T)
    
    while(answer==0) answer <- menu(c('yes','no'),'Would you like to use a subset of Fatty Acids?',graphics=T)
    
    if (answer==1){
      #   PR.RDA <- capscale(dista~as.factor(prey.ix),comm=preys)
      #   plot(preys%*%PR.RDA$CCA$v[,1:2],pch=1:n.preys,col=2:(n.preys+1))
      #   points(predators%*%PR.RDA$CCA$v[,1:2],pch=16)
      #cat('please select lower right and upper left corner for legend','\n','(can be outside of plot region)')
      #legend(locator(2),(unique(preys.ix)),xpd=T,pch=1:n.preys,col=2:(n.preys+1))
      
      nv <- select.list(title='please choose the fatty acids to use (at least 3)',choices = colnames(preys)[sv$ix],graphics=T,multiple=T)
      
      #nv <- readline(prompt = "please enter number of variables for analysis \n")
      six <- match(nv,colnames(preys))
      n.fats <- length(six)
      m.fats=n.fats-1
    } else {six = 1:n.fats}
    
    if (dev.cur()!=1)
      dev.off()
    
    # get prey means - loop
    mprey <- matrix(,n.preys,n.fats)
    fcc.mean <- rep(NA,n.preys)
    fcc.var <- rep(NA,n.preys)
    var_c <-  matrix(,n.preys,n.fats)
    
    for (i in 1:n.preys){
      
      if (is.null(dim(fat.cont))){
        mprey[i,] <- clo(apply(preys[prey.ix==unique(prey.ix)[i],six]*mean_c[match(preys.ix[prey.ix==unique(prey.ix)[i]],rownames(mean_c)),six],2,function(x){exp(weighted.mean(log(x),w=fat.cont[prey.ix==unique(prey.ix)[i]]))}))
        var_c[i,] <- (sd_c[rownames(mean_c) %in% preys.ix[prey.ix==unique(prey.ix)[i]],six])^2
        fcc.mean[i] <- exp(mean(log(fat.cont[prey.ix==unique(prey.ix)[i]])))
        fcc.var[i] <- exp(sd(log((fat.cont[prey.ix==unique(prey.ix)[i]]))))^2
        
      } else # combine means and variance
      {
        mprey[i,] <- clo(apply(preys[prey.ix==unique(prey.ix)[i],six]*mean_c[match(preys.ix[prey.ix==unique(prey.ix)[i]],rownames(mean_c)),six],2,function(x){exp(weighted.mean(log(x),w=fc.mean[match(preys.ix[prey.ix==unique(prey.ix)[i]],rownames(mean_c))]))}))
        
        fcc.mean[i] <- exp(mean(log(fc.mean[rownames(mean_c) %in% preys.ix[prey.ix==unique(prey.ix)[i]]])))
        fcc.var[i] <- (mean(fc.sd[rownames(mean_c) %in% preys.ix[prey.ix==unique(prey.ix)[i]]]))^2
        var_c[i,] <- (colMeans(sd_c[rownames(mean_c) %in% preys.ix[prey.ix==unique(prey.ix)[i]],six]))^2
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
    
    datas.FA <- list(fc_mean=fcc.mean,fc_tau=1/fcc.var,n.fats=n.fats,m.fats=m.fats,R=R,Rnot=Rnot,preym=preym,preds=preds,ni=ni,mean_c=mean_c,tau_c=1/var_c)
    
    if(length(datas)<=1){
      datas <- list(n.preys = n.preys,n.preds=n.preds,prey.ix=prey.ix,SC=SC,datas.FA=datas.FA,datas.SI=NULL,even=NULL)
    } else {
      datas$datas.FA = datas.FA
      datas$SC = SC
      datas$n.preys = n.preys
      datas$n.preds=n.preds
      datas$prey.ix=prey.ix
    }
    guiSet('datas',datas)
    return(datas)
  }

resetsc <- function(datas=NULL){if (is.null(datas)){datas <- guiGetSafe('datas')};if(length(datas)>1){datas$SC=F;guiSet('datas',datas);return(datas)}}

run_MCMC <- function(datas=NULL,nIter=10000,nBurnin=1000,nChains=1,nThin=10,Data.Type='Fatty.Acid.Profiles',Analysis.Type='Population.proportions',even=0.1,plott=T){
    # have three types here: FA, SI and combined, then methods dispatch based on type of arg
    
    if(missing(datas)) {datas = guiGetSafe('datas');GUI = T} else {GUI=F}
    
    
    if(length(datas)<=1){warning('No data processed yet')} else {   
      
      datas$even=even
      
      if (Data.Type == 'Combined.Analysis')
      {
        class(datas) <- 'combined'
      } else if (Data.Type=='Fatty.Acid.Profiles')
      {
        class(datas) <- 'FA'
      } else if (Data.Type=='Stable.Isotopes')
      {
        class(datas) <- 'SI'
      }
      
      Covs = guiGetSafe('Covs')
      
      if(any(is.na(Covs)) & Analysis.Type == 'Analysis.with.Covariates'){stop('analysis with covariates selected, but no covariates entered.')}
      
      outputs <- switch(Analysis.Type,
                        Population.proportions = .Poppropanalysis(datas,nIter,nBurnin,nChains,nThin,plott),
                        Individual.proportions = .PopandIndprops(datas,nIter,nBurnin,nChains,nThin,plott),
                        Analysis.with.Covariates = .AnalysiswithCov(datas,Covs,nIter,nBurnin,nChains,nThin,plott)
      )
      guiSet('MCMCout',outputs)
      if(GUI==F) return(outputs)
    }
  }

jagger <- function(sysfile,nBurnin,nIter,nThin,i){
  
  load('jagsdata.Rdata')
  JM <- jags.model(file=sysfile,data=jagsdata,inits=list(.RNG.name="base::Mersenne-Twister",.RNG.seed=i),n.chains=1)
  
  cat('\n','proceeding to burn-in phase','\n')
  update(JM,n.iter=nBurnin)
  cat('\n','sampling from parameter distributions','\n')
  res<- coda.samples(model=JM,variable.names='prop',n.iter=nIter,thin=nThin)
  
  return(res)
   
}

spawn <- function(sysfile,nBurnin,nIter,nThin){
  
  options(useFancyQuotes=F)
  # spawn seperate R slave processes for each chain and run jagger in them
  for (i in 1:nChains){
    outfile <- paste('MCout',i,'.Rdata',sep='')
        
    cmd <- paste("Rscript -e 'library(FASTIN,quietly=T,verbose=F);options(warn=-1);res",i,"<- jagger(",dQuote(sysfile),",",eval(nBurnin),",",eval(nIter),",",eval(nThin),",",i,");save.image(file=",dQuote(outfile),");print(",dQuote("all done"),")'",sep='')
    
    if(i<nChains)
    {
      system(cmd,wait=F)
    }
    else {
      system(cmd,wait=T)
    }
    
  }
  
  res <- eval(parse(text=eval(load('MCout1.Rdata'))))
  
  if (nChains>1) {
    for(i in 2:nChains)
      res[[i]] <- eval(parse(text=eval(load(paste('MCout',i,'.Rdata',sep='')))))[[1]]
  }
  
  return(res)
  
}

addcovs <- function(Grps=NULL,Covs=NULL){
    
    if (length(Covs)>0 & length(Grps)==0)
    {
      Covs <- cbind(rep(1,nrow(Covs)),Covs)
      n.covs <- ncol(Covs)
      guiSet('Covs',Covs)
    } else if (length(Covs)==0 & length(Grps)>0) 
    {
      Grp.names <- unlist(unique(Grps)) 
      
      for (i in 1:ncol(Grps)){
        vg <- as.vector(Grps[,i])
        Grps[,i] <- as.factor(vg)
      }
      
      Covs <- model.matrix(attr(model.frame(1:nrow(Grps)~.,data=Grps),'terms'),data=Grps)[,]
      colnames(Covs) <- Grp.names[length(Grp.names):1]
      guiSet('Covs',Covs)
      
    } else if (length(Covs)>0 & length(Grps)>0) 
    {
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

diags <- function(MCMCout=NULL,accuracy=0.01,proba=0.95,quant=0.025){
    # do RL diag

    cat('\n','\n')
    cat('#################################','\n')
    cat('Raftery-Lewis diagnostics','\n')
    cat('#################################','\n','\n')
    
    #if(length(MCMCout)>1){
        class(MCMCout) <- 'mcmc.list'
        rd <- raftery.diag(MCMCout, q=quant, r=accuracy, s=proba)
   # } else {
    #    class(MCMCout) <- 'mcmc'
     #   rd <- raftery.diag(MCMCout, q=quant, r=accuracy, s=proba)
      #  rd <- list(rd)
    #}
   
    print(rd)
    if(rd[[1]]$resmatrix[1]=='Error'){
    cat('\n','Based on these diagnostics you should repeat the pilot MCMC with at least ',rd[[1]]$resmatrix[2],' \n','iterations to calculate diagnostics',' \n',' \n')
   } else{
    rit <- max(unlist(lapply(rd,function(x){x$resmatrix[,2]})))
    thin <- max(unlist(lapply(rd,function(x){x$resmatrix[,4]})))
    cat('\n','Based on these diagnostics you should repeat the MCMC with ',rit,' iterations','\n','and a thinning interval of ',round(thin),' ,if these values are higher than the values','\n','used to produce these diagnostics','\n','\n')
    }

    # do GR diag
    if(length(rd)>1){
       cat('\n','\n')
       cat('#################################','\n')
       cat('Gelman-Rubin diagnostics','\n')
       cat('#################################','\n','\n')

       print(gelman.diag(MCMCout,transform=T))
        
       cat('\n','Both univariate upper C.I. and multivariate psrf','\n','should be close to 1 if the chains converged','\n','\n','\n')
    }
    
}

FASTIN <- function(){
  #require(tcltk)   # this is needed - but leads to crashes...
  #require(fgui)
  
  # gui helper functions
    run_MCMC <- function(nIter=10000,nBurnin=1000,nChains=1,nThin=10,Data.Type='Fatty.Acid.Profiles',Analysis.Type='Population.proportions',even=0.1){
    # have three types here: FA, SI and combined, then methods dispatch based on type of arg
    
    datas = guiGetSafe('datas')    
    
    if(length(datas)<=1){stop('No data processed yet')} else {   
      
      datas$even=even
      
      if (Data.Type == 'Combined.Analysis')
      {
        class(datas) <- 'combined'
      } else if (Data.Type=='Fatty.Acid.Profiles')
      {
        class(datas) <- 'FA'
      } else if (Data.Type=='Stable.Isotopes')
      {
        class(datas) <- 'SI'
      }
      
      Covs = guiGetSafe('Covs')
      
      if(any(is.na(Covs)) & Analysis.Type == 'Analysis.with.Covariates'){stop('analysis with covariates selected, but no covariates entered.')}
      
      outputs <- switch(Analysis.Type,
                        Population.proportions = .Poppropanalysis(datas,nIter,nBurnin,nChains,nThin),
                        Individual.proportions = .PopandIndprops(datas,nIter,nBurnin,nChains,nThin),
                        Analysis.with.Covariates = .AnalysiswithCov(datas,Covs,nIter,nBurnin,nChains,nThin)
      )
      guiSet('MCMCout',outputs)
     
    }
  }
  
    
    pnorm_even <- function(even=0.1){p=2*(1-pnorm(log(95)/2,0,sqrt(1/even)));return(p)}
  
    dispsummaries <- function(Display.Summary=NULL){output <- guiGetSafe('MCMCout') ; 
                              if(class(output)=='pop_props'|class(output)=='ind_props'|class(output)=='cov_props')
                              {summary(output)}
  }

    dispdiags <- function(accuracy=0.01,proba=0.95,quant=0.025){
        output <- guiGetSafe('MCMCout') 
        if(class(output)=='pop_props'|class(output)=='ind_props'|class(output)=='cov_props')
            {
                diags(MCMCout=output,accuracy,proba,quant)
            }
    }
    
    plotoutputs<- function(){output <- guiGetSafe('MCMCout') 
                           if(class(output)=='pop_props'|class(output)=='ind_props'|class(output)=='cov_props')
                           {plot(output)}
  }
  
    saveoutputs <- function(Path="MCMCout.Rdata"){output <- guiGetSafe('MCMCout');save(output,file=Path)}
  
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
      Covnames <- names(Covs)
      Grps <- read.csv(Groups,header=T)
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
  
    SaveData <- function(Path="datas.Rdata"){datas <- guiGetSafe('datas');stopifnot(length(datas)>1);save(datas,file=Path)}
    LoadData <- function(Path=NULL){ load(Path); guiSet('datas',datas)}
  
    guiSet( "LIST_WIDTH", 50)
    guiSet( "ENTRY_WIDTH", 10)
  
  #Gui - tried some meaningful indentation here, but still not quite right...
    output <- gui(.fastin, title = 'FASTIN main menu',
                 argCommand=list(add.covs=guiNestedF(addcovs,"add.covs",  argFilter=list(Groups="{{} {.csv}}",Covariates="{{} {.csv}}"),
                                                     argText=c(Covariates = "Add Covariates (optional)",Groups = "Add Groups (optional)"),cancelButton=F,exec='Add'),
                     Save.Outputs=guiNestedF(saveoutputs,"Save.Outputs",argText = list(Path='Choose filename'),cancelButton=F,exec='save'),
                     Save.Data=guiNestedF(SaveData,"Save.Data",argText = list(Path='Choose filename'),cancelButton=F,exec='save'),
                     Load.Data=guiNestedF(LoadData,"Load.Data",argFilter=list(Path= "{{} {.Rdata}}"),cancelButton=F,exec='load'),
                     Display.Summaries=dispsummaries,
                     groupings = resetsc,
                     Plot.Outputs=plotoutputs,
                     SI.data = guiNestedF(addSI,"SI.data",
                                                      argFilter=list(predators.SI= "{{} {.csv}}",preys.SI= "{{} {.csv}}",Frac.Coeffs.mean="{{} {.csv}}",Frac.Coeffs.var="{{} {.csv}}"), 
                                                      title = 'Stable Isotope data entry form',
                                                      exec = "Add Stable Isotope data", closeOnExec = TRUE,
                                                      output = NULL,
                                                      argText=c(predators.SI="Load predator(s) stable isotope data (csv)",
                                                                preys.SI="Load prey stable isotope data (csv)",
                                                                Frac.Coeffs.mean="Load prey and stable isotope specific fractionation means (csv)",
                                                                Frac.Coeffs.var="Load prey and stable isotope specific fractionation sd (csv)",
                                                                FC.mean="Mean fractionation coefficients (use R's c() notation)",
                                                                FC.var="SD of fractionation coefficients (use R's c() notation)",
                                                                R.diag.SI = "Diagonal of the prior for predator (co)-variance matrix (>0, smaller value is less informative)"),
                                                      cancelButton=F), 
                     FA.data = guiNestedF(addFA,"FA.data",argFilter=list(predators.FA= "{{} {.csv}}",preys.FA= "{{} {.csv}}",fat.conts = "{{} {.csv}}",Conv.Coeffs.mean="{{} {.csv}}",Conv.Coeffs.var="{{} {.csv}}"), 
                                                      #argEdit = list(CC.mean=NULL,CC.var=NULL,R.diag=NULL),
                                                      title = 'Fatty Acid Profile data entry form',
                                                      argText=c(predators.FA="Load predator(s) fatty acid data (csv)",
                                                                preys.FA="Load prey fatty acid data (csv)",
                                                                Conv.Coeffs.mean="Load prey and fatty acid specific conversion means (csv)",
                                                                Conv.Coeffs.var="Load prey and fatty acid specific conversion sd (csv)",
                                                                fat.conts = "Prey fat content (csv)",
                                                                FC.mean="Mean fat content (use R's c() notation)",
                                                                FC.var="SD of fat content (use R's c() notation)", 
                                                                CC.mean="Mean conversion coefficients (use R's c() notation)",
                                                                CC.var="SD of conversion coefficients (use R's c() notation)",
                                                                R.diag = "Diagonal of the prior for predator (co)-variance matrix (>0, smaller value is less informative)"),
                                                      exec = "Add Fatty Acid data", cancelButton=F,closeOnExec = TRUE,output = NULL),  
                     MCMC = guiNestedF(run_MCMC,"MCMC",
                                                   title = 'Run FASTIN analysis',
                                                   exec = "Run MCMC", closeOnExec = TRUE,
                                                   output = NULL,
                                                   argCommand=list(even = guiNestedF(pnorm_even,'even',exec='Set new prior',argSlider=list(even=c(0.001,2,0.001),output=c(0,1,0.001)),callback=guiExec,argText=c(output = "Probability that prop(prey x) = 0.95*prop(preys other than x)"))),
                                                   argOption = list(Data.Type = c('Fatty.Acid.Profiles','Stable.Isotopes','Combined.Analysis'), defaultChoice=1,Analysis.Type = c('Population.proportions','Individual.proportions','Analysis.with.Covariates'), defaultChoice=1),                                                   
                                                   argText=c(nIter = 'Number of MCMC iterations',
                                                             nBurnin = 'Number of MCMC iterations to discard (burn-in)',
                                                             nChains = 'Number of Markov Chains',
                                                             nThin = 'Thinning interval of MCMC chains',
                                                             even= 'Prior eveness of proportions'
                                                   ),cancelButton=F),
                     disp.diags = guiNestedF(dispdiags,"disp.diags",
                                                   title = 'Show Convergence Diagnositcs', 
                                                   exec = "Display diagnostics in the R console",
                                                   argText=c(proba = 'probability',
                                                       quant='quantile to estiamte'
                                                   ),cancelButton=F)
                 ),                                 
                 argText=c(SI.data='Add Stable Isotope data', 
                           FA.data='Add Fatty Acid profiles',
                           groupings = 'Reset previous prey grouping',
                           Save.Data= "Save data",
                           Load.Data = 'Load previously saved dataset',
                           add.covs='Add covariates and/or groups',
                           MCMC = 'Run Bayesian analysis (MCMC)',
                     disp.diags = 'Show convergence diagnostics'
                 ),
                 exec=NULL,output=NULL,argGridOrder=c(1,1,1,2,2,2,3,3,4,4,4), argGridSticky=rep("w",length(formals(.fastin)))
  )
  #detach("package:fgui", unload=TRUE)
  #detach("package:tcltk", unload=TRUE)
  # even=0.1
  # guiSet("even",even)
  
  return(output)
  
}
