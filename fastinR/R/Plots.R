#' Plot Non-metric multidimensional scaling (NMDS) or (clr transformed) diet data
#' 
#' Diet data produced by \code{\link{add_FA}} and/or \code{\link{add_FA}} is projected onto 2 dimesions using NMDS, where Fatty Acid data is first clr transformed and then concatenated with SI data to produce the plot.
#' @param datas Diet data produced by \code{\link{add_FA}} and/or \code{\link{add_SI}}
#' @seealso \code{\link{add_FA}},\code{\link{add_SI}},\code{\link{select_vars}},\code{\link{var_select_plot}},\code{\link{run_MCMC}}
#' @author Philipp Neubauer
#' @references Neubauer,.P. and Jensen, O.P. (in prep)
#' @export
dataplot <- function(datas=NULL){
  
  # check if GUI is being used
  if(exists('GUI',envir=.GlobalEnv)){
    GUI <- get('GUI',envir=.GlobalEnv)
    if(GUI) datas <- guiGetSafe('datas')
  } else {
    GUI=F
  }
  
  preya <- {};preya.SI <- {}
  preda <- {};preda.SI <- {}
  if(!is.null(datas$datas.FA$preys)) {
    preya=cbind(preya,clr(datas$datas.FA$preys*((datas$datas.FA$mean_c/datas$datas.FA$tau_c)[datas$prey.ix,])))
    preda=cbind(preda,clr(datas$datas.FA$preds.FA))     
  }
  if(!is.null(datas$datas.SI$preys.SI)) {
    if (ifelse(!is.null(datas$prey.ix),all(datas$prey.ix == datas$prey.ix.SI),TRUE)) {
      preya=cbind(preya,as.matrix(datas$datas.SI$preys.SI+datas$datas.SI$mean_cs[datas$prey.ix.SI,]))
      preda=cbind(preda,as.matrix(datas$datas.SI$preds.SI))
    } else {
      preya.SI=cbind(preya.SI,as.matrix(datas$datas.SI$preys.SI+datas$datas.SI$mean_cs[datas$prey.ix.SI,]))
      preda.SI=cbind(preda.SI,as.matrix(datas$datas.SI$preds.SI))
      dista.SI <- dist(rbind(preya.SI,preda.SI))
    }
  }  
  
  names(preya)  <- names(preda)
  
  dista <- dist(rbind(preya,preda))
  if (!is.null(preya.SI)) dista = sqrt(dista^2 + dista.SI^2)
  mds <- metaMDS(dista,trymax=100)
  
  externalDevice<-FALSE
  if (!is.function(options()$device) & GUI==T){
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
  points(pl,'sites',pch=cbind(as.numeric(factor(datas$prey.ix,levels=unique(datas$prey.ix))),rep(16,datas$n.preds)),col=cbind(1+as.numeric(factor(datas$prey.ix,levels=unique(datas$prey.ix))),rep(1,datas$n.preds)))
  legend('topleft',c('Predators',unique(datas$prey.ix)),xpd=T,pch=c(16,1:datas$n.preys),col=c(1,2:(datas$n.preys+1)))
  
}

#' @name plot
#' @title Plot density or denstrip plots of posterior diet proportions
#' @description A pairwise plot to inspect correlations between prey items in the MCMC estiamtion, and a density/denstrip plot of diet proportion posterior distributions. To be run with output from \code{\link{run_MCMC}}
#' @S3method plot pop_props
#' @S3method plot cov_props
#' @S3method plot cov_props
#' @param x MCMC output from \code{\link{run_MCMC}}, containing diet proportion MCMC chains
#' @param save Either a string to be used as prefix for saved plots, or FALSE for disabling saving to file.
#' @param density If TRUE (default), density plots are drawn, if FALSE, denstrip plots drawn instead.
#' @details If plots are saved they are not drawn at the same time. That may change in the future...
#' @references Neubauer.P. and Jensen, O.P. (in prep)
#' @author Philipp Neubauer
#' @seealso \code{\link{run_MCMC}},\code{\link{diags}}
NULL

#' @method plot pop_props
#' @export
plot.pop_props <- function(x,save="fastinR_MCMC_",density=T,...){
  
  if(save!=F){sava <- menu(title='save plots?',choices = c('yes','no'),graphics=T)}else{sava=0}
  
  # check gui
  if(exists('GUI',envir=.GlobalEnv)) {
    GUI <- get('GUI',envir=.GlobalEnv)
  } else {
    GUI <- F
  }
  
  preya.names <- x$prey.names
  
  outs={}
  for (k in 1:x$nChains){
    outp <- matrix(unlist(x[[k]]),ncol=ncol(x[[1]]),byrow=F)
    outs <- rbind(outs,outp)
  }
  
  outs <- as.data.frame(outs)
  colnames(outs) <- preya.names
  
  
  if (sava==1) {
    pdf(paste(save,"correlations.pdf",sep=''))
  } else {
          if (!is.function(options()$device) & GUI==T & GUI==T){
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
              
            }
          }
  }
plot(outs,main='Correlation of proportion estimates')
if (sava==1) dev.off()


if(density==F){
  mp=reshape::melt.data.frame(outs)
  if (sava==1){pdf(paste(save,"pop_proportions.pdf",sep=''))}else{par(ask=T)}
  rpp = bwplot(variable~value,data=mp,xlim=c(0,1),xlab=list(label=expression("Posterior distribution of diet proportions"),cex=1)
               ,panel = function(x, y) { 
                 grid::grid.segments(1,0,0,0)
                 xlist <- split(x, factor(y))
                 for (i in seq(along=xlist))
                   denstrip::panel.denstrip(x=xlist[[i]], at=i,colmax=(1),ticks=quantile(xlist[[i]],c(0.025,0.975)),mticks=quantile(xlist[[i]],0.5),tlen=1.2,mlen=1.2,nwd=1.5)
               },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T,labels=preya.names)))
  print(rpp)
  par(ask=F)
  trellis.focus("panel", 1, 1,highlight=F)
  panel.abline(h=seq(1.5,ncol(outs),1),col=1,lty=2)
  trellis.unfocus()
  if (sava==1) dev.off()
} else {
  
  xx<-reshape::melt.data.frame(outs)
  
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
  ggs_density <- function(D, pos=0,title) {
    # Manage subsetting a family of parameters
    
    dims = dim(D)
    D<-as.data.frame(D)
    
    f <- ggplot(aes(x=value, fill=as.factor(variable), col=as.factor(variable)),data=D)
    
    f <- f + geom_density(alpha=0.3,adjust=1.5,data=D) + 
      geom_vline(aes(xintercept=pos),linetype=1:length(pos),data=as.data.frame(pos),colour=1)+
      geom_vline(aes(xintercept=0),colour=1)+
      geom_hline(aes(yintercept=0),colour=1)+
      scale_x_continuous(limits=c(0,1))+
      xlab('Proportion')+
      ylab('Posterior Density')+
      ggtitle(title)+    
      theme(panel.grid.major=element_line(colour = NA),
            panel.grid.minor=element_line(colour = NA),
            panel.background=element_rect(fill=NA),
            #legend.position = "none",
            legend.title=element_blank(),
            axis.ticks = element_line(colour=1),
            axis.text = element_text(colour=1,size = rel(0.8)),
            axis.title = element_text(size = rel(1)))+
      geom_rug(alpha=0.1)
    
    f
  }
  
  if(sava==1) {pdf(paste(save,'_densities.pdf',sep=''))}else{par(ask=T)}
    multiplot(ggs_density(xx,title=''))
    
    if(sava==1) dev.off()
    
  }
# turn off external device if using one
}

#' @method plot ind_props
#' @export
plot.ind_props <- function(x,save="fastinR_MCMC_",density=T,...){
  
  if(save!=F){sava <- menu(title='save plots?',choices = c('yes','no'),graphics=T)}else{sava=0}
  
  # check gui
  if(exists('GUI',envir=.GlobalEnv)) {
    GUI <- get('GUI',envir=.GlobalEnv)
  } else {
    GUI <- F
  }
  
  preya.names <- x$prey.names
  
  outs={}
  for (k in 1:x$nChains){
    outp <- matrix(unlist(x[[k]]),ncol=ncol(x[[1]]),byrow=F)
    outs <- rbind(outs,outp)
  }
  outs <- as.data.frame(outs)
  
  colnames(outs) <- colnames(x[[1]])
  popix <- grep('pop',colnames(outs))
  names(outs)[popix] <- preya.names
  
  if (sava==1)  {
    pdf(paste(save,"correlations.pdf",sep=''))
  } else {
          if (!is.function(options()$device) & GUI==T){
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
              
            }
          }
  }
  plot(outs[,popix],main='Correlation of proportion estimates')
  if (sava==1) dev.off()
  
  # population proportions
  if(density==F){
    # draw popualtion posteriors
    mp=reshape::melt.data.frame(outs[,popix])
    if (sava==1){pdf(paste(save,"pop_proportions.pdf",sep=''))}else{par(ask=T)}
    rpp = bwplot(variable~value,data=mp,xlim=c(0,1),xlab=list(label=expression("Posterior distribution of population diet proportions"),cex=1)
                 ,panel = function(x, y) { 
                   grid::grid.segments(1,0,0,0)
                   xlist <- split(x, factor(y))
                   for (i in seq(along=xlist))
                     denstrip::panel.denstrip(x=xlist[[i]], at=i,colmax=(1),ticks=quantile(xlist[[i]],c(0.025,0.975)),mticks=quantile(xlist[[i]],0.5),tlen=1.2,mlen=1.2,nwd=1.5)
                 },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T,labels=preya.names)))
    print(rpp)
    par(ask=F)
    trellis.focus("panel", 1, 1,highlight=F)
    panel.abline(h=seq(1.5,ncol(outs),1),col=1,lty=2)
    trellis.unfocus()
    if (sava==1) dev.off()
    
  } else
  {
    multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
      
      
      # Make a list from the ... arguments and plotlist
      plots <- c(list(...), plotlist)
      
      numPlots = length(plots)
      
      # If layout is NULL, then use 'cols' to determine layout
      if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
      }
      
      if (numPlots==1) {
        print(plots[[1]])
        
      } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
          # Get the i,j matrix positions of the regions that contain this subplot
          matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
          
          print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                          layout.pos.col = matchidx$col))
        }
      }
    }
    
    xx<-reshape::melt.data.frame(outs[,popix])
    
    ggs_density <- function(D, pos=0,title) {
      # Manage subsetting a family of parameters
      
      dims = dim(D)
      D<-as.data.frame(D)
      
      f <- ggplot(aes(x=value,fill=as.factor(variable), col=as.factor(variable)),data=D)
      
      f <- f + geom_density(alpha=0.3,adjust=1.5,data=D) + 
        geom_vline(aes(xintercept=pos),linetype=1:length(pos),data=as.data.frame(pos),colour=1)+
        geom_vline(aes(xintercept=0),colour=1)+
        geom_hline(aes(yintercept=0),colour=1)+
        scale_x_continuous(limits=c(0,1))+
        xlab('Proportion')+
        ylab('Posterior Density')+
        ggtitle(title)+    
        theme(panel.grid.major=element_line(colour = NA),
              panel.grid.minor=element_line(colour = NA),
              panel.background=element_rect(fill=NA),
              #legend.position = "none",
              legend.title=element_blank(),
              axis.ticks = element_line(colour=1),
              axis.text = element_text(colour=1,size = rel(0.8)),
              axis.title = element_text(size = rel(1)))+
        geom_rug(alpha=0.1)
      
      f
    }
    
    if(sava==1) {pdf(paste(save,'_densities.pdf',sep=''))}else{par(ask=T)}
    multiplot(ggs_density(xx,title='Population proportions'))
    if(sava==1) dev.off()
    
  }
  
  # draw individual posteriors
  
  indix <- (max(popix)+1):ncol(outs)
  n.preys <- length(popix)
  n.preds <- length(indix)/n.preys
  orders <- order(apply(outs[,indix[1:n.preds]],2,median))
  this.prey <- which.max(tapply(apply(outs[,indix],2,median),rep(1:n.preys,each=n.preds),median))
  orders <- order(apply(outs[,indix[(1:n.preds)+(this.prey-1)*n.preds]],2,median))
  indix <- rep(orders,each=n.preys)+c(0,n.preds*1:(n.preys-1))+max(popix)
  indix2 <- order(names(outs[,indix]))
  
  if(density==F){
    mp=reshape::melt.data.frame(outs[,indix])
    
    if (sava==1){pdf(paste(save,"ind_proportions.pdf",sep=''))}else{par(ask=T)}
    rpp = bwplot(variable~value,data=mp,xlim=c(0,1),xlab=list(label=expression("Posterior distribution of individual diet proportions"),cex=1),auto.key = list(text=preya.names, points = F,
                                                                                                                                                               rectangles = F,
                                                                                                                                                               lines = FALSE,
                                                                                                                                                               col=1:n.preds,alpha=1,space='right'),
                 ,panel = function(x, y) { 
                   grid::grid.segments(1,0,0,0)
                   xlist <- split(x, factor(y))           
                   for (i in seq(along=xlist)){
                     denstrip::panel.denstrip(x=xlist[[i]], at=i,colmin="transparent",colmax=((i-1)%%n.preys)+1,ticks=quantile(xlist[[i]],c(0.025,0.975)),mticks=quantile(xlist[[i]],0.5),tlen=1.2,mlen=1.2,nwd=1.5,width=1/n.preds*n.preys)
                   }
                 },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T,labels=sprintf(paste0("Predator %d"),orders),at=seq(n.preys/2+0.5,n.preys*n.preds,n.preys))))
    print(rpp)
    par(ask=F)
    trellis.focus("panel", 1, 1)
    panel.abline(h=seq(n.preys+0.5,n.preys*n.preds,n.preys),col=1,lty=2)
    trellis.unfocus()
    if (sava==1) dev.off()
    
  } else {
    
    colnames(outs) <-NA
    outr <- data.frame()
    a=1
    for (i in 1:n.preds){
      outr <- rbind(outr,cbind(rep(i,nrow(outs)),outs[,indix[a:(a+n.preys-1)]]))
      a = a+n.preys
    }
    
    colnames(outr)[2:(n.preys+1)] <- preya.names    
    outr[,1] <- as.factor(outr[,1])
    attr(outr[,1],'levels') <- sprintf(paste0("Predator %d"),orders)
    
    xx<-reshape::melt.data.frame(outr)
    
    yy<-split(xx,factor(xx[,1]))
    
    ggs_density <- function(D, pos=0,title) {
      # Manage subsetting a family of parameters
      
      
      dims = dim(D)
      #   if (is.null(dims)){
      D<-as.data.frame(D)
      # }
      
      # Plot
      
      f <- ggplot(aes(x=value,fill=as.factor(variable), col=as.factor(variable)),data=D)
      
      f <- f + geom_density(alpha=0.3,adjust=1.5,data=D,trim=F)+ 
        geom_vline(aes(xintercept=pos),linetype=1:length(pos),data=as.data.frame(pos),colour=1)+
        geom_vline(aes(xintercept=0),colour=1)+
        geom_hline(aes(yintercept=0),colour=1)+
        scale_x_continuous(limits=c(0,1))+
        xlab('Proportion')+
        ylab('Posterior Density')+
        ggtitle(title)+    
        theme(panel.grid.major=element_line(colour = NA),
              panel.grid.minor=element_line(colour = NA),
              panel.background=element_rect(fill=NA),
              #legend.position = "none",
              legend.title=element_blank(),
              axis.ticks = element_line(colour=1),
              axis.text = element_text(colour=1,size = rel(0.8)),
              axis.title = element_text(size = rel(1)))+
        geom_rug(alpha=0.1)
      
      return(f)
    }
    multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
      # Make a list from the ... arguments and plotlist
      plots <- c(list(...), plotlist)
      
      numPlots = length(plots)
      
      # If layout is NULL, then use 'cols' to determine layout
      if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
      }
      
      if (numPlots==1) {
        print(plots[[1]])
        
      } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
          # Get the i,j matrix positions of the regions that contain this subplot
          matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
          
          print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                          layout.pos.col = matchidx$col))
        }
      }
    }
    g={}
    for (i in seq(along=yy)) g[[i]] <- ggs_density(yy[[i]],title=levels(outr[,1])[i])
    
    if(sava==1) pdf('ind_densities.pdf')
    multiplot(plotlist=g,cols=2)
    if(sava==1) dev.off()
  }
}

#' @method plot cov_props
#' @export
plot.cov_props <- function(x,save="fastinR_MCMC_",density=T,...){
  
  if(save!=F){sava <- menu(title='save plots?',choices = c('yes','no'),graphics=T)}else{sava=0}
   
  # check gui
  if(exists('GUI',envir=.GlobalEnv)) {
    GUI <- get('GUI',envir=.GlobalEnv)
  } else {
    GUI <- F
  }
  
  preya.names <- x$prey.names
  Covs <- x$Covs
  
  outs={}
  for (k in 1:x$nChains){
    outp <- matrix(unlist(x[[k]]),ncol=ncol(x[[1]]),byrow=F)
    outs <- rbind(outs,outp)
  }
  outs <- as.data.frame(outs)
  
  colnames(outs) <- colnames(x[[1]])
  
  
  cidx <- apply(Covs,2,function(x){any(x!=0 & x!=1)})
  #number of groups and covariates
  nGr <- sum(cidx==F)
  Gridx <- which(cidx==F)
  
  popix <- grep('pop',colnames(outs))
  n.preys <- length(popix)/nGr
  names(outs)[popix] <- preya.names
  
  
  betaix <- grep('beta',colnames(outs))
  
  this.eff <- outs[,popix[1:(nGr*n.preys)]]
  k=1
  if (sava==1)  {
    pdf(paste(save,"correlations.pdf",sep=''))
  } else {
          if (!is.function(options()$device) & GUI==T){
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
              
            }
          }
  }
  for (n in 1:nGr) {colnames(this.eff)[k:(k+n.preys-1)] <- sprintf(paste('Group',n,'/ %s'),preya.names); k=k+n.preys}
  plot( this.eff,main='Correlation of proportion estimates')
  if (sava==1) dev.off()
  
  # population proportions
  if(density==F){
    # draw popualtion posteriors
    mp=reshape::melt.data.frame(this.eff)
    if (sava==1){pdf(paste(save,"pop_proportions.pdf",sep=''))}else{par(ask=T)}
    rpp = bwplot(variable~value,data=mp,xlim=c(0,1),xlab=list(label=expression("Posterior distribution of population diet proportions"),cex=1)
                 ,panel = function(x, y) { 
                   grid::grid.segments(1,0,0,0)
                   xlist <- split(x, factor(y))
                   for (i in seq(along=xlist))
                     denstrip::panel.denstrip(x=xlist[[i]], at=i,colmax=((i-1)%%n.preys)+1,ticks=quantile(xlist[[i]],c(0.025,0.975)),mticks=quantile(xlist[[i]],0.5),tlen=1.2,mlen=1.2,nwd=1.5)
                 },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T)))
    print(rpp)
    par(ask=F)
    trellis.focus("panel", 1, 1)
    panel.abline(h=seq(n.preys+0.5,ncol(this.eff),n.preys),col=1,lty=2)
    trellis.unfocus()
    if (sava==1) dev.off()
  } else
  {
    multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
      
      # Make a list from the ... arguments and plotlist
      plots <- c(list(...), plotlist)
      
      numPlots = length(plots)
      
      # If layout is NULL, then use 'cols' to determine layout
      if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
      }
      
      if (numPlots==1) {
        print(plots[[1]])
        
      } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
          # Get the i,j matrix positions of the regions that contain this subplot
          matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
          
          print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                          layout.pos.col = matchidx$col))
        }
      }
    }
    
    xx<-reshape::melt.data.frame(this.eff)
    nf <- sapply(strsplit(as.character(xx[,1]),' / '),I)[1,]
    yy <- split(reshape::melt.data.frame(this.eff),factor(nf))
    ggs_density <- function(D, pos=0,title) {
      # Manage subsetting a family of parameters
      
      dims = dim(D)
      D<-as.data.frame(D)
      
      f <- ggplot(aes(x=value,fill=as.factor(variable), col=as.factor(variable)),data=D)
      
      f <- f + geom_density(alpha=0.3,adjust=1.5,data=D) + 
        geom_vline(aes(xintercept=pos),linetype=1:length(pos),data=as.data.frame(pos),colour=1)+
        geom_vline(aes(xintercept=0),colour=1)+
        geom_hline(aes(yintercept=0),colour=1)+
        scale_x_continuous(limits=c(0,1))+
        xlab('Proportion')+
        ylab('Posterior Density')+
        ggtitle(title)+    
        theme(panel.grid.major=element_line(colour = NA),
              panel.grid.minor=element_line(colour = NA),
              panel.background=element_rect(fill=NA),
              #legend.position = "none",
              legend.title=element_blank(),
              axis.ticks = element_line(colour=1),
              axis.text = element_text(colour=1,size = rel(0.8)),
              axis.title = element_text(size = rel(1)))+
        geom_rug(alpha=0.1)
      
      f
    }
    
    g={}
    for (i in seq(along=yy)) g[[i]] <- ggs_density(yy[[i]],title=names(yy)[i])
    
    
    if(sava==1) {pdf(paste(save,'_densities.pdf',sep=''))}else{par(ask=T)}
    multiplot(plotlist=g)
    if(sava==1) dev.off()
    
  }
  
  # draw individual posteriors
  
  indix <- (max(popix)+1):ncol(outs)
  n.preds <- length(indix)/n.preys
  # check which prey to ordewr by (highest proportion)
  this.prey <- which.max(tapply(apply(outs[,indix],2,median),rep(1:n.preys,each=n.preds),median))
  orders <- order(apply(outs[,indix[(1:n.preds)+(this.prey-1)*n.preds]],2,median))
  indix <- rep(orders,each=n.preys)+c(0,n.preds*1:(n.preys-1))+max(popix)
  indix2 <- order(names(outs[,indix]))
  
  if(density==F){
    mp=reshape::melt.data.frame(outs[,indix])
    if (sava==1){pdf(paste(save,"ind_proportions.pdf",sep=''))}else{par(ask=T)}
    rpp = bwplot(variable~value,data=mp,xlim=c(0,1),xlab=list(label=expression("Posterior distribution of individual diet proportions"),cex=1),auto.key = list(text=preya.names, points = F,
                                                                                                                                                               rectangles = F,
                                                                                                                                                               lines = FALSE,
                                                                                                                                                               col=1:n.preys,alpha=1,space='right'),
                 ,panel = function(x, y) { 
                   grid::grid.segments(1,0,0,0)
                   xlist <- split(x, factor(y))
                   k=0
                   for (i in seq(along=xlist)){
                     k=k+1
                     denstrip::panel.denstrip(x=xlist[[i]], at=i,colmin="transparent",colmax=((i-1)%%n.preys)+1,ticks=quantile(xlist[[i]],c(0.025,0.975)),mticks=quantile(xlist[[i]],0.5),tlen=1.2,mlen=1.2,nwd=1.5,width=1/n.preds*n.preys)
                   }
                 },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T,labels=sprintf(paste0("Predator %d"),orders),at=seq(n.preys/2+0.5,n.preys*n.preds,n.preys))))
    print(rpp)
    par(ask=F)
    trellis.focus("panel", 1, 1,highlight=F)
    panel.abline(h=seq(n.preys+0.5,n.preys*n.preds,n.preys),col=1,lty=2)
    trellis.unfocus()
    if (sava==1) dev.off()
    
  } else {
    
    colnames(outs) <-NA
    outr <- data.frame()
    a=1
    for (i in 1:n.preds){
      outr <- rbind(outr,cbind(rep(i,nrow(outs)),outs[,indix[a:(a+n.preys-1)]]))
      a = a+n.preys
    }
    
    colnames(outr)[2:(n.preys+1)] <- preya.names    
    outr[,1] <- as.factor(outr[,1])
    attr(outr[,1],'levels') <- sprintf(paste0("Predator %d"),orders)
    
    xx<-reshape::melt.data.frame(outr)
    
    yy<-split(xx,factor(xx[,1]))
    
    ggs_density <- function(D, pos=0,title) {
      # Manage subsetting a family of parameters
      dims = dim(D)
      #   if (is.null(dims)){
      D<-as.data.frame(D)
      # }
      
      # Plot
      
      f <- ggplot(aes(x=value,fill=as.factor(variable), col=as.factor(variable)),data=D)
      
      f <- f + geom_density(alpha=0.3,adjust=1.5,data=D,trim=F)+ 
        geom_vline(aes(xintercept=pos),linetype=1:length(pos),data=as.data.frame(pos),colour=1)+
        geom_vline(aes(xintercept=0),colour=1)+
        geom_hline(aes(yintercept=0),colour=1)+
        scale_x_continuous(limits=c(0,1))+
        xlab('Proportion')+
        ylab('Posterior Density')+
        ggtitle(title)+    
        theme(panel.grid.major=element_line(colour = NA),
              panel.grid.minor=element_line(colour = NA),
              panel.background=element_rect(fill=NA),
              #legend.position = "none",
              legend.title=element_blank(),
              axis.ticks = element_line(colour=1),
              axis.text = element_text(colour=1,size = rel(0.8)),
              axis.title = element_text(size = rel(1)))+
        geom_rug(alpha=0.1)
      
      return(f)
    }
    multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
      
      
      # Make a list from the ... arguments and plotlist
      plots <- c(list(...), plotlist)
      
      numPlots = length(plots)
      
      # If layout is NULL, then use 'cols' to determine layout
      if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
      }
      
      if (numPlots==1) {
        print(plots[[1]])
        
      } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
          # Get the i,j matrix positions of the regions that contain this subplot
          matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
          
          print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                          layout.pos.col = matchidx$col))
        }
      }
    }
    g={}
    for (i in seq(along=yy)) g[[i]] <- ggs_density(yy[[i]],title=levels(outr[,1])[i])
    
    if(sava==1) pdf('ind_densities.pdf')
    multiplot(plotlist=g,cols=2)
    if(sava==1) dev.off()
  }        
}

#' @name multiplot
#' @title Plot density or denstrip plots from multiple MCMC runs
#' @description A density/denstrip plot produced for multiple MCMC runs from \code{\link{run_MCMC}}. Useful to compare across runs with different data subsets or data types. The analyses needs to be of the same type, for now only Analysis.Type = Population.Proportions in \code{\link{run_MCMC}} works.
#' @param MCMCouts A (named) list of objects produced by  \code{\link{run_MCMC}}
#' @param density If TRUE (default), density plots are drawn, if FALSE, denstrip plots drawn instead. 
#' @param save Either a string to be used as prefix for saved plots, or FALSE for disabling saving to file. 
#' @details If plots are saved they are not drawn at the same time.
#' @references Neubauer.P. and Jensen, O.P. (in prep) 
#' @author Philipp Neubauer
#' @seealso \code{\link{run_MCMC}},\code{\link{diags}}
#' @export
multiplot <- function(MCMCouts,save="fastinR_Multiplot",density=T){
  
  if(save!=F){sava <- menu(title='save plots?',choices = c('yes','no'),graphics=T)}else{sava=0}
  
  # check gui
  if(exists('GUI',envir=.GlobalEnv)) {
    GUI <- get('GUI',envir=.GlobalEnv)
  } else {
    GUI <- F
  }
  
  if (sava==0) {
    
    if (!is.function(options()$device) & GUI==T){
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
        
      }
    }
  }
  
  preya.names <- (MCMCouts[[1]]$prey.names)
  n.preys <- length(preya.names)
    
  if(density==F){
    
    outs <- {}
    for (i in 1:length(MCMCouts)){
      
      x <- MCMCouts[[i]]
      outa={}
      for (k in 1:x$nChains){
        outp <- matrix(unlist(x[[k]]),ncol=ncol(x[[1]]),byrow=F)
        outa <- rbind(outa,outp)
      }
      outs <- cbind(outs,outa)
    }
    
    outs <- as.data.frame(outs)
    
    
    if(is.null(names(MCMCouts))){MCnames <- sprintf(paste('MCMC %s'),1:length(MCMCouts))}else{MCnames <-names(MCMCouts)}
    
    k=1
    for (n in 1:length(MCMCouts)) {colnames(outs)[k:(k+n.preys-1)] <- sprintf(paste(MCnames[n],'/ %s'),preya.names); k=k+n.preys}
    
    # draw popualtion posteriors
    mp=reshape::melt.data.frame(outs)
    if (sava==1) pdf(paste(save,"_proportions.pdf",sep=''))
    rpp = bwplot(variable~value,data=mp,xlim=c(0,1),xlab=list(label=expression("Posterior distribution of population diet proportions"),cex=1)
                 ,panel = function(x, y) { 
                   grid::grid.segments(1,0,0,0)
                   xlist <- split(x, factor(y))
                   for (i in seq(along=xlist))
                     denstrip::panel.denstrip(x=xlist[[i]], at=i,colmax=((i-1)%%n.preys)+1,ticks=quantile(xlist[[i]],c(0.025,0.975)),mticks=quantile(xlist[[i]],0.5),tlen=1.2,mlen=1.2,nwd=1.5)
                 },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T)))
    print(rpp)
    
    trellis.focus("panel", 1, 1,highlight=F)
    panel.abline(h=seq(n.preys+0.5,ncol(outs),n.preys),col=1,lty=2)
    trellis.unfocus()
    if (sava==1) dev.off()
  } else 
    {
    
    if(is.null(names(MCMCouts))){MCnames <- sprintf(paste('MCMC %s'),1:length(MCMCouts))}else{MCnames <-names(MCMCouts)}
    
    outs <- data.frame()
    for (i in 1:length(MCMCouts)){
      
      x <- MCMCouts[[i]]
      outa={}
      for (k in 1:x$nChains){
        outp <- matrix(unlist(x[[k]]),ncol=ncol(x[[1]]),byrow=F)
        outa <- rbind(outa,outp)
      }
      outs <- rbind(outs,cbind(rep(i,nrow(outa)),(outa)))
    }
    
    outs[,1] <- as.factor(outs[,1])
    attr(outs[,1],'levels') <- MCnames

    
    xx<-reshape::melt.data.frame(outs)
    
    yy<-split(xx,factor(xx$variable))
        
    ggs_density <- function(D, pos=0,title) {
      # Manage subsetting a family of parameters
      
      dims = dim(D)
      #   if (is.null(dims)){
      D<-as.data.frame(D)
      # }
      
      # Plot
      
      f <- ggplot(aes(x=value, fill=as.factor(V1), col=as.factor(V1)),data=D)
      
      f <- f + geom_density(alpha=0.3,adjust=1.5,data=D) + 
        geom_vline(aes(xintercept=pos),linetype=1:length(pos),data=as.data.frame(pos),colour=1)+
        geom_vline(aes(xintercept=0),colour=1)+
        geom_hline(aes(yintercept=0),colour=1)+
        scale_x_continuous(limits=c(0,1))+
        xlab('Proportion')+
        ylab('Posterior Density')+
        ggtitle(title)+    
        theme(panel.grid.major=element_line(colour = NA),
              panel.grid.minor=element_line(colour = NA),
              panel.background=element_rect(fill=NA),
              #legend.position = "none",
              legend.title=element_blank(),
              axis.ticks = element_line(colour=1),
              axis.text = element_text(colour=1,size = rel(0.8)),
              axis.title = element_text(size = rel(1)))+
        geom_rug(alpha=0.1)
      
      return(f)
    }
    multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
      
      # Make a list from the ... arguments and plotlist
      plots <- c(list(...), plotlist)
      
      numPlots = length(plots)
      
      # If layout is NULL, then use 'cols' to determine layout
      if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
      }
      
      if (numPlots==1) {
        print(plots[[1]])
        
      } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
          # Get the i,j matrix positions of the regions that contain this subplot
          matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
          
          print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                          layout.pos.col = matchidx$col))
        }
      }
    }
    g={}
    for (i in seq(along=yy)) g[[i]] <- ggs_density(yy[[i]],title=preya.names[i])
    
    if(sava==T) pdf(paste(save,"_densities.pdf",sep=''))
    multiplot(plotlist=g,cols=2)
    if(sava==T) dev.off()
    
    par(ask=T)
    
    ggs_double_violin <- function(D,num_cat,between_cat,within_cat,x.title,y.title,trans=NULL,breaks=NULL,ylims=NA,scallab=NULL){
      
      require(ggplot2)
      
      cbPalette <-c( "#ADD8E6","#87CEEB","#6495ED",  "#4169E1")
      num <- which(names(D)== substitute(num_cat))
      bet_cat <- which(names(D)== substitute(between_cat))
      with_cat <- which(names(D)== substitute(within_cat))
      
      if (is.null(scallab)) scallab = within_cat
      
      D <- data.frame(num = D[,num],bet_cat= D[,bet_cat],with_cat= D[,with_cat])
      
      if (!is.null(breaks)){
        m <- ggplot()+ geom_violin(data=D,aes(x=bet_cat,y=num),col=NA,fill=NA)+
          geom_rect(aes(xmin=0.5,xmax=1.5,ymin=0,ymax=1000),col=NA,fill='grey',alpha=0.8)+
          geom_rect(aes(xmin=2.5,xmax=3.5,ymin=0,ymax=1000),col=NA,fill='grey',alpha=0.8)+
          geom_violin(data=D,aes(x=bet_cat,y=num,fill=with_cat,col=with_cat))+
          scale_y_continuous(trans=paste(trans),breaks=breaks,limits=ylims)+
          ylab(y.title)+
          xlab(paste(x.title))+
          scale_fill_manual(substitute(scallab),values=cbPalette)+
          scale_colour_manual(substitute(scallab),values=cbPalette)+                          
          theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_text(vjust = -0.3),
                axis.title.y = element_text(vjust = 0.3),
                legend.title= element_text(size=14),
                legend.text = element_text(size=14),
                axis.text = element_text(colour='black',size=14),
                axis.title=element_text(colour='black',size=16)
                
          )
      } else {
        m <- ggplot()+ geom_violin(data=D,aes(x=bet_cat,y=num),col=NA,fill=NA)+
          geom_rect(aes(xmin=0.5,xmax=1.5,ymin=-Inf,ymax=Inf),col=NA,fill='grey',alpha=0.8)+
          geom_rect(aes(xmin=2.5,xmax=3.5,ymin=-Inf,ymax=Inf),col=NA,fill='grey',alpha=0.8)+
          geom_violin(data=D,aes(x=bet_cat,y=num,fill=with_cat,col=with_cat))+
          ylab(y.title)+
          xlab(paste(x.title))+
          scale_fill_manual(substitute(scallab),values=cbPalette)+
          scale_colour_manual(substitute(scallab),values=cbPalette)+
          theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_text(vjust = -0.3),
                axis.title.y = element_text(vjust = 0.3),
                legend.title= element_text(size=14),
                legend.text = element_text(size=14),
                axis.text = element_text(colour='black',size=14),
                axis.title=element_text(colour='black',size=16)
          )
      }
      
      
      
      m
    }
    
    
    names(outs)[2:(n.preys+1)] <- preya.names
    xx<-reshape::melt.data.frame(outs)
    
    names(xx) <- c('Method','Prey','draws')
    
    if(sava==T) pdf(paste(save,"_violin.pdf",sep=''))
    g <- ggs_double_violin(xx,draws,'Prey','Method','Prey','Density',trans=NULL,breaks=NULL,ylims=NA,scallab='Method')
    print(g)
    if(sava==T) dev.off()
    
  }
  
  par(ask=F)
}

#' Revert MCMC output class to coda's mcmc.list and use plot.mcmc to display chains. To be run with output from \code{\link{run_MCMC}}
#' @param x MCMC output from \code{\link{run_MCMC}}, containing diet proportion MCMC chains
#' @references Neubauer.P. and Jensen, O.P. (in prep) 
#' @author Philipp Neubauer
#' @seealso \code{\link{run_MCMC}},\code{\link{diags}}
#' @export
MCMCplot <- function(x){
  
  # check gui
  if(exists('GUI',envir=.GlobalEnv)) {
    GUI <- get('GUI',envir=.GlobalEnv)
  } else {
    GUI <- F
  }
  
  y  <- vector("list", x$nChains)
  for (i in 1:x$nChains) y[[i]] <- x[[i]]
  class(y) <- 'mcmc.list'
  
  if (!is.function(options()$device) & GUI==T){
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
   
    }
  }
  
  plot(y)
  
}