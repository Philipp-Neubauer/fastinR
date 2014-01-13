#' @name DietProportionPlot
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
plot.pop_props <- function(x,save="FASTIN_MCMC_",density=T,...){
  
  if(save!=F){sava <- menu(title='save plots?',choices = c('yes','no'),graphics=T)}else{sava=0}
  
  preya.names <- x$prey.names
  
  outs={}
  for (k in 1:x$nChains){
    outp <- matrix(unlist(x[[k]]),ncol=ncol(x[[1]]),byrow=F)
    outs <- rbind(outs,outp)
  }
  colnames(outs) <- preya.names
  outs <- as.data.frame(outs)
  
  if (sava==1) {
    pdf(paste(save,"correlations.pdf",sep=''))
  } else {externalDevice<-FALSE
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
    
    f <- f + geom_density(alpha=0.3,adjust=1.2,data=D) + 
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
plot.ind_props <- function(x,save="FASTIN_MCMC_",density=T,...){
  
  if(save!=F){sava <- menu(title='save plots?',choices = c('yes','no'),graphics=T)}else{sava=0}
  
  preya.names <- x$prey.names
  
  outs={}
  for (k in 1:x$nChains){
    outp <- matrix(unlist(x[[k]]),ncol=ncol(x[[1]]),byrow=F)
    outs <- rbind(outs,outp)
  }
  outs <- as.data.frame(outs)
  
  colnames(outs) <- colnames(x[[1]])
  popix <- grep('pop',colnames(outs))
  
  if (sava==1)  {
    pdf(paste(save,"correlations.pdf",sep=''))
  } else {externalDevice<-FALSE
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
      
      f <- f + geom_density(alpha=0.3,adjust=1.2,data=D) + 
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
  orders <- order(apply(outs[,indix[1:datas$n.preds]],2,median))
  this.prey <- which.max(tapply(apply(outs[,indix],2,median),rep(1:datas$n.preys,each=datas$n.preds),median))
  orders <- order(apply(outs[,indix[(1:datas$n.preds)+(this.prey-1)*datas$n.preds]],2,median))
  indix <- rep(orders,each=datas$n.preys)+c(0,datas$n.preds*1:(datas$n.preys-1))+max(popix)
  indix2 <- order(names(outs[,indix]))
  
  if(density==F){
    mp=reshape::melt.data.frame(outs[,indix])
    
    if (sava==1){pdf(paste(save,"ind_proportions.pdf",sep=''))}else{par(ask=T)}
    rpp = bwplot(variable~value,data=mp,xlim=c(0,1),xlab=list(label=expression("Posterior distribution of individual diet proportions"),cex=1),auto.key = list(text=preya.names, points = F,
                                                                                                                                                               rectangles = F,
                                                                                                                                                               lines = FALSE,
                                                                                                                                                               col=1:datas$n.preds,alpha=1,space='right'),
                 ,panel = function(x, y) { 
                   grid::grid.segments(1,0,0,0)
                   xlist <- split(x, factor(y))           
                   for (i in seq(along=xlist)){
                     denstrip::panel.denstrip(x=xlist[[i]], at=i,colmin="transparent",colmax=((i-1)%%datas$n.preys)+1,ticks=quantile(xlist[[i]],c(0.025,0.975)),mticks=quantile(xlist[[i]],0.5),tlen=1.2,mlen=1.2,nwd=1.5,width=1/datas$n.preds*datas$n.preys)
                   }
                 },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T,labels=sprintf(paste0("Predator %d"),orders),at=seq(datas$n.preys/2+0.5,datas$n.preys*datas$n.preds,datas$n.preys))))
    print(rpp)
    par(ask=F)
    trellis.focus("panel", 1, 1)
    panel.abline(h=seq(datas$n.preys+0.5,datas$n.preys*datas$n.preds,datas$n.preys),col=1,lty=2)
    trellis.unfocus()
    if (sava==1) dev.off()
    
  } else {
    
    colnames(outs) <-NA
    outr <- data.frame()
    a=1
    for (i in 1:datas$n.preds){
      outr <- rbind(outr,cbind(rep(i,nrow(outs)),outs[,indix[a:(a+datas$n.preys-1)]]))
      a = a+datas$n.preys
    }
    
    colnames(outr)[2:(datas$n.preys+1)] <- preya.names    
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
      
      f <- f + geom_density(alpha=0.3,adjust=1.2,data=D,trim=F)+ 
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
plot.cov_props <- function(x,save="FASTIN_MCMC_",density=T,...){
  
  if(save!=F){sava <- menu(title='save plots?',choices = c('yes','no'),graphics=T)}else{sava=0}
    
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
  
  betaix <- grep('beta',colnames(outs))
  
  this.eff <- outs[,popix[1:(nGr*datas$n.preys)]]
  k=1
  if (sava==1)  {
    pdf(paste(save,"correlations.pdf",sep=''))
  } else {externalDevice<-FALSE
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
  }
  for (n in 1:nGr) {colnames(this.eff)[k:(k+datas$n.preys-1)] <- sprintf(paste('Group',n,'/ %s'),preya.names); k=k+datas$n.preys}
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
                     denstrip::panel.denstrip(x=xlist[[i]], at=i,colmax=((i-1)%%datas$n.preys)+1,ticks=quantile(xlist[[i]],c(0.025,0.975)),mticks=quantile(xlist[[i]],0.5),tlen=1.2,mlen=1.2,nwd=1.5)
                 },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T)))
    print(rpp)
    par(ask=F)
    trellis.focus("panel", 1, 1)
    panel.abline(h=seq(datas$n.preys+0.5,ncol(this.eff),datas$n.preys),col=1,lty=2)
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
      
      f <- f + geom_density(alpha=0.3,adjust=1.2,data=D) + 
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
  # check which prey to ordewr by (highest proportion)
  this.prey <- which.max(tapply(apply(outs[,indix],2,median),rep(1:datas$n.preys,each=datas$n.preds),median))
  orders <- order(apply(outs[,indix[(1:datas$n.preds)+(this.prey-1)*datas$n.preds]],2,median))
  indix <- rep(orders,each=datas$n.preys)+c(0,datas$n.preds*1:(datas$n.preys-1))+max(popix)
  indix2 <- order(names(outs[,indix]))
  
  if(density==F){
    mp=reshape::melt.data.frame(outs[,indix])
    if (sava==1){pdf(paste(save,"ind_proportions.pdf",sep=''))}else{par(ask=T)}
    rpp = bwplot(variable~value,data=mp,xlim=c(0,1),xlab=list(label=expression("Posterior distribution of individual diet proportions"),cex=1),auto.key = list(text=preya.names, points = F,
                                                                                                                                                               rectangles = F,
                                                                                                                                                               lines = FALSE,
                                                                                                                                                               col=1:datas$n.preys,alpha=1,space='right'),
                 ,panel = function(x, y) { 
                   grid::grid.segments(1,0,0,0)
                   xlist <- split(x, factor(y))
                   k=0
                   for (i in seq(along=xlist)){
                     k=k+1
                     denstrip::panel.denstrip(x=xlist[[i]], at=i,colmin="transparent",colmax=((i-1)%%datas$n.preys)+1,ticks=quantile(xlist[[i]],c(0.025,0.975)),mticks=quantile(xlist[[i]],0.5),tlen=1.2,mlen=1.2,nwd=1.5,width=1/datas$n.preds*datas$n.preys)
                   }
                 },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T,labels=sprintf(paste0("Predator %d"),orders),at=seq(datas$n.preys/2+0.5,datas$n.preys*datas$n.preds,datas$n.preys))))
    print(rpp)
    par(ask=F)
    trellis.focus("panel", 1, 1,highlight=F)
    panel.abline(h=seq(datas$n.preys+0.5,datas$n.preys*datas$n.preds,datas$n.preys),col=1,lty=2)
    trellis.unfocus()
    if (sava==1) dev.off()
    
  } else {
    
    colnames(outs) <-NA
    outr <- data.frame()
    a=1
    for (i in 1:datas$n.preds){
      outr <- rbind(outr,cbind(rep(i,nrow(outs)),outs[,indix[a:(a+datas$n.preys-1)]]))
      a = a+datas$n.preys
    }
    
    colnames(outr)[2:(datas$n.preys+1)] <- preya.names    
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
      
      f <- f + geom_density(alpha=0.3,adjust=1.2,data=D,trim=F)+ 
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
#' @details If plots are saved they are not drawn at the same time.
#' @references Neubauer.P. and Jensen, O.P. (in prep) 
#' @author Philipp Neubauer
#' @seealso \code{\link{run_MCMC}},\code{\link{diags}}
#' @export
multiplot <- function(MCMCouts,density=T){
  
  sava <- menu(title='save plots?',choices = c('yes','no'),graphics=T)
  
  if (sava==0) {
    externalDevice<-FALSE
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
  }
  
  preya.names <- unique(MCMCouts[[i]]$prey.ix)
    
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
    for (n in 1:length(MCMCouts)) {colnames(outs)[k:(k+datas$n.preys-1)] <- sprintf(paste(MCnames[n],'/ %s'),preya.names); k=k+datas$n.preys}
    
    # draw popualtion posteriors
    mp=reshape::melt.data.frame(outs)
    if (sava==1) pdf("Multiplot_proportions.pdf")
    rpp = bwplot(variable~value,data=mp,xlim=c(0,1),xlab=list(label=expression("Posterior distribution of population diet proportions"),cex=1)
                 ,panel = function(x, y) { 
                   grid::grid.segments(1,0,0,0)
                   xlist <- split(x, factor(y))
                   for (i in seq(along=xlist))
                     denstrip::panel.denstrip(x=xlist[[i]], at=i,colmax=((i-1)%%datas$n.preys)+1,ticks=quantile(xlist[[i]],c(0.025,0.975)),mticks=quantile(xlist[[i]],0.5),tlen=1.2,mlen=1.2,nwd=1.5)
                 },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T)))
    print(rpp)
    
    trellis.focus("panel", 1, 1,highlight=F)
    panel.abline(h=seq(datas$n.preys+0.5,ncol(outs),datas$n.preys),col=1,lty=2)
    trellis.unfocus()
    if (sava==1) dev.off()
  } else {
    
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
      
      f <- f + geom_density(alpha=0.3,adjust=1.2,data=D) + 
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
    
    if(sava==T) pdf('multiplot_densities.pdf')
    multiplot(plotlist=g,cols=2)
    if(sava==T) dev.off()
  }}

#' Revert MCMC output class to coda's mcmc.list and use plot.mcmc to display chains. To be run with output from \code{\link{run_MCMC}}
#' @param x MCMC output from \code{\link{run_MCMC}}, containing diet proportion MCMC chains
#' @references Neubauer.P. and Jensen, O.P. (in prep) 
#' @author Philipp Neubauer
#' @seealso \code{\link{run_MCMC}},\code{\link{diags}}
#' @export
MCMCplot <- function(x){
  
  y  <- vector("list", x$nChains)
  for (i in 1:x$nChains) y[[i]] <- x[[i]]
  class(y) <- 'mcmc.list'
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
  
  plot(y)
  
}