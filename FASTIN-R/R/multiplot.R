multiplot <- function(MCMCouts,density=F){
   
    sava <- menu(title='save plots?',choices = c('yes','no'),graphics=T)
    datas <-  guiGetSafe('datas')
    preya.names <- unique(datas$prey.ix)

     if(density==F){
    
    outs <- {}
    for (i in 1:length(MCMCouts)){

        x <- MCMCouts[[i]]
        outa={}
        for (k in 1:length(x)){
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
                denstrip::panel.denstrip(x=xlist[[i]], at=i,colmax=((i-1)%%datas$n.preys)+1,mticks=mean(xlist[[i]]))
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
        for (k in 1:length(x)){
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
  require(ggplot2)
  
  dims = dim(D)
#   if (is.null(dims)){
    D<-as.data.frame(D)
# }
   
  # Plot
 
    f <- ggplot(aes(x=value, col=as.factor(V1)),data=D)
 
  f <- f + geom_density(alpha=0.3,adjust=4,data=D) + 
    geom_vline(aes(xintercept=pos),linetype=1:length(pos),data=as.data.frame(pos),colour=1)+
    geom_vline(aes(xintercept=0),colour=1)+
    geom_hline(aes(yintercept=0),colour=1)+
    coord_cartesian(xlim=c(0,1))+
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
  require(grid)

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

