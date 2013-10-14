ggs_density <- function(D, pos=0) {
  # Manage subsetting a family of parameters
  require(ggplot2)
  require(reshape2)
  
  dims = dim(D)
#   if (is.null(dims)){
    D<-as.data.frame(D)
# }
    
  D<-melt(D,variable.name='Chain')
  names(D)[1] <- 'Chain'
  

  # Plot
 
    f <- ggplot(aes(x=value, linetype=as.factor(Chain)),data=D)
 
  f <- f + geom_density(alpha=0.3,adjust=2,data=D) + 
    geom_vline(aes(xintercept=pos),linetype=1:length(pos),data=as.data.frame(pos),colour=1)+
    geom_vline(aes(xintercept=0),colour=1)+
    geom_hline(aes(yintercept=0),colour=1)+
    coord_cartesian(xlim=c(0,max(D$value)))+
    theme(panel.grid.major=element_line(colour = NA),
  panel.grid.minor=element_line(colour = NA),
  panel.background=element_rect(fill=NA),
          legend.position = "none",
          axis.ticks = element_line(colour=1),
          axis.text = element_text(colour=1,size = rel(1.2)),
          axis.title = element_text(size = rel(1.5)))+
     geom_rug(alpha=0.1)
  
  return(f)
}
