#helper_functs


DD <- function(x){
  preya.names <- x$prey.names
  outs={}
  for (k in 1:x$nChains){
    outp <- matrix(unlist(x[[k]]),ncol=ncol(x[[1]]),byrow=F)
    outs <- rbind(outs,outp)
  }
  outs <- as.data.frame(outs)
  colnames(outs) <- preya.names
  D <-reshape::melt.data.frame(outs)
}


ggs_density <- function(D, pos=0,gp,title) {
  # Manage subsetting a family of parameters
  
  dims = dim(D)
  D<-as.data.frame(D)
  f <- ggplot(aes(x=value, fill=as.factor(variable), col=as.factor(variable)),data=D)+
    scale_colour_manual(values=c("#E69F00", "#56B4E9", "#009E73","#6A3D9A"))+
    scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73","#6A3D9A"))+
    geom_density(alpha=0.3,adjust=1.5,data=D) +
    geom_vline(aes(xintercept=pos),col=c("#E69F00", "#56B4E9", "#009E73","#6A3D9A"),alpha=1,data=as.data.frame(pos),linewidth=3) +
    geom_vline(aes(xintercept=value),col=rep(c("#6A3D9A"),each=5),alpha=1,data=as.data.frame(melt(gp)[16:20,]),linewidth=3,linetype=2) +
    geom_vline(aes(xintercept=0),colour=1)+
    geom_hline(aes(yintercept=0),colour=1)+      
    scale_x_continuous(limits=c(0,1))+
    xlab('Proportion')+
    ylab('Posterior Density')+ 
    theme(panel.grid.major=element_line(colour = NA),
          panel.grid.minor=element_line(colour = NA),
          panel.background=element_rect(fill=NA),
          #legend.position = "none",
          legend.title=element_blank(),
          axis.ticks = element_line(colour=1),
          axis.text = element_text(colour=1,size = rel(0.8)),
          axis.title = element_text(size = rel(1)))
  #geom_rug(alpha=0.1)
  
  f
}


  bad_prop <- read.csv('even_bad_props.csv',h=F)[,2:5]
good_prop <- read.csv('even_good_props.csv',h=F)[,2:5]
uneven_prop <- read.csv('uneven_good_props.csv',h=F)[,2:5]

D = DD(even.good.mcmc)
ggs_density(D, pos = colMeans(good_prop),good_prop,title)

D = DD(uneven.good.mcmc)
ggs_density(D, pos = colMeans(uneven_prop),uneven_prop,title)

D = DD(even.bad.mcmc)
ggs_density(D, pos = colMeans(bad_prop),bad_prop,title)
@

