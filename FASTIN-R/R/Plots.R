plot.pop_props <- function(x,...){
    
    sava <- menu(title='save plots?',choices = c('yes','no'),graphics=T)
   
    preya.names <- unique(guiGetSafe('datas')$prey.ix)

     outs={}
        for (k in 1:length(x)){
            outp <- matrix(unlist(x[[k]]),ncol=ncol(x[[1]]),byrow=F)
            outs <- rbind(outs,outp)
        }
    colnames(outs) <- preya.names
    outs <- as.data.frame(outs)
    
    if (sava==1) pdf("MCMC_correlations.pdf")
    plot(outs,main='Correlation of proportion estimates')
    if (sava==1) dev.off()
   

    mp=reshape::melt.data.frame(outs)
    if (sava==1) pdf("Post_pop_proportions.pdf")
    par(ask=T)
    rpp = bwplot(variable~value,data=mp,xlim=c(0,1),xlab=list(label=expression("Posterior distribution of diet proportions"),cex=1)
        ,panel = function(x, y) { 
            grid::grid.segments(1,0,0,0)
            xlist <- split(x, factor(y))
            for (i in seq(along=xlist))
                denstrip::panel.denstrip(x=xlist[[i]], at=i,colmax=(1),mticks=mean(xlist[[i]]))
        },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T,labels=preya.names)))
    print(rpp)
    par(ask=F)
    trellis.focus("panel", 1, 1,highlight=F)
    panel.abline(h=seq(1.5,ncol(outs),1),col=1,lty=2)
    trellis.unfocus()
    if (sava==1) dev.off()
}
plot.ind_props <- function(x,...){
    
    sava <- menu(title='save plots?',choices = c('yes','no'),graphics=T)
    datas <-  guiGetSafe('datas')
    preya.names <- unique(datas$prey.ix)

    outs={}
        for (k in 1:length(x)){
            outp <- matrix(unlist(x[[k]]),ncol=ncol(x[[1]]),byrow=F)
            outs <- rbind(outs,outp)
        }
    outs <- as.data.frame(outs)
    
    colnames(outs) <- colnames(x[[1]])
    popix <- grep('pop',colnames(outs))
    
    if (sava==1) pdf("MCMC_correlations.pdf")
    par(ask=T)
    plot(outs[,popix],main='Correlation of proportion estimates')
    if (sava==1) dev.off()
   
    # draw popualtion posteriors
    mp=reshape::melt.data.frame(outs[,popix])
    if (sava==1) pdf("Post_pop_proportions.pdf")
    rpp = bwplot(variable~value,data=mp,xlim=c(0,1),xlab=list(label=expression("Posterior distribution of population diet proportions"),cex=1)
        ,panel = function(x, y) { 
            grid::grid.segments(1,0,0,0)
            xlist <- split(x, factor(y))
            for (i in seq(along=xlist))
                denstrip::panel.denstrip(x=xlist[[i]], at=i,colmax=(1),mticks=mean(xlist[[i]]))
        },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T,labels=preya.names)))
    print(rpp)
    par(ask=F)
    trellis.focus("panel", 1, 1,highlight=F)
    panel.abline(h=seq(1.5,ncol(outs),1),col=1,lty=2)
    trellis.unfocus()
    if (sava==1) dev.off()
    
     # draw individual posteriors
    
     indix <- (max(popix)+1):ncol(outs)
     orders <- order(apply(outs[,indix[1:datas$n.preds]],2,median))
     this.prey <- which.max(tapply(apply(outs[,indix],2,median),rep(1:datas$n.preys,each=datas$n.preds),median))
     orders <- order(apply(outs[,indix[(1:datas$n.preds)+(this.prey-1)*datas$n.preds]],2,median))
     indix <- rep(orders,each=datas$n.preys)+c(0,datas$n.preds*1:(datas$n.preys-1))+max(popix)
     indix2 <- order(names(outs[,indix]))
     mp=reshape::melt.data.frame(outs[,indix])
par(ask=T)
    if (sava==1) pdf("Post_ind_proportions.pdf")
    rpp = bwplot(variable~value,data=mp,xlim=c(0,1),xlab=list(label=expression("Posterior distribution of individual diet proportions"),cex=1),auto.key = list(text=preya.names, points = F,
               rectangles = F,
               lines = FALSE,
               col=1:datas$n.preds,alpha=1,space='right'),
        ,panel = function(x, y) { 
            grid::grid.segments(1,0,0,0)
            xlist <- split(x, factor(y))           
            for (i in seq(along=xlist)){
                denstrip::panel.denstrip(x=xlist[[i]], at=i,colmin="transparent",colmax=((i-1)%%datas$n.preys)+1,mticks=mean(xlist[[i]]),width=1/datas$n.preds*datas$n.preys)
            }
        },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T,labels=sprintf(paste0("Predator %d"),orders),at=seq(datas$n.preys/2+0.5,datas$n.preys*datas$n.preds,datas$n.preys))))
    print(rpp)
    par(ask=F)
    trellis.focus("panel", 1, 1)
    panel.abline(h=seq(datas$n.preys+0.5,datas$n.preys*datas$n.preds,datas$n.preys),col=1,lty=2)
    trellis.unfocus()
    if (sava==1) dev.off()
}
plot.cov_props <- function(x,...){
   
    sava <- menu(title='save plots?',choices = c('yes','no'),graphics=T)
    datas <-  guiGetSafe('datas')
    preya.names <- unique(datas$prey.ix)

     outs={}
        for (k in 1:length(x)){
            outp <- matrix(unlist(x[[k]]),ncol=ncol(x[[1]]),byrow=F)
            outs <- rbind(outs,outp)
        }
    outs <- as.data.frame(outs)
    
    colnames(outs) <- colnames(x[[1]])
    
    Covs <- guiGetSafe('Covs') 
    cidx <- apply(Covs,2,function(x){any(x!=0 & x!=1)})
                                        #number of groups and covariates
    nGr <- sum(cidx==F)
    Gridx <- which(cidx==F)

    outs <- MCMCdatas$MCMC
    popix <- grep('pop',colnames(outs))

    betaix <- grep('beta',colnames(outs))

    par(ask=T)
    this.eff <- outs[,popix[1:(nGr*datas$n.preys)]]
    k=1
    if (sava==1) pdf("MCMC_correlations.pdf")
    for (n in 1:nGr) {colnames(this.eff)[k:(k+datas$n.preys-1)] <- sprintf(paste('Group',n,'/ %s'),preya.names); k=k+datas$n.preys}
    plot( this.eff,main='Correlation of proportion estimates')
    if (sava==1) dev.off()
   

    # draw popualtion posteriors
    mp=reshape::melt.data.frame(this.eff)
    if (sava==1) pdf("Post_pop_proportions.pdf")
    rpp = bwplot(variable~value,data=mp,xlim=c(0,1),xlab=list(label=expression("Posterior distribution of population diet proportions"),cex=1)
        ,panel = function(x, y) { 
            grid::grid.segments(1,0,0,0)
            xlist <- split(x, factor(y))
            for (i in seq(along=xlist))
                denstrip::panel.denstrip(x=xlist[[i]], at=i,colmax=((i-1)%%datas$n.preys)+1,mticks=mean(xlist[[i]]))
        },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T)))
    print(rpp)
    par(ask=F)
    trellis.focus("panel", 1, 1)
    panel.abline(h=seq(datas$n.preys+0.5,ncol(this.eff),datas$n.preys),col=1,lty=2)
    trellis.unfocus()
    if (sava==1) dev.off()
    
     # draw individual posteriors
  
     indix <- (max(popix)+1):ncol(outs)
    # check which prey to ordewr by (highest proportion)
   this.prey <- which.max(tapply(apply(outs[,indix],2,median),rep(1:datas$n.preys,each=datas$n.preds),median))
     orders <- order(apply(outs[,indix[(1:datas$n.preds)+(this.prey-1)*datas$n.preds]],2,median))
     indix <- rep(orders,each=datas$n.preys)+c(0,datas$n.preds*1:(datas$n.preys-1))+max(popix)
     indix2 <- order(names(outs[,indix]))
     mp=reshape::melt.data.frame(outs[,indix])
    par(ask=T)
    if (sava==1) pdf("Post_ind_proportions.pdf")
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
                denstrip::panel.denstrip(x=xlist[[i]], at=i,colmin="transparent",colmax=((i-1)%%datas$n.preys)+1,mticks=mean(xlist[[i]]),width=1/datas$n.preds*datas$n.preys)
            }
        },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T,labels=sprintf(paste0("Predator %d"),orders),at=seq(datas$n.preys/2+0.5,datas$n.preys*datas$n.preds,datas$n.preys))))
    print(rpp)
    par(ask=F)
    trellis.focus("panel", 1, 1,highlight=F)
    panel.abline(h=seq(datas$n.preys+0.5,datas$n.preys*datas$n.preds,datas$n.preys),col=1,lty=2)
    trellis.unfocus()
    if (sava==1) dev.off()
}
