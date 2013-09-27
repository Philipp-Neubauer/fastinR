plot.pop_props <- function(MCMCdatas){
    preya.names <- unique(guiGetSafe('datas')$prey.ix)
    outs <- MCMCdatas$MCMC
    plot(outs[,1],t='l',col=1,xlab='MCMC iteration',ylab='Proportion',ylim=c(0,1))
    for (i in 1:ncol(outs)) lines(outs[,i],col=i)

    par(ask=T)
    plot(outs,main='Correlation of proportion estimates')

    require(lattice)

    mp=reshape::melt(outs)

    rpp = bwplot(variable~value,data=mp,xlim=c(0,1),xlab=list(label=expression("Posterior distribution of diet proportions"),cex=1)
        ,panel = function(x, y) { 
            grid::grid.segments(1,0,0,0)
            xlist <- split(x, factor(y))
            for (i in seq(along=xlist))
                denstrip::panel.denstrip(x=xlist[[i]], at=i,colmax=(1),mticks=mean(xlist[[i]]))
        },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T,labels=preya.names)))
    print(rpp)
    
    trellis.focus("panel", 1, 1)
    panel.abline(h=seq(1.5,ncol(outs),1),col=1,lty=2)
    trellis.unfocus()

}
plot.ind_props <- function(MCMCdatas){
    datas <-  guiGetSafe('datas')
    preya.names <- unique(datas$prey.ix)
    outs <- MCMCdatas$MCMC
    popix <- grep('pop',colnames(MCMCdatas$MCMC))
 
    plot(outs[,popix[1]],t='l',col=1,xlab='MCMC iteration',ylab='Proportion',ylim=c(0,1))
    for (i in 1:ncol(outs[,popix])) lines(outs[,popix[i]],col=i)

    par(ask=T)
    plot(outs[,popix],main='Correlation of proportion estimates')

    require(lattice)

    # draw popualtion posteriors
    mp=reshape::melt.data.frame(outs[,popix])

    rpp = bwplot(variable~value,data=mp,xlim=c(0,1),xlab=list(label=expression("Posterior distribution of population diet proportions"),cex=1)
        ,panel = function(x, y) { 
            grid::grid.segments(1,0,0,0)
            xlist <- split(x, factor(y))
            for (i in seq(along=xlist))
                denstrip::panel.denstrip(x=xlist[[i]], at=i,colmax=(1),mticks=mean(xlist[[i]]))
        },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T,labels=preya.names)))
    print(rpp)
    par(ask=F)
    trellis.focus("panel", 1, 1)
    panel.abline(h=seq(1.5,ncol(outs),1),col=1,lty=2)
    trellis.unfocus()

     # draw individual posteriors
     popix <- grep('pop',colnames(MCMCdatas$MCMC))
     indix <- (max(popix)+1):ncol(MCMCdatas$MCMC)
     orders <- order(apply(MCMCdatas$MCMC[,indix[1:datas$n.preds]],2,median))
     indix <- rep(orders,each=3)+c(0,datas$n.preds,datas$n.preds*2)+max(popix)
     indix2 <- order(names(outs[,indix]))
     mp=reshape::melt.data.frame(outs[,indix])
par(ask=T)
    rpp = bwplot(variable~value,data=mp,xlim=c(0,1),xlab=list(label=expression("Posterior distribution of individual diet proportions"),cex=1),auto.key = list(text=preya.names, points = F,
               rectangles = F,
               lines = FALSE,
               col=1:3,alpha=1),
        ,panel = function(x, y) { 
            grid::grid.segments(1,0,0,0)
            xlist <- split(x, factor(y))
            k=0
            for (i in seq(along=xlist)){
                k=k+1
                denstrip::panel.denstrip(x=xlist[[i]], at=i,colmax=((i-1)%%datas$n.preys)+1,mticks=mean(xlist[[i]]))
            }
        },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T,labels=sprintf(paste0("Predator %d"),1:datas$n.preds),at=seq(datas$n.preys/2+0.5,datas$n.preys*datas$n.preds,datas$n.preys))))
    print(rpp)
    par(ask=F)
    trellis.focus("panel", 1, 1)
    panel.abline(h=seq(datas$n.preys+0.5,datas$n.preys*datas$n.preds,datas$n.preys),col=1,lty=2)
    trellis.unfocus()
    
}
plot.cov_props <- function(MCMCdatas){
    datas <-  guiGetSafe('datas')
    preya.names <- unique(datas$prey.ix)

    Covs <- guiGetSafe('Covs') 
    cidx <- apply(Covs,2,function(x){any(x!=0 & x!=1)})
                                        #number of groups and covariates
    nGr <- sum(cidx==F)
    Gridx <- which(cidx==F)

    outs <- MCMCdatas$MCMC
    popix <- grep('pop',colnames(MCMCdatas$MCMC))

    betaix <- grep('beta',colnames(MCMCdatas$MCMC))
    
    plot(outs[, betaix[1]],t='l',col=1,xlab='MCMC iteration',ylab='beta',ylim=range(outs[, betaix]),main='')
    for (i in 1:ncol(outs[, betaix])) lines(outs[, betaix[i]],col=i)

    par(ask=T)
    this.eff <- outs[,popix[1:(nGr*datas$n.preys)]]
    k=1
    for (n in 1:nGr) {colnames(this.eff)[k:(k+datas$n.preys-1)] <- sprintf(paste('Group',n,'/ %s'),preya.names); k=k+datas$n.preys}
    plot( this.eff,main='Correlation of proportion estimates')

    require(lattice)

    # draw popualtion posteriors
    mp=reshape::melt.data.frame(this.eff)

    rpp = bwplot(variable~value,data=mp,xlim=c(0,1),xlab=list(label=expression("Posterior distribution of population diet proportions"),cex=1)
        ,panel = function(x, y) { 
            grid::grid.segments(1,0,0,0)
            xlist <- split(x, factor(y))
            for (i in seq(along=xlist))
                denstrip::panel.denstrip(x=xlist[[i]], at=i,colmax=(i%%3)+1,mticks=mean(xlist[[i]]))
        },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T)))
    print(rpp)
    par(ask=F)
    trellis.focus("panel", 1, 1)
    panel.abline(h=seq(datas$n.preys+0.5,ncol(this.eff),datas$n.preys),col=1,lty=2)
    trellis.unfocus()

     # draw individual posteriors
     popix <- grep('pop',colnames(MCMCdatas$MCMC))
     indix <- (max(popix)+1):ncol(MCMCdatas$MCMC)
     orders <- order(apply(MCMCdatas$MCMC[,indix[1:datas$n.preds]],2,median))
     indix <- rep(orders,each=3)+c(0,datas$n.preds,datas$n.preds*2)+max(popix)
     indix2 <- order(names(outs[,indix]))
     mp=reshape::melt.data.frame(outs[,indix])
    par(ask=T)
    rpp = bwplot(variable~value,data=mp,xlim=c(0,1),xlab=list(label=expression("Posterior distribution of individual diet proportions"),cex=1),auto.key = list(text=preya.names, points = F,
               rectangles = F,
               lines = FALSE,
               col=1:3,alpha=1),
        ,panel = function(x, y) { 
            grid::grid.segments(1,0,0,0)
            xlist <- split(x, factor(y))
            k=0
            for (i in seq(along=xlist)){
                k=k+1
                denstrip::panel.denstrip(x=xlist[[i]], at=i,colmax=((i-1)%%datas$n.preys)+1,mticks=mean(xlist[[i]]))
            }
        },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T,labels=sprintf(paste0("Predator %d"),1:datas$n.preds),at=seq(datas$n.preys/2+0.5,datas$n.preys*datas$n.preds,datas$n.preys))))
    print(rpp)
    par(ask=F)
    trellis.focus("panel", 1, 1)
    panel.abline(h=seq(datas$n.preys+0.5,datas$n.preys*datas$n.preds,datas$n.preys),col=1,lty=2)
    trellis.unfocus()
}
