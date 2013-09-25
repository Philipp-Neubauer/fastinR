plot.pop_props <- function(MCMCdatas){
    preya.names <- guiGetSafe("prey.names")
    outs <- MCMCdatas$MCMC
    plot(outs[,1],t='l',col=1,xlab='MCMC iteration',ylab='Proportion')
    for (i in 1:ncol(outs)) lines(outs[,i],col=i)

    par(ask=T)
    plot(outs,main='Correlation of proportion estimates')
}
plot.ind_props <- function(MCMCdatas){}
plot.cov_props <- function(MCMCdatas){}
