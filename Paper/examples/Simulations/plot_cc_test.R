load('./examples/Simulations/cc.test.Rdata')
cc.test <- matrix(unlist(cc.test),ncol=3,byrow=T)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,upper, x, lower, angle=90, code=3, length=length, ...)
}


### first plot densities

par(mar=c(5, 4, 4, 4))
h=density(rnorm(10000,1,0.5),adjust=2)
plot((h$y/(2*max(h$y)))+5,h$x,t='l',lty=2,xlim=c(0,6),ylim=c(0,3),axes=F,xlab='',ylab='')
title(ylab=expression(kappa))
h=density(rnorm(10000,1,0.4),adjust=2)
lines((h$y/(2*max(h$y)))+4,h$x,lty=2)
h=density(rnorm(10000,1,0.3),adjust=2)
lines((h$y/(2*max(h$y)))+3,h$x,lty=2)
h=density(rnorm(10000,1,0.2),adjust=2)
lines((h$y/(2*max(h$y)))+2,h$x,lty=2)
h=density(rnorm(10000,1,0.1),adjust=2)
lines((h$y/(2*max(h$y)))+1,h$x,lty=2)
axis(2,0:3,cex=1.2)

par(new=T,mar=c(5, 4, 4, 4))

cmeans <- aggregate(cc.test[,2:3],list(cc.test[,1]),mean)
plot(seq(0,0.5,0.1)-0.005,cmeans[,2],axes=F,ylim=c(-2,4),xlim=c(0,0.6),pch=16,cex=1.2,xlab='',ylab='')
axis(1,at=seq(0,0.5,0.1),cex=1.2)
title(xlab=expression(sigma[kappa]))
axis(4,cex=1.2)
mtext(expression(Delta[pi]),4,line=2.5)
mtext(expression(kappa),2,line=2.5)

qs <-aggregate(cc.test[,2:3],list(cc.test[,1]),function(x){quantile(x,c(0.025,0.975),na.rm=T)})
error.bar(seq(0,0.5,0.1)-0.005,cmeans[,2],qs[,2][,2],qs[,2][,1],lwd=1.2)

points(seq(0.,0.5,0.1)+0.005,cmeans[,3],axes=NULL,pch=1,cex=1.2)
#qs <- apply((cc.test[,2]),1,function(x){quantile(x,c(0.05,0.95),na.rm=T)})
error.bar(seq(0,0.5,0.1)+0.005,cmeans[,3],qs[,3][,2],qs[,3][,1],lwd=1.2)

abline(h=0,lty=3)
