# prelims ----
require(FASTIN)

# simulations for CC var dependence --------------

n.fats=12
n.preys=3
nsamples=rep(50,n.preys)

# sep sets (roughly,probabilistically) the separation of prey items in multivariate (compositional) space
sep=3
n.preds =1

testa3 <- array(,c(5,100,2))

for(r in 1:100){
    ks=0;k=0.45;
  for (k in seq(0.05,0.45,0.1)){
  ks=ks+1
    As = matrix(rlnorm(n.preys*n.fats,rnorm(n.preys*n.fats,sep)),n.preys,n.fats,byrow=T)
    
    preys <- array(,c(nsamples[1],n.fats,n.preys))
    for (i in 1:n.preys){
      preys[,,i] <- rdirichlet(nsamples[i],As[i,])
    }
    
    preys.ix={};preya={}
    for (i in 1:n.preys){
      preys.ix <- c(preys.ix,rep(i,nsamples[i]))
      preya<- rbind(preya,preys[,,i])
    }
    
    #preyRDA <- rda(alr(preya) ~ as.factor(preys.ix))
    
    #randomly draw fat content
    fc_mean <- (sample(40,n.preys))/10+1
    fc_mean
    # set variance of fat content to 0.4 - this is loosely based on literature
    fc_tau <- rep(1/0.4,n.preys)
    
# draw mean conversion coeffs - cvar is the variance of mean conversion coeffs around 1; 0.2 seems reasonable  (c.f., Rosen & Tollit histogram)
cvar = k
mean_css=matrix(rlnorm(n.fats*n.preys,0,cvar),n.preys,n.fats) # rtnorm is a truncated normal, with lower truncation set to 0
mean_c=mean_css # just so I can modify mean_cs later without loosing the actual mean
var_coeffs = matrix((0.03)^2,n.preys,n.fats)  # set the variance to some constant - pretty close, at least for FA whose CC isn't too far from 1 (Iverson 2004 graphs)
tau_coeffs = 1/var_coeffs

# draw individual coeffs for each predator (e.g., may depend on individual diet status - fasting vs low fat diet versus hight fat diet) from a distribution about mean coeffs
#cs <- t(matrix(rlnorm(n.preys*n.fats,log(mean_c),sqrt(var_coeffs)),n.fats,n.preys))
#hist(cs)
# these are the proportion that each predator eats of each prey
# population 'mean', Dirichlet alpha notation, hard to interpret but computationally efficient.
# v(1,2) determines the variance of the alpha coeffs - the higher this variance, the lower the among predator variation in diet proportions.
v1 = 0.5; v2 = 0.7
Q = rlnorm(n.preys,rlnorm(n.preys,0,v1),v2)
# individual proportions
props=rdirichlet(n.preds,Q)
# show simulated proportions
props

# the mixture in matrix algebra (just a matrix product)
# prey mean is the geometric mean of each column
mprey <-  clo(t(apply(preys,3,function(x){apply(x,2,function(y){exp(mean(log(y)))})})))
preda <-  clo(props%*%(fc_mean*mprey*mean_c))
  
# normalize and transform
preds <- unclass(data.matrix(alr(preda)))

# transform prey stats
preym <-  alr(mprey)[,]

m.fats = (n.fats-1)
# get sums of squares for each prey matrix
R <- array(,c(m.fats,m.fats,n.preys))
for (i in 1:n.preys){
  R[,,i]=(cov(alr(preys[,,i]))*(nsamples[i]-1))
  #R[,,i] <- s$u %*% solve(diag(s$d)) %*% t(s$v)
}

# # set SS matrix for wishart prior on proportions
# p_tau=0.2
# 2*(1-pnorm(log(95)/2,0,sqrt(1/p_tau))) # chance that p1=95*p2
# S = diag(p_tau,n.preys)
# SS = diag(1,n.preys)
# set uninformative prior SS matrix for wishart prior alr transformed predator data
Rnot =diag(1,m.fats)

#zeros=rep(0,n.preys)


initials=list(list(
  fc = fc_mean,
  fracs=matrix(mean_c,n.preys,n.fats),
  ps = rep(1/n.preys,n.preys),
  #pmean=rep(0,n.preys),
  #A = Q,#rep(1/n.preys,n.preys),#
  #pprec=diag(1,n.preys), # precision on proportion covariance
  prey.means=preym,
  predprec = diag(1,m.fats),
  prey.precs = array(0.1,c(n.preys,m.fats,m.fats))
  
))

# use real mean_c first, then do again with CC set to one
for(tries in 1:2){
# look at how much the assumption of a constant c biases the results
if (tries == 2){
  mean_c=matrix(1,n.preys,n.fats) # set c to one for all CC
  var_coeffs = matrix(k^2,n.preys,n.fats)  # set the variance to some constant - pretty close, at least for FA whose CC isn't too far from 1 (Iverson 2004 graphs)
  tau_coeffs = 1/var_coeffs}

############# don't try this unless you chose to simualte relatively few FA (else it will take literally forever)

#bugsInits(initials,1,'FattyInits.R')

#bugsData(c('R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','nsamples','preds','preym'),'datas.R')
datas=list('R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','nsamples','preds','preym')

vars = c('prop')

# compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...

Fattysz<- BRugsFit('Fatty_dir_Single_prop.R', datas, inits=initials, numChains = 1, vars,
         nBurnin = 1000, nIter = 10000, nThin = 10, coda = T,
         DIC = F, working.directory = getwd(), digits = 4, seed=NULL,
         BRugsVerbose = getOption("BRugsVerbose"))


Fatty= as.data.frame((Fattysz)[[1]])

prop.ix <- grep('prop',names(Fatty))
#hist(Fattyz[,prop.ix[1]],20)
#plot(Fatty[,prop.ix[1]])
# 
# pred.ix <- grep('prey',names(Fatty))
# plot(Fatty[,pred.ix[1]],t='l')
# postpreys = matrix(colMeans(Fatty[,pred.ix]),3,n.fats,byrow=T)
# 
# g=1
# pps <- ggs_density(cbind(Fatty[,prop.ix[g]],Fattzy[,prop.ix[g]]),c(mean(Fatty[,prop.ix[g]]),mean(Fattzy[,prop.ix[g]]),props[g]))
# pps

postprops = matrix(colMeans(Fatty[,prop.ix]),n.preds,n.preys,byrow=T)
#postprops2 = matrix(apply(Fatty[,prop.ix],2,function(x){density(x,adjust=2)$x[density(x,adjust=2)$y==max(density(x,adjust=2)$y)]}),n.preds,n.preys,byrow=T)

#post.props - props

testa3[ks,r,tries] <-dist(alr(rbind(postprops,props)))

#### log normal model

mean_c=log(mean_css) # log scale mean
var_coeffs = matrix((0.03)^2,n.preys,n.fats)  # set the variance to some constant - pretty close, at least for FA whose CC isn't too far from 1 (Iverson 2004 graphs)
tau_coeffs = 1/var_coeffs

# look at how much the assumption of a constant c biases the results
if (tries == 2){
  mean_c=matrix(0,n.preys,n.fats) # set c to one for all CC
  var_coeffs = matrix(k,n.preys,n.fats)  # set the variance to some constant - pretty close, at least for FA whose CC isn't too far from 1 (Iverson 2004 graphs)
  tau_coeffs = 1/var_coeffs}

Fattysz<- BRugsFit('Fatty_dir_Single_prop_ln.R', datas, inits=initials, numChains = 1, vars,
                   nBurnin = 1000, nIter = 10000, nThin = 10, coda = T,
                   DIC = F, working.directory = getwd(), digits = 4, seed=NULL,
                   BRugsVerbose = getOption("BRugsVerbose"))


Fatty= as.data.frame((Fattysz)[[1]])

prop.ix <- grep('prop',names(Fatty))
#hist(Fattyz[,prop.ix[1]],20)
#plot(Fatty[,prop.ix[1]])
# 
# pred.ix <- grep('prey',names(Fatty))
# plot(Fatty[,pred.ix[1]],t='l')
# postpreys = matrix(colMeans(Fatty[,pred.ix]),3,n.fats,byrow=T)
# 
# g=1
#   pps <- ggs_density(cbind(Fatty[,prop.ix[1]],Fatty[,prop.ix[2]],Fatty[,prop.ix[3]]),c(mean(Fatty[,prop.ix[1]]),mean(Fatty[,prop.ix[2]]),mean(Fatty[,prop.ix[3]]),props))
#   pps

postprops = matrix(colMeans(Fatty[,prop.ix]),n.preds,n.preys,byrow=T)

testa3[ks,r,tries+2] <-dist(alr(rbind(postprops,props)))
}}}#this one for tries

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

### first plot densities

par(mar=c(5, 4, 4, 4))
h=density(rnorm(10000,1,0.45),adjust=2)
plot((h$y/(2*max(h$y)))+4,h$x,t='l',lty=2,xlim=c(0,5),ylim=c(0,3),axes=F,xlab='',ylab='')
title(ylab=expression(kappa))
h=density(rnorm(10000,1,0.35),adjust=2)
lines((h$y/(2*max(h$y)))+3,h$x,lty=2)
h=density(rnorm(10000,1,0.25),adjust=2)
lines((h$y/(2*max(h$y)))+2,h$x,lty=2)
h=density(rnorm(10000,1,0.15),adjust=2)
lines((h$y/(2*max(h$y)))+1,h$x,lty=2)
h=density(rnorm(10000,1,0.05),adjust=2)
lines((h$y/(2*max(h$y))),h$x,lty=2)
axis(2,0:3,cex=1.2)

par(new=T,mar=c(5, 4, 4, 4))

plot(seq(0.05,0.45,0.1)-0.015,rowMeans(testa3[,,1],na.rm=T),axes=N,ylim=c(-2,4),xlim=c(0.05,0.55),pch=16,cex=1.2)
axis(1,at=seq(0.05,0.45,0.1),cex=1.2)
title(xlab=expression(sigma[kappa]))
axis(4,cex=1.2)
mtext(expression(Delta[pi]),4,line=2.5)
qs <- apply((testa3[,,1]),1,function(x){quantile(x,c(0.05,0.95),na.rm=T)})
error.bar(seq(0.05,0.45,0.1)-0.015,rowMeans(testa3[,,1]),qs[2,],qs[1,],lwd=1.2)

points(seq(0.05,0.45,0.1)-0.005,rowMeans(testa3[,,2],na.rm=T),axes=NULL,pch=1,cex=1.2)
qs <- apply((testa3[,,2]),1,function(x){quantile(x,c(0.05,0.95),na.rm=T)})
error.bar(seq(0.05,0.45,0.1)-0.005,rowMeans(testa3[,,2]),qs[2,],qs[1,],lwd=1.2)

points(seq(0.05,0.45,0.1)+0.005,rowMeans(testa3[,,3],na.rm=T),axes=NULL,pch=17,cex=1.2)
qs <- apply((testa3[,,3]),1,function(x){quantile(x,c(0.05,0.95),na.rm=T)})
error.bar(seq(0.05,0.45,0.1)+0.005,rowMeans(testa3[,,3],na.rm=T),qs[2,],qs[1,],lwd=1.2)

points(seq(0.05,0.45,0.1)+0.015,rowMeans(testa3[,,3],na.rm=T),axes=NULL,pch=2,cex=1.2)
qs <- apply((testa3[,,4]),1,function(x){quantile(x,c(0.05,0.95),na.rm=T)})
error.bar(seq(0.05,0.45,0.1)+0.015,rowMeans(testa3[,,4],na.rm=T),qs[2,],qs[1,],lwd=1.2)

abline(h=0,lty=3)

# simulations for number FA dependence --------------

n.preys=3
nsamples=rep(50,n.preys)

# sep sets (roughly,probabilistically) the separation of prey items in multivariate (compositional) space
sep=3
n.preds =1

# set CC var t0 max in the last sim
k=0.45

# set up test array
testa4 <- array(,c(5,50,2))

for(r in 1:50){
  ks=0;
  for (n.fats in seq(5,25,5)){
    ks=ks+1
    As = matrix(rlnorm(n.preys*n.fats,rnorm(n.preys*n.fats,sep)),n.preys,n.fats,byrow=T)
    
    preys <- array(,c(nsamples[1],n.fats,n.preys))
    for (i in 1:n.preys){
      preys[,,i] <- rdirichlet(nsamples[i],As[i,])
    }
    
    preys.ix={};preya={}
    for (i in 1:n.preys){
      preys.ix <- c(preys.ix,rep(i,nsamples[i]))
      preya<- rbind(preya,preys[,,i])
    }
    
    #preyRDA <- rda(alr(preya) ~ as.factor(preys.ix))
    
    #randomly draw fat content
    fc_mean <- (sample(40,n.preys))/10+1
    fc_mean
    # set variance of fat content to 0.4 - this is loosely based on literature
    fc_tau <- rep(1/0.4,n.preys)
    
    # draw mean conversion coeffs - cvar is the variance of mean conversion coeffs around 1; 0.2 seems reasonable  (c.f., Rosen & Tollit histogram)
    cvar = k
    mean_css=matrix(rlnorm(n.fats*n.preys,0,cvar),n.preys,n.fats) # rtnorm is a truncated normal, with lower truncation set to 0
    mean_c=mean_css # just so I can modify mean_cs later without loosing the actual mean
    var_coeffs = matrix((0.03)^2,n.preys,n.fats)  # set the variance to some constant - pretty close, at least for FA whose CC isn't too far from 1 (Iverson 2004 graphs)
    tau_coeffs = 1/var_coeffs
    
    # draw individual coeffs for each predator (e.g., may depend on individual diet status - fasting vs low fat diet versus hight fat diet) from a distribution about mean coeffs
    #cs <- t(matrix(rlnorm(n.preys*n.fats,log(mean_c),sqrt(var_coeffs)),n.fats,n.preys))
    #hist(cs)
    # these are the proportion that each predator eats of each prey
    # population 'mean', Dirichlet alpha notation, hard to interpret but computationally efficient.
    # v(1,2) determines the variance of the alpha coeffs - the higher this variance, the lower the among predator variation in diet proportions.
    v1 = 0.5; v2 = 0.7
    Q = rlnorm(n.preys,rlnorm(n.preys,0,v1),v2)
    # individual proportions
    props=rdirichlet(n.preds,Q)
    # show simulated proportions
    props
    
    # the mixture in matrix algebra (just a matrix product)
    # prey mean is the geometric mean of each column
    mprey <-  clo(t(apply(preys,3,function(x){apply(x,2,function(y){exp(mean(log(y)))})})))
    preda <-  clo(props%*%(fc_mean*mprey*mean_c))
    
    # normalize and transform
    preds <- unclass(data.matrix(alr(preda)))
    
    # transform prey stats
    preym <-  alr(mprey)[,]
    
    m.fats = (n.fats-1)
    # get sums of squares for each prey matrix
    R <- array(,c(m.fats,m.fats,n.preys))
    for (i in 1:n.preys){
      R[,,i]=(cov(alr(preys[,,i]))*(nsamples[i]-1))
      #R[,,i] <- s$u %*% solve(diag(s$d)) %*% t(s$v)
    }
    
    # # set SS matrix for wishart prior on proportions
    # p_tau=0.2
    # 2*(1-pnorm(log(95)/2,0,sqrt(1/p_tau))) # chance that p1=95*p2
    # S = diag(p_tau,n.preys)
    # SS = diag(1,n.preys)
    # set uninformative prior SS matrix for wishart prior alr transformed predator data
    Rnot =diag(1,m.fats)
    
    #zeros=rep(0,n.preys)
    
    
    initials=list(list(
      fc = fc_mean,
      fracs=matrix(mean_c,n.preys,n.fats),
      ps = rep(1/n.preys,n.preys),
      #pmean=rep(0,n.preys),
      #A = Q,#rep(1/n.preys,n.preys),#
      #pprec=diag(1,n.preys), # precision on proportion covariance
      prey.means=preym,
      predprec = diag(1,m.fats),
      prey.precs = array(0.1,c(n.preys,m.fats,m.fats))
      
    ))
    
    # use real mean_c first, then do again with CC set to one
    for(tries in 1:2){
      # look at how much the assumption of a constant c biases the results
      if (tries == 2){
        mean_c=matrix(1,n.preys,n.fats) # set c to one for all CC
        var_coeffs = matrix(k^2,n.preys,n.fats)  # set the variance to some constant - pretty close, at least for FA whose CC isn't too far from 1 (Iverson 2004 graphs)
        tau_coeffs = 1/var_coeffs}
      
      ############# don't try this unless you chose to simualte relatively few FA (else it will take literally forever)
      
      #bugsInits(initials,1,'FattyInits.R')
      
      #bugsData(c('R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','nsamples','preds','preym'),'datas.R')
      datas=list('R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','nsamples','preds','preym')
      
      vars = c('prop')
      
      # compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...
      
      Fattysz<- BRugsFit('Fatty_dir_Single_prop.R', datas, inits=initials, numChains = 1, vars,
                         nBurnin = 1000, nIter = 10000, nThin = 10, coda = T,
                         DIC = F, working.directory = getwd(), digits = 4, seed=NULL,
                         BRugsVerbose = getOption("BRugsVerbose"))
      
      
      Fatty= as.data.frame((Fattysz)[[1]])
      
      prop.ix <- grep('prop',names(Fatty))
      #hist(Fattyz[,prop.ix[1]],20)
      #plot(Fatty[,prop.ix[1]])
      # 
      # pred.ix <- grep('prey',names(Fatty))
      # plot(Fatty[,pred.ix[1]],t='l')
      # postpreys = matrix(colMeans(Fatty[,pred.ix]),3,n.fats,byrow=T)
      # 
      # g=1
      # pps <- ggs_density(cbind(Fatty[,prop.ix[g]],Fattzy[,prop.ix[g]]),c(mean(Fatty[,prop.ix[g]]),mean(Fattzy[,prop.ix[g]]),props[g]))
      # pps
      
      postprops = matrix(colMeans(Fatty[,prop.ix]),n.preds,n.preys,byrow=T)
      #postprops2 = matrix(apply(Fatty[,prop.ix],2,function(x){density(x,adjust=2)$x[density(x,adjust=2)$y==max(density(x,adjust=2)$y)]}),n.preds,n.preys,byrow=T)
      
      #post.props - props
      
      testa4[ks,r,tries] <-dist(alr(rbind(postprops,props)))
    }}}#this one for tries

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

### first plot densities

#### plot for increasing number of vars

plot(seq(5,25,5)-0.5,rowMeans(testa4[,,1],na.rm=T),axes=F,xlab='',ylab='',xlim=c(4.4,25.6),ylim=c(0,4),pch=16,cex=1.2)
axis(1,at=seq(5,25,5),cex=1.2)
title(xlab=expression('Fatty Acids'))
axis(2,cex=1.2)
mtext(expression(Delta[pi]),2,line=2.5)
qs <- apply((testa4[,,1]),1,function(x){quantile(x,c(0.05,0.95),na.rm=T)})
error.bar(seq(5,25,5)-0.5,rowMeans(testa4[,,1]),qs[2,],qs[1,],lwd=1.2)

points(seq(5,25,5)+0.5,rowMeans(testa4[,,2],na.rm=T),axes=NULL,pch=1,cex=1.2)
qs <- apply((testa4[,,2]),1,function(x){quantile(x,c(0.05,0.95),na.rm=T)})
error.bar(seq(5,25,5)+0.5,rowMeans(testa4[,,2]),qs[2,],qs[1,],lwd=1.2)


mp=melt(testa4[,,1])
boxplot(mp$value~as.factor(mp$Var1-0.5),at=seq(1,13,3),xlim=c(0,15),axes=F)

mp=melt(testa4[,,2])
boxplot(mp$value~as.factor(mp$Var1+0.5),at=seq(2,14,3),add=T,axes=F,lty=2)
axis(1,at=seq(1.5,13.5,3),lab=seq(5,25,5))
axis(2)
title(xlab='number of fatty acids',ylab=expression(Delta[pi]))

# simulations for entropy dependence --------------

n.preys=3
nsamples=rep(50,n.preys)

# sep sets (roughly,probabilistically) the separation of prey items in multivariate (compositional) space
sep=3
n.preds =1
n.fats=12

# set CC var t0 max in the last sim
k=0.45

# maximum entropy to standardize
props=rep(1/n.preys,n.preys)
maxent <- -sum(props*log(props))

# set up test array
testa5 <- array(,c(5,20,2))
this.ent <- matrix(,5,20)

for(r in 1:20){
  ks=0;
  for (p in seq(0.5,2,0.5)){
    ks=ks+1
    As = matrix(rlnorm(n.preys*n.fats,rnorm(n.preys*n.fats,sep)),n.preys,n.fats,byrow=T)
    
    preys <- array(,c(nsamples[1],n.fats,n.preys))
    for (i in 1:n.preys){
      preys[,,i] <- rdirichlet(nsamples[i],As[i,])
    }
    
    preys.ix={};preya={}
    for (i in 1:n.preys){
      preys.ix <- c(preys.ix,rep(i,nsamples[i]))
      preya<- rbind(preya,preys[,,i])
    }
    
    #preyRDA <- rda(alr(preya) ~ as.factor(preys.ix))
    
    #randomly draw fat content
    fc_mean <- (sample(40,n.preys))/10+1
    fc_mean
    # set variance of fat content to 0.4 - this is loosely based on literature
    fc_tau <- rep(1/0.4,n.preys)
    
    # draw mean conversion coeffs - cvar is the variance of mean conversion coeffs around 1; 0.2 seems reasonable  (c.f., Rosen & Tollit histogram)
    cvar = k
    mean_css=matrix(rlnorm(n.fats*n.preys,0,cvar),n.preys,n.fats) # rtnorm is a truncated normal, with lower truncation set to 0
    mean_c=mean_css # just so I can modify mean_cs later without loosing the actual mean
    var_coeffs = matrix((0.03)^2,n.preys,n.fats)  # set the variance to some constant - pretty close, at least for FA whose CC isn't too far from 1 (Iverson 2004 graphs)
    tau_coeffs = 1/var_coeffs
    
    # draw individual coeffs for each predator (e.g., may depend on individual diet status - fasting vs low fat diet versus hight fat diet) from a distribution about mean coeffs
    #cs <- t(matrix(rlnorm(n.preys*n.fats,log(mean_c),sqrt(var_coeffs)),n.fats,n.preys))
    #hist(cs)
    # these are the proportion that each predator eats of each prey
    # population 'mean', Dirichlet alpha notation, hard to interpret but computationally efficient.
    # v(1,2) determines the variance of the alpha coeffs - the higher this variance, the lower the among predator variation in diet proportions.
    v1 = 0.5; v2 = 0.7
    Q = rlnorm(n.preys,c(rep(0,n.preys-1),p),v2)
    rlnorm(n.preys,c(rep(0,n.preys-1),p),v2)# individual proportions
    props=rdirichlet(n.preds,Q)
    # show simulated proportions
    #props
    
    this.ent[ks,r] = -sum(props*log(props))
    
    # the mixture in matrix algebra (just a matrix product)
    # prey mean is the geometric mean of each column
    mprey <-  clo(t(apply(preys,3,function(x){apply(x,2,function(y){exp(mean(log(y)))})})))
    preda <-  clo(props%*%(fc_mean*mprey*mean_c))
    
    # normalize and transform
    preds <- unclass(data.matrix(alr(preda)))
    
    # transform prey stats
    preym <-  alr(mprey)[,]
    
    m.fats = (n.fats-1)
    # get sums of squares for each prey matrix
    R <- array(,c(m.fats,m.fats,n.preys))
    for (i in 1:n.preys){
      R[,,i]=(cov(alr(preys[,,i]))*(nsamples[i]-1))
      #R[,,i] <- s$u %*% solve(diag(s$d)) %*% t(s$v)
    }
    
    # # set SS matrix for wishart prior on proportions
    # p_tau=0.2
    # 2*(1-pnorm(log(95)/2,0,sqrt(1/p_tau))) # chance that p1=95*p2
    # S = diag(p_tau,n.preys)
    # SS = diag(1,n.preys)
    # set uninformative prior SS matrix for wishart prior alr transformed predator data
    Rnot =diag(1,m.fats)
    
    #zeros=rep(0,n.preys)
    
    
    initials=list(list(
      fc = fc_mean,
      fracs=matrix(mean_c,n.preys,n.fats),
      ps = rep(1/n.preys,n.preys),
      #pmean=rep(0,n.preys),
      #A = Q,#rep(1/n.preys,n.preys),#
      #pprec=diag(1,n.preys), # precision on proportion covariance
      prey.means=preym,
      predprec = diag(1,m.fats),
      prey.precs = array(0.1,c(n.preys,m.fats,m.fats))
      
    ))
    
    # use real mean_c first, then do again with CC set to one
    for(tries in 1:2){
      # look at how much the assumption of a constant c biases the results
      if (tries == 2){
        mean_c=matrix(1,n.preys,n.fats) # set c to one for all CC
        var_coeffs = matrix(k^2,n.preys,n.fats)  # set the variance to some constant - pretty close, at least for FA whose CC isn't too far from 1 (Iverson 2004 graphs)
        tau_coeffs = 1/var_coeffs}
      
      ############# don't try this unless you chose to simualte relatively few FA (else it will take literally forever)
      
      #bugsInits(initials,1,'FattyInits.R')
      
      #bugsData(c('R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','nsamples','preds','preym'),'datas.R')
      datas=list('R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','nsamples','preds','preym')
      
      vars = c('prop')
      
      # compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...
      
      Fattysz<- BRugsFit('Fatty_dir_Single_prop.R', datas, inits=initials, numChains = 1, vars,
                         nBurnin = 1000, nIter = 10000, nThin = 10, coda = T,
                         DIC = F, working.directory = getwd(), digits = 4, seed=NULL,
                         BRugsVerbose = getOption("BRugsVerbose"))
      
      
      Fatty= as.data.frame((Fattysz)[[1]])
      
      prop.ix <- grep('prop',names(Fatty))
      #hist(Fattyz[,prop.ix[1]],20)
      #plot(Fatty[,prop.ix[1]])
      # 
      # pred.ix <- grep('prey',names(Fatty))
      # plot(Fatty[,pred.ix[1]],t='l')
      # postpreys = matrix(colMeans(Fatty[,pred.ix]),3,n.fats,byrow=T)
      # 
      # g=1
      # pps <- ggs_density(cbind(Fatty[,prop.ix[g]],Fattzy[,prop.ix[g]]),c(mean(Fatty[,prop.ix[g]]),mean(Fattzy[,prop.ix[g]]),props[g]))
      # pps
      
      postprops = matrix(colMeans(Fatty[,prop.ix]),n.preds,n.preys,byrow=T)
      #postprops2 = matrix(apply(Fatty[,prop.ix],2,function(x){density(x,adjust=2)$x[density(x,adjust=2)$y==max(density(x,adjust=2)$y)]}),n.preds,n.preys,byrow=T)
      
      #post.props - props
      
      testa5[ks,r,tries] <-dist(alr(rbind(postprops,props)))
    }}}#this one for tries


## convert for plots
y<-melt(log(testa5))
x<-as.vector(this.ent/maxent)
ent.test <- as.data.frame(cbind(y,x))
iCC.lm <- lm(value~x+as.factor(Var3) ,data=ent.test )

c.lim<-as.data.frame(predict(iCC.lm, ent.test , level=0.95, interval="prediction"))
ent.test <- as.data.frame(cbind(ent.test,c.lim))

# PLOT WITH REGRESSION LINE, CONFIDENCE INTERVAL AND PREDICTION INTERVAL
p0 <- ggplot(ent.test , aes(x=x , y=value,shape=factor(Var3) )) + 
  labs(x='E',y=expression(log(Delta[pi])))+
  geom_point(cex=3) + facet_grid(. ~ Var3)+
  geom_smooth(method = 'lm',fullrange=F, alpha = 0.2) +
  geom_ribbon(aes(y = fit, ymin = lwr, ymax = upr),fill=1,alpha = 0.1,linetype=0) +
  scale_fill_manual('Interval', values = c( 'black','black')) +
  scale_linetype_manual(values = 1:2)+
  scale_colour_manual(values = c( 'black','black'))+
theme(panel.grid.major=element_line(colour = 1),
      panel.grid.minor=element_line(colour = NA),
      panel.background=element_rect(fill=NA),
      legend.position =  c(0.85, 0.85),
      axis.ticks = element_line(colour=1),
      axis.text = element_text(colour=1,size = rel(1.5)),
      axis.title = element_text(size = rel(1.5)))


# simulations for posterior p values dependence --------------

n.preys=3
nsamples=rep(50,n.preys)

# sep sets (roughly,probabilistically) the separation of prey items in multivariate (compositional) space
sep=3
n.preds =1
n.fats=12

# set CC var t0 max in the last sim
k=0.45

# maximum entropy to standardize
props=rep(1/n.preys,n.preys)
maxent <- -sum(props*log(props))

# set up test array
testa_p <- array(,c(100,n.preys,2))
this.ent_p <- matrix(,100,n.preys)
ks=0;

for(r in 1:20){

  for (p in seq(0.5,2,0.5)){
    ks=ks+1
    As = matrix(rlnorm(n.preys*n.fats,rnorm(n.preys*n.fats,sep)),n.preys,n.fats,byrow=T)
    
    preys <- array(,c(nsamples[1],n.fats,n.preys))
    for (i in 1:n.preys){
      preys[,,i] <- rdirichlet(nsamples[i],As[i,])
    }
    
    preys.ix={};preya={}
    for (i in 1:n.preys){
      preys.ix <- c(preys.ix,rep(i,nsamples[i]))
      preya<- rbind(preya,preys[,,i])
    }
    
    #preyRDA <- rda(alr(preya) ~ as.factor(preys.ix))
    
    #randomly draw fat content
    fc_mean <- (sample(40,n.preys))/10+1
    fc_mean
    # set variance of fat content to 0.4 - this is loosely based on literature
    fc_tau <- rep(1/0.4,n.preys)
    
    # draw mean conversion coeffs - cvar is the variance of mean conversion coeffs around 1; 0.2 seems reasonable  (c.f., Rosen & Tollit histogram)
    cvar = k
    mean_css=matrix(rlnorm(n.fats*n.preys,0,cvar),n.preys,n.fats) # rtnorm is a truncated normal, with lower truncation set to 0
    mean_c=mean_css # just so I can modify mean_cs later without loosing the actual mean
    var_coeffs = matrix((0.03)^2,n.preys,n.fats)  # set the variance to some constant - pretty close, at least for FA whose CC isn't too far from 1 (Iverson 2004 graphs)
    tau_coeffs = 1/var_coeffs
    
    # draw individual coeffs for each predator (e.g., may depend on individual diet status - fasting vs low fat diet versus hight fat diet) from a distribution about mean coeffs
    #cs <- t(matrix(rlnorm(n.preys*n.fats,log(mean_c),sqrt(var_coeffs)),n.fats,n.preys))
    #hist(cs)
    # these are the proportion that each predator eats of each prey
    # population 'mean', Dirichlet alpha notation, hard to interpret but computationally efficient.
    # v(1,2) determines the variance of the alpha coeffs - the higher this variance, the lower the among predator variation in diet proportions.
    v1 = 0.5; v2 = 0.7
    Q = rlnorm(n.preys,c(rep(0,n.preys-1),p),v2)
    rlnorm(n.preys,c(rep(0,n.preys-1),p),v2)# individual proportions
    props=rdirichlet(n.preds,Q)
    # show simulated proportions
    #props
    
    this.ent_p[ks,] = abs(props-0.5)
    
    # the mixture in matrix algebra (just a matrix product)
    # prey mean is the geometric mean of each column
    mprey <-  clo(t(apply(preys,3,function(x){apply(x,2,function(y){exp(mean(log(y)))})})))
    preda <-  clo(props%*%(fc_mean*mprey*mean_c))
    
    # normalize and transform
    preds <- unclass(data.matrix(alr(preda)))
    
    # transform prey stats
    preym <-  alr(mprey)[,]
    
    m.fats = (n.fats-1)
    #for simplicity
    ni <- nsamples
    # get sums of squares for each prey matrix
    R <- array(,c(m.fats,m.fats,n.preys))
    for (i in 1:n.preys){
      R[,,i]=(cov(alr(preys[,,i]))*(nsamples[i]-1))
      #R[,,i] <- s$u %*% solve(diag(s$d)) %*% t(s$v)
    }
    
    # # set SS matrix for wishart prior on proportions
    # p_tau=0.2
    # 2*(1-pnorm(log(95)/2,0,sqrt(1/p_tau))) # chance that p1=95*p2
    # S = diag(p_tau,n.preys)
    # SS = diag(1,n.preys)
    # set uninformative prior SS matrix for wishart prior alr transformed predator data
    Rnot =diag(1,m.fats)
    
    #zeros=rep(0,n.preys)
    
    
    initials=list(list(
      fc = fc_mean,
      fracs=matrix(mean_c,n.preys,n.fats),
      ps = rep(1/n.preys,n.preys),
      #pmean=rep(0,n.preys),
      #A = Q,#rep(1/n.preys,n.preys),#
      #pprec=diag(1,n.preys), # precision on proportion covariance
      prey.means=preym,
      predprec = diag(1,m.fats),
      prey.precs = array(0.1,c(n.preys,m.fats,m.fats))
      
    ))
    
    # use real mean_c first, then do again with CC set to one
    for(tries in 1:2){
      # look at how much the assumption of a constant c biases the results
      if (tries == 2){
        mean_c=matrix(1,n.preys,n.fats) # set c to one for all CC
        var_coeffs = matrix(k^2,n.preys,n.fats)  # set the variance to some constant - pretty close, at least for FA whose CC isn't too far from 1 (Iverson 2004 graphs)
        tau_coeffs = 1/var_coeffs}
      
      ############# don't try this unless you chose to simualte relatively few FA (else it will take literally forever)
      
      #bugsInits(initials,1,'FattyInits.R')
      
      #bugsData(c('R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','nsamples','preds','preym'),'datas.R')
      datas=list('R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','ni','preds','preym')
      
      vars = c('prop')
      
      # compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...
      
      Fattysz<- BRugsFit('Fatty_dir_Single_prop.R', datas, inits=initials, numChains = 1, vars,
                         nBurnin = 1000, nIter = 10000, nThin = 10, coda = T,
                         DIC = F, working.directory = getwd(), digits = 4, seed=NULL,
                         BRugsVerbose = getOption("BRugsVerbose"))
      
      
      Fatty= as.data.frame((Fattysz)[[1]])
      
      prop.ix <- grep('prop',names(Fatty))
      #hist(Fattyz[,prop.ix[1]],20)
      #plot(Fatty[,prop.ix[1]])
      # 
      # pred.ix <- grep('prey',names(Fatty))
      # plot(Fatty[,pred.ix[1]],t='l')
      # postpreys = matrix(colMeans(Fatty[,pred.ix]),3,n.fats,byrow=T)
      # 
      # g=1
      # pps <- ggs_density(cbind(Fatty[,prop.ix[g]],Fattzy[,prop.ix[g]]),c(mean(Fatty[,prop.ix[g]]),mean(Fattzy[,prop.ix[g]]),props[g]))
      # pps
      
      #postprops = abs(rowMeans(apply(Fatty[,prop.ix],1,function(x){x-props})>0)-0.5)
      #postprops2 = matrix(apply(Fatty[,prop.ix],2,function(x){density(x,adjust=2)$x[density(x,adjust=2)$y==max(density(x,adjust=2)$y)]}),n.preds,n.preys,byrow=T)
      
      #post.props - props
      
      testa_p[ks,,tries] <-abs(rowMeans(apply(Fatty[,prop.ix],1,function(x){x-props})>0)-0.5)
    }}}#this one for tries

plot(this.ent_p,testa_p[,,1],xlab=abss)


### try with variable selection and transformation ###########


n.preys=4
nsamples=rep(26,n.preys)

## re-simualte
nsims <- 50
test_VS <- matrix(,6,nsims)
this.nv <- matrix(,2,nsims)
this.dist <- rep(NA,nsims)
this.Ent <- rep(NA,nsims)
this.full.time <- rep(NA,nsims)
this.vs.time <- rep(NA,nsims)

for (k in 38:nsims){

  n.fats.org=n.fats=20
  
# sep sets (roughly,probabilistically) the separation of prey items in multivariate (compositional) space
sep=runif(1,min=0.5,max=4)
As = matrix(rlnorm(n.preys*n.fats,rnorm(n.fats,sep)),n.preys,n.fats,byrow=T)

  #invs <- list(NULL)
  
preys <- array(,c(nsamples[1],n.fats,n.preys))
for (i in 1:n.preys){
  preys[,,i] <- rdirichlet(nsamples[i],As[i,])
  #invs[[i]] <- try(solve(cov(alr(preys[,,i]))),T)
  }


  
preys.ix={};preya={}
for (i in 1:n.preys){
  preys.ix <- c(preys.ix,rep(i,nsamples[i]))
  preya<- rbind(preya,preys[,,i])
}

require('robCompositions')

pr_DB$set_entry(FUN = aDist, names = c("test", "adist"))

dists <- matrix(,nsamples*n.preys,nsamples*n.preys)
for (i in 1:sum(nsamples)){
  for (j in i:sum(nsamples)){
    dists[j,i] <- aDist(preya[i,],preya[j,])
  }}
dista <- as.dist(dists)

# do groupwise dists
gdists <- matrix(,n.preys,n.preys)
for (i in 1:n.preys){
  for (j in i:n.preys){
    gdists[j,i] <- dist(t(apply(preya[preys.ix==i,],2,gmean)),t(apply(preya[preys.ix==j,],2,gmean)),method='adist')
  }}

this.dist[k] <- mean(gdists[gdists>0],na.rm=T)

detach("package:robCompositions", unload=TRUE)

preyRDA <- capscale(dista ~ as.factor(preys.ix),comm=preya)

# note how the proportion of (constrained) RDA axes compared to unconstraned PC axes 
# increases with sep - meaning the among prey variance increases relative to within prey variance.
# summary(preyRDA,display = NULL)
# plot(preyRDA,t='n')
# points(preyRDA,col=preys.ix)
# text(preyRDA,'sp',col=4,pch=3)
# 
# plot(cumsum(sort(clo(rowSums(t(t(cbind(preyRDA$CCA$v)*c(preyRDA$CCA$eig))^2))),decreasing =T)),ylab='cumulative proportion')
# 
# 
# plot(cumsum(sort(clo(rowSums(t(t(cbind(preyRDA$CCA$v,preyRDA$CA$v))*c(preyRDA$CCA$eig,preyRDA$CA$eig))^2)),decreasing =T)),ylab='cumulative proportion')
  
#randomly draw fat content
fc_mean <- (sample(40,n.preys))/10+1
fc_mean
# set variance of fat content to 0.4 - this is loosely based on literature
fc_tau <- rep(1/0.4,n.preys)

# draw mean conversion coeffs - cvar is the variance of mean conversion coeffs around 1; 0.2 seems reasonable  (c.f., Rosen & Tollit histogram)
cvar = 0.45
mean_css=matrix(rlnorm(n.fats*n.preys,0,cvar),n.preys,n.fats.org) # rtnorm is a truncated normal, with lower truncation set to 0
mean_c=mean_css # just so I can modify mean_cs later without loosing the actual mean
var_coeffs = matrix((0.03)^2,n.preys,n.fats.org)  # set the variance to some constant - pretty close, at least for FA whose CC isn't too far from 1 (Iverson 2004 graphs)
tau_coeffs = 1/var_coeffs
#hist(mean_c)
# draw individual coeffs for each predator (e.g., may depend on individual diet status - fasting vs low fat diet versus hight fat diet) from a distribution about mean coeffs
#cs <- t(matrix(rlnorm(n.preys*n.fats,log(mean_c),sqrt(var_coeffs)),n.fats,n.preys))
#hist(cs)
# these are the proportion that each predator eats of each prey
# population 'mean', Dirichlet alpha notation, hard to interpret but computationally efficient.
# v(1,2) determines the variance of the alpha coeffs - the higher this variance, the lower the among predator variation in diet proportions.
v1 = 1; v2 = 0.7
Q = rlnorm(n.preys,rlnorm(n.preys,0,v1),v2)
# individual proportions
n.preds =1    
props=rdirichlet(n.preds,Q)
    # show simulated proportions
    #props
this.Ent[k] <- entropy(props)    

# the mixture in matrix algebra (just a matrix product)
# prey mean is the geometric mean of each column

  mprey <-  clo(t(apply(preys,3,function(x){apply(x,2,function(y){exp(mean(log(y)))})})))
  ##gmprey <-  exp(rowMeans(colMeans(log(preys[,six,]))))
  
  preda = clo(props%*%(fc_mean*mprey*mean_css))
  
  # normalize and transform
  preds <-unclass(alr(preda))
  # transform prey stats
  preym <-  (alr(mprey))[,]
 
  m.fats = n.fats-1
  
  # get sums of squares for each prey matrix
  R <- array(,c(m.fats,m.fats,n.preys))
  ni<-rep(NA,n.preys)
  for (i in 1:n.preys){
    ni[i] <- max(n.fats+1,(nsamples[i]-1))
    R[,,i]=cov(alr(preys[,,i]))*ni[i]
    #R[,,i] <- s$u %*% solve(diag(s$d)) %*% t(s$v)
  }
  
  # set uninformative prior SS matrix for wishart prior alr transformed predator data
  Rnot =diag(1,m.fats)
  
  # # set SS matrix for wishart prior on proportions
  p_tau=0.01
  2*(1-pnorm(log(95)/2,0,sqrt(1/p_tau))) # chance that p1=95*p2
  S = diag(p_tau,n.preys)
  SS = diag(1,n.preys)
  
  zeros=rep(0,n.preys)
  
  initials=list(list(
    fc = fc_mean,
    fracs=matrix(mean_c,n.preys,n.fats),
    #pnorm = matrix(1/n.preys,n.preds,n.preys),
    pmean=rep(1,n.preys),
    #A = Q,#rep(1/n.preys,n.preys),#
    #pprec=diag(1,n.preys), # precision on proportion covariance
    prey.means=preym,
    predprec = diag(1,m.fats),
    prey.precs = array(1,c(n.preys,m.fats,m.fats) )
    
  ))
  
  ############# don't try this unless you chose to simualte relatively few FA (else it will take literally forever)
  
  bugsInits(initials,1,'FattyInits.R')
  
  bugsData(c('ni','SS','S','zeros','R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','preds','preym'),'datas.R')
  datas=list('ni','SS','S','zeros','R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','preds','preym')
  
  vars = c('prop')
  
  # compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...
  #Fatty_dir_Single_prop
  ptm <- proc.time()
  Fattysz<- BRugsFit('Fatty_dir_Single_prop.R', data='datas.R', inits=initials, numChains = 1, vars,
                     nBurnin = 1000, nIter = 10000, nThin = 10, coda = T,
                     DIC = F, working.directory = getwd(), digits = 4, seed=NULL,
                     BRugsVerbose = getOption("BRugsVerbose"))
  this.full.time[k] = proc.time() - ptm
  
  Fatty= as.data.frame((Fattysz)[[1]])
  
  prop.ix <- grep('prop',names(Fatty))
  #hist(Fattyz[,prop.ix[1]],20)
  #plot(Fatty[,prop.ix[1]],t='l')
  # 
    g=1
   pps <- ggs_density(Fatty[,prop.ix],c(colMeans(props)))
  # pps <- ggs_density(cbind(Fatty1[,prop.ix[g]],Fatty2[,prop.ix[g]],Fatty3[,prop.ix[g]],Fatty4[,prop.ix[g]],Fatty5[,prop.ix[g]]),c(mean(Fatty1[,prop.ix[g]]),mean(Fatty2[,prop.ix[g]]),mean(Fatty3[,prop.ix[g]]),mean(Fatty4[,prop.ix[g]]),mean(Fatty5[,prop.ix[g]]),props[g]))
  #   pps
  
  postprops = matrix(colMeans(Fatty[,prop.ix]),n.preds,n.preys,byrow=T)
  #postprops2 = matrix(apply(Fatty[,prop.ix],2,function(x){density(x,adjust=2)$x[density(x,adjust=2)$y==max(density(x,adjust=2)$y)]}),n.preds,n.preys,byrow=T)
  
  test_VS[1,k] <- dist(alr(rbind(postprops,props)))
  
  # or ignore them
  mean_c = mean_css*0+1
  tau_coeffs = mean_c 
  
  bugsData(c('R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','ni','preds','preym'),'datas.R')
  
  Fattysz<- BRugsFit('Fatty_dir_Single_prop.R', data='datas.R', inits=initials, numChains = 1, vars,
                     nBurnin = 1000, nIter = 10000, nThin = 10, coda = T,
                     DIC = F, working.directory = getwd(), digits = 4, seed=NULL,
                     BRugsVerbose = getOption("BRugsVerbose"))
  
  
  Fatty= as.data.frame((Fattysz)[[1]])
  
  prop.ix <- grep('prop',names(Fatty))
  
  postprops = matrix(colMeans(Fatty[,prop.ix]),n.preds,n.preys,byrow=T)
  
  test_VS[2,k] <- dist(alr(rbind(postprops,props)))
  
  
  
## variable transformation  -----
 
  
  nv <- max(which(cumsum(c(preyRDA$CCA$eig,preyRDA$CA$eig))/(sum(preyRDA$CCA$eig)+sum(preyRDA$CA$eig))>0.99)[1],3)
  this.nv[1,k]<- nv
  
  
preym <-  matrix(,n.preys,nv)
for (p in 1:n.preys){
  
  preym[p,] <- clo(mprey[p,]*mean_css[p,])%*%cbind(preyRDA$CCA$v,preyRDA$CA$v[,min(1,(nv-(n.preys-1))):(nv-(n.preys-1))])
}

preds = clo((fc_mean*props) %*%(mprey*mean_css))
preds = t(apply((preds),1,function(x){t(x)%*%cbind(preyRDA$CCA$v,preyRDA$CA$v[,min(1,(nv-(n.preys-1))):(nv-(n.preys-1))])}))

#preds/(props%*%preym)
#mean(preds/(clo(fc_mean*props)%*%(preym)))

# plot(preym[,1:2])
# points(t(preds[1:2]),col=4,pch=16)
# points(t((clo(props)%*%(preym))[1:2]),col=5)

n.fats=m.fats=nv
#preyst <- preys[,1:2,]
# get sums of squares for each prey matrix
R <- array(,c(m.fats,m.fats,n.preys))
for (i in 1:n.preys){
  R[,,i]=(cov(t(apply(clo(preys[,,i]),1,function(x){clo(t(x)*mean_c[i,])%*%cbind(preyRDA$CCA$v,preyRDA$CA$v[,min(1,(nv-(n.preys-1))):(nv-(n.preys-1))])})))*(nsamples[i]-1))
  #R[,,i] <- s$u %*% solve(diag(s$d)) %*% t(s$v)
  #preyst[,,i] <- t(apply(preys[,,i],1,function(x){clo(t(x))%*%preyRDA$CCA$v}))
}

# set uninformative prior SS matrix for wishart prior alr transformed predator data
Rnot =diag(1e-6,m.fats)
# # set SS matrix for wishart prior on proportions
p_tau=0.05
2*(1-pnorm(log(95)/2,0,sqrt(1/p_tau))) # chance that p1=95*p2
S = diag(p_tau,n.preys)
SS = diag(1,n.preys)

zeros=rep(0,n.preys)

#zeros=rep(0,n.preys)

iinitials=list(list(
  #fc = fc_mean,
  #fracs=matrix(mean_c,n.preys,n.fats),
  #pnorm = matrix(1/n.preys,n.preds,n.preys),
  pmean=rep(0,n.preys),
  #A = Q,#rep(1/n.preys,n.preys),#
  pprec=diag(1,n.preys), # precision on proportion covariance
  prey.means=preym,
  predprec = diag(1,m.fats),
  prey.precs = array(1,c(n.preys,m.fats,m.fats))
  
))


# bugsInits(iinitials,1,'FattsInits.R')
# 
# bugsData(c('SS','S','zeros','R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','nsamples','preds','preym'),'datass.R')

datass=list('SS','S','zeros','fc_mean','fc_tau','R','Rnot','n.preys','n.preds','n.fats','m.fats','nsamples','preds','preym')

vars = c('prop')


# compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...

Fattysz<- BRugsFit('Fatty_VS.R', datass, inits=iinitials, numChains = 1, vars,
                   nBurnin = 10000, nIter = 10000, nThin = 10, coda = T,
                   DIC = F, working.directory = getwd(), digits = 4, seed=NULL,
                   BRugsVerbose = getOption("BRugsVerbose"))

Fatty= as.data.frame((Fattysz)[[1]])

prop.ix <- grep('prop',names(Fatty))
#hist(Fattyz[,prop.ix[1]],20)
#plot(Fatty[,prop.ix[2]],t='l')
# 
# g=1
# pps <- ggs_density(Fatty[,prop.ix],c(colMeans(Fatty[,prop.ix]),(props)))
# pps

postprops = matrix(colMeans(Fatty[,prop.ix]),n.preds,n.preys,byrow=T)

# pred.ix <- grep('pred',names(Fatty))
# plot(data.matrix(Fatty[,pred.ix[1:2]]))
# points(t(preds[,1:2]),col=2,pch=18,cex=4)
# points(t((clo(fc_mean*props)%*%(preym))[,1:2]),col=2,pch=19,cex=3)
# 
# preys.ix <- grep('prey',names(Fatty))
# points(col=3,rbind(data.matrix(Fatty[,preys.ix[1:2]]),data.matrix(Fatty[,preys.ix[6:7]]),data.matrix(Fatty[,preys.ix[11:12]])))
# points(rbind(preyst[,,1],preyst[,,2],preyst[,,3]),col=4)
# points(preym[,1:2],col=5,pch=17)

test_VS[3,k] <- dist(alr(rbind(postprops,props)))

## now do it again but without mean_c

Rnot =diag(0.00001,m.fats)

R <- array(,c(m.fats,m.fats,n.preys))
for (i in 1:n.preys){
  R[,,i]=(cov(t(apply(clo(preys[,,i]),1,function(x){clo(t(x))%*%cbind(preyRDA$CCA$v,preyRDA$CA$v[,min(1,(nv-(n.preys-1))):(nv-(n.preys-1))])})))*(nsamples[i]-1))
  #R[,,i] <- s$u %*% solve(diag(s$d)) %*% t(s$v)
  #preyst[,,i] <- t(apply(preys[,,i],1,function(x){clo(t(x))%*%preyRDA$CCA$v}))
}
R=R

preym <-  matrix(,n.preys,nv)
for (p in 1:n.preys){
  
  preym[p,] <- clo(fc_mean[p]*mprey[p,])%*%cbind(preyRDA$CCA$v,preyRDA$CA$v[,min(1,(nv-(n.preys-1))):(nv-(n.preys-1))])
}

Fattysz<- BRugsFit('Fatty_VS.R', datass, inits=iinitials, numChains = 1, vars,
                   nBurnin = 10000, nIter = 10000, nThin = 10, coda = T,
                   DIC = F, working.directory = getwd(), digits = 4, seed=NULL,
                   BRugsVerbose = getOption("BRugsVerbose"))

Fatty= as.data.frame((Fattysz)[[1]])

prop.ix <- grep('prop',names(Fatty))

postprops = matrix(colMeans(Fatty[,prop.ix]),n.preds,n.preys,byrow=T)

test_VS[4,k] <- dist(alr(rbind(postprops,props)))

### try with variable selection ###########

#plot(cumsum(sort(clo(rowSums(t(t(cbind(preyRDA$CCA$v)*c(preyRDA$CCA$eig))^2))),decreasing =T)),ylab='cumulative proportion')
#nv <- max(which(cumsum(sort(clo(rowSums(t(t(cbind(preyRDA$CCA$v)*c(preyRDA$CCA$eig))^2))),decreasing =T))>0.99)[1],3)
  nv <- max(which(cumsum(sort(clo(rowSums(t(t(cbind(preyRDA$CCA$v,preyRDA$CA$v))*c(preyRDA$CCA$eig,preyRDA$CA$eig))^2)),decreasing =T))>0.99)[1],3)
this.nv[2,k] <-nv

#plot(cumsum(sort(clo(rowSums(t(t(cbind(preyRDA$CCA$v,preyRDA$CA$v))*c(preyRDA$CCA$eig,preyRDA$CA$eig))^2)),decreasing =T)),ylab='cumulative proportion')

# plot contibution to axes of each variable
#sv = sort(clo(rowSums(t(t(preyRDA$CCA$v)*preyRDA$CCA$eig)^2)),decreasing =T,index.return=T)
sv = sort(clo(rowSums(t(t(cbind(preyRDA$CCA$v,preyRDA$CA$v))*c(preyRDA$CCA$eig,preyRDA$CA$eig))^2)),decreasing =T,index.return=T)
  # number of variable to keep
six <- sv$ix[1:nv]
# 
# vs <- var_select(preya,preys.ix)  
#  
#   nv <- max(vs$nv,3)  
# if(vs$nv<3)
# {six <- vs$six[1:3] }else{
#   six <- vs$vars}
#  
#   six <-vs$six[1:nv]
# sorted.vars = sort(clo(rowSums(t(t(cbind(preyRDA$CCA$v,preyRDA$CA$v))*c(preyRDA$CCA$eig,preyRDA$CA$eig))^2)),decreasing =T,index.return=T)
# # number of variable to keep
# nv <- 10
# six <- sorted.vars$ix[1:nv]

mprey <-  clo(t(apply(preys,3,function(x){apply(x,2,function(y){exp(mean(log(y)))})})))
##gmprey <-  exp(rowMeans(colMeans(log(preys[,six,]))))

preda = clo(props%*%(fc_mean*mprey*mean_css))

# normalize and transform
preds <-unclass(t(alr(preda[,six])))
# transform prey stats
preym <-  (alr(mprey[,six]))[,]

# plot these 
# PCA
# scores = princomp(alr(preya[,six]))$scores[,1:2]
# plot(scores,col=preys.ix)

mean_c = mean_css[,six]
tau_coeffs = 1/var_coeffs[,six]

n.fats = nv
m.fats = n.fats-1 
# get sums of squares for each prey matrix
R <- array(,c(m.fats,m.fats,n.preys))
for (i in 1:n.preys){
  R[,,i]=cov(alr(preys[,six,i]))*(nsamples[i]-1)
  #R[,,i] <- s$u %*% solve(diag(s$d)) %*% t(s$v)
}

# set uninformative prior SS matrix for wishart prior alr transformed predator data
Rnot =diag(0.1,m.fats)

# # set SS matrix for wishart prior on proportions
p_tau=0.2
2*(1-pnorm(log(95)/2,0,sqrt(1/p_tau))) # chance that p1=95*p2
S = diag(p_tau,n.preys)
SS = diag(1,n.preys)

zeros=rep(0,n.preys)

initials=list(list(
  fc = fc_mean,
  fracs=matrix(mean_c,n.preys,n.fats),
  #pnorm = matrix(1/n.preys,n.preds,n.preys),
  pmean=rep(0,n.preys),
  #A = Q,#rep(1/n.preys,n.preys),#
  #pprec=diag(1,n.preys), # precision on proportion covariance
  prey.means=preym,
  predprec = diag(1,m.fats),
  prey.precs = array(0.1,c(n.preys,m.fats,m.fats) )
))
  
  ############# don't try this unless you chose to simualte relatively few FA (else it will take literally forever)
  
  #bugsInits(initials,1,'FattyInits.R')
  
  #bugsData(c('R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','ni','preds','preym'),'datas.R')
  datas=list('R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','ni','preds','preym')
  
  vars = c('prop')
  
  # compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...
  ptm <- proc.time()                    #Fatty_dir_Single_prop
  Fattysz<- BRugsFit('Fatty_dir_Single_prop.R', data=datas, inits=initials, numChains = 1, vars,
                     nBurnin = 1000, nIter = 10000, nThin = 10, coda = T,
                     DIC = F, working.directory = getwd(), digits = 4, seed=NULL,
                     BRugsVerbose = getOption("BRugsVerbose"))
  this.vs.time[k] = proc.time() - ptm
  
  Fatty= as.data.frame((Fattysz)[[1]])
  
  prop.ix <- grep('prop',names(Fatty))
  #hist(Fattyz[,prop.ix[1]],20)
  #plot(Fatty[,prop.ix[1]],t='l')
  # 
#   g=1
# pps <- ggs_density(Fatty[,prop.ix],c(colMeans(Fatty[,prop.ix]),colMeans(props)))
# pps <- ggs_density(cbind(Fatty1[,prop.ix[g]],Fatty2[,prop.ix[g]],Fatty3[,prop.ix[g]],Fatty4[,prop.ix[g]],Fatty5[,prop.ix[g]]),c(mean(Fatty1[,prop.ix[g]]),mean(Fatty2[,prop.ix[g]]),mean(Fatty3[,prop.ix[g]]),mean(Fatty4[,prop.ix[g]]),mean(Fatty5[,prop.ix[g]]),props[g]))
#   pps
  
  postprops = matrix(colMeans(Fatty[,prop.ix]),n.preds,n.preys,byrow=T)
  #postprops2 = matrix(apply(Fatty[,prop.ix],2,function(x){density(x,adjust=2)$x[density(x,adjust=2)$y==max(density(x,adjust=2)$y)]}),n.preds,n.preys,byrow=T)
  
test_VS[5,k] <- dist(alr(rbind(postprops,props)))
  

# or ignore them
mean_c = mean_css[,six]*0+1
tau_coeffs = mean_c 

bugsData(c('S','zeros','R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','ni','preds','preym'),'datas.R')

Fattysz<- BRugsFit('Fatty_dir_Single_prop.R', data='datas.R', inits=initials, numChains = 1, vars,
                   nBurnin = 1000, nIter = 10000, nThin = 10, coda = T,
                   DIC = F, working.directory = getwd(), digits = 4, seed=NULL,
                   BRugsVerbose = getOption("BRugsVerbose"))


Fatty= as.data.frame((Fattysz)[[1]])

prop.ix <- grep('prop',names(Fatty))

postprops = matrix(colMeans(Fatty[,prop.ix]),n.preds,n.preys,byrow=T)

test_VS[6,k] <- dist(alr(rbind(postprops,props)))

}

par(mar=c(5,5,4,2))
boxplot(t(test_VS[c(5,6),1:40])-t(test_VS[c(1,2),1:40]),axes=F)
axis(1,c(1,2),c(expression(kappa[k]),expression(kappa[u])))
abline(h=0,cex=1.2,lty=3)
ggs_density(t(test_VS[,1:53]))
ggs_density(Fatty[,prop.ix],c(props,postprops[,1:4]))


#### try anova formualtion ------

### first set of predtors
v1 = 0.5; v2 = 0.7
Q = rlnorm(n.preys,rlnorm(n.preys,0,v1),v2)
# individual proportions
props=rdirichlet(n.preds,Q)
# show simulated proportions
props


# second set
v1 = 0.5; v2 = 0.7
Q = rlnorm(n.preys,rlnorm(n.preys,0,v1),v2)
# individual proportions
props2=rdirichlet(n.preds,Q)
# show simulated proportions
props <- rbind(props,props2 )

# define size dummy vector
sizec <- c(rep(0,n.preds/2),rep(1,n.preds/2))

preda = clo(props%*%(fc_mean*mprey*mean_css))

# normalize and transform
preds <-(alr(preda[,six]))[,]
# transform prey stats
preym <-  (alr(mprey[,six]))[,]

# plot these 
# PCA
scores = princomp(alr(preya[,six]))$scores[,1:2]
plot(scores,col=preys.ix)

mean_c = mean_css[,six]
tau_coeffs = 1/var_coeffs[,six]

n.fats = nv
m.fats = n.fats-1 
# get sums of squares for each prey matrix
R <- array(,c(m.fats,m.fats,n.preys))
for (i in 1:n.preys){
  R[,,i]=cov(alr(preys[,six,i]))*(nsamples[i]-1)
  #R[,,i] <- s$u %*% solve(diag(s$d)) %*% t(s$v)
}

# set uninformative prior SS matrix for wishart prior alr transformed predator data
Rnot =diag(1,m.fats)

# # set SS matrix for wishart prior on proportions
p_tau=0.2
2*(1-pnorm(log(95)/2,0,sqrt(1/p_tau))) # chance that p1=95*p2
S = diag(p_tau,n.preys)
SS = diag(1,n.preys)

zeros=rep(0,n.preys)

initials=list(list(
  fc = fc_mean,
  fracs=matrix(mean_c,n.preys,n.fats),
  pnorm = matrix(1/n.preys,n.preds,n.preys),
  pmean=rep(0,n.preys),beta=rep(0,n.preys),
  #A = Q,#rep(1/n.preys,n.preys),#
  pprec=diag(1,n.preys), # precision on proportion covariance
  prey.means=preym,
  predprec = diag(1,m.fats),
  prey.precs = array(0.1,c(n.preys,m.fats,m.fats) )
))

############# don't try this unless you chose to simualte relatively few FA (else it will take literally forever)

bugsInits(initials,1,'FattyInits.R')

bugsData(c('sizec','S','SS','zeros','R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','nsamples','preds','preym'),'datas.R')
datas=list('sizec','S','SS','zeros','R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','nsamples','preds','preym')

vars = c('prop','beta')

# compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...

Fattysz<- BRugsFit('Fatty_dir_pop_and_Ind_anova.R', datas, inits=initials, numChains = 1, vars,
                   nBurnin = 1000, nIter = 10000, nThin = 10, coda = T,
                   DIC = F, working.directory = getwd(), digits = 4, seed=NULL,
                   BRugsVerbose = getOption("BRugsVerbose"))


Fattya= as.data.frame((Fattysz)[[1]])

beta.ix <- grep('beta',names(Fattya))
#hist(Fattyz[,prop.ix[1]],20)
plot(Fattya[,beta.ix[4]])

prop.ix <- grep('prop',names(Fattya))

postprops = matrix(colMeans(Fattya[,prop.ix]),n.preds,n.preys,byrow=T)

## now incorporate SI data ----

require(MASS)

isos=2
n_SI=rep(50,n.preys)
prey.m=matrix(c(-15,-20,-24,7,10,2),2,n.preys,byrow=T)

prey.cov = array(,c(2,2,n.preys))
prey.cov[,,1]=diag(1,2)
prey.cov[,,2]=diag(1.4,2)
prey.cov[,,3]=diag(0.7,2)

preys.SI <- array(NA,c(n_SI[1],isos,n.preys))
for (i in 1:n.preys){
  preys.SI[,,i] <-mvrnorm(n_SI[i],prey.m[,i],prey.cov[,,i])
}

# use publsihed data from MCCutchan et al 2003 for C and N
mean_cs = c(0.4,2.3)
tau_cs =1/c(0.2,0.61)^2 # 1.61

# histograms in MCCutchan et al 2003 look very normal
cs.SI <- t(matrix(rnorm(n.preds*isos,mean_cs,sqrt(1/tau_cs)),2,n.preds))

# the mixture in matrix algebra (just a matrix product)
preym.SI <-  t(apply(preys.SI,3,colMeans))
preds.SI <-  props %*% preym.SI+cs.SI


## first some data and inits

Rnot_SI = diag(1,isos)

# get sums of squares for each prey matrix
R_SI <- array(,c(isos,isos,n.preys))
for (i in 1:n.preys){
  R_SI[,,i]=cov(preys.SI[,,i])*(n_SI[i]-1)  
}

joint_initials=list(list(
 
  pnorm = matrix(0,n.preds,n.preys),
  pmean=rep(0,n.preys),
  pprec=diag(10,n.preys), # precision on proportion covariance
  
  # fatty acid inits
  fc = fc_mean,
  fracs=matrix(mean_c,n.preys,n.fats),
  prey.means=preym,
  predprec = diag(0.01,m.fats),
  prey.precs = array(1,c(n.preys,m.fats,m.fats)),
  
  #isotope inits
  cs=array(mean_cs,c(n.preds,n.preys,isos)),
  prey.means_SI=preym.SI,
  predprec_SI = diag(0.01,isos),
  prey.precs_SI = array(1,c(n.preys,isos,isos))
))

###only to write out data and inits to file

#bugsInits(joint_initials,1,'JointInits.R')
#bugsData(list('SS','R','R_SI','fc_mean','fc_tau','mean_c','tau_coeffs','mean_cs','tau_cs','Rnot','Rnot_SI','S','zeros','n.preys','n.preds','n.fats','m.fats','isos','nsamples','n_SI','preds','preym','preds.SI','preym.SI'),'jointdata.R')

joint_datas=list('SS','R','R_SI','fc_mean','fc_tau','mean_c','tau_coeffs','mean_cs','tau_cs','Rnot','Rnot_SI','S','zeros','n.preys','n.preds','n.fats','m.fats','isos','nsamples','n_SI','preds','preym','preds.SI','preym.SI')

# RUN MODEL

vars = c('prop','pmean')

Joint <- BRugsFit('Fatty_and_SI.R', joint_datas, inits=joint_initials, numChains = 1, vars,
                  nBurnin = 1000, nIter = 10000, nThin = 1, coda = T,
                  DIC = F, working.directory = getwd(), digits = 3)

Joint = as.data.frame((Joint)[[1]])
names(Joint)

prop.ix <- grep('prop',names(Joint))
hist(Joint[,prop.ix[1]],20)
plot(Joint[,prop.ix[3]])
post.props = matrix(colMeans(Joint[,prop.ix]),n.preds,n.preys,byrow=T)
post.props - props

pm.ix <- grep('pm',names(Joint))
hist(Joint[,pm.ix[1]],20)
plot(Joint[,pm.ix[1]])

colMeans(props)
exp(colMeans(Joint[,pm.ix]))/sum(exp(colMeans(Joint[,pm.ix])))

