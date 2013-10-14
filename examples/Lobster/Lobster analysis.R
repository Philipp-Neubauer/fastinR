# simulate fatty acid proportions in 3 prey species ----
require(MCMCpack)
require(compositions)
require(BRugs)
require(vegan)
require(msm)
require(reshape2)
require(ggplot2)

source('~/Work/Dropbox/Projects private files/Fatty Acid model/Multiplot.R')
source('gg_density.R')

gmean <- function(x){exp(mean(log(x)))}

# try with Lobster Reserve data -----------

RLR <- read.csv("/home/philbobsqp/Work/Dropbox/Projects private files/Fatty Acid model/Lobster diet/Reserve_FA.csv",header=F)

# Lobster means and SD
RLR.means <- clo(RLR[2:nrow(RLR),1])
RLR.sd <- RLR[2:nrow(RLR),2]*sqrt(RLR[1,1])/100

# Prey means and SD
RLR.prey.means <- apply(RLR[2:nrow(RLR),seq(3,ncol(RLR),2)],2,clo)
RLR.prey.sd <- t(apply(data.matrix(RLR[2:nrow(RLR),seq(4,ncol(RLR),2)]),1,function(x){x*sqrt(data.matrix(RLR[1,seq(4,ncol(RLR),2)]))}))/100

# pool data 

RLR.means.new <- cbind(RLR.prey.means[,1:2],apply(RLR.prey.means[,3:4],1,gmean),apply(RLR.prey.means[,5:7],1,gmean),apply(RLR.prey.means[,8:9],1,gmean),RLR.prey.means[,10])
# just use mean sd - numbers are comparable and it's just an illustration...
RLR.sd.new <- cbind(RLR.prey.sd[,1:2],rowMeans(RLR.prey.sd[,3:4]),rowMeans(RLR.prey.sd[,5:7]),rowMeans(RLR.prey.sd[,8:9]),RLR.prey.sd[,10])

# replace zeros in sds
for (i in 1:nrow(RLR.sd.new))
  RLR.sd.new[i,RLR.sd.new[i,]==0] <- 1e-8#min(RLR.sd.new[i,RLR.sd.new[i,]!=0])

n.preys=ncol(RLR.means.new)
#get fat contents
fc.mean.new = rep(1,n.preys)# c(1.4,1.3,4.8,1,0.5,0.3)
fc.sd.new = c(fc.mean.new /1)

RLR.means.new[RLR.means.new==0]=0.0001

nsample = as.numeric(RLR[1,seq(3,ncol(RLR-1),2)])
nsamples = c(nsample[1:2],sum(nsample[3:4]),sum(nsample[5:7]),sum(nsample[8:9]),nsample[10])

RLR.sim <- rep(NA,nrow(RLR.means.new)+1)
for (i in 1:ncol(RLR.means.new))
  RLR.sim <- rbind(RLR.sim,cbind(rep(i,50),abs(mvrnorm(50,RLR.means.new[,i],diag(RLR.sd.new[,i])^2))))
RLR.sim <- (RLR.sim[-1,])
RLR.sim[,2:ncol(RLR.sim)] <- clo(RLR.sim[,2:ncol(RLR.sim)])


require('robCompositions')
dists <- matrix(,nrow(RLR.sim),nrow(RLR.sim))
for (i in 1:nrow(RLR.sim)){
  for (j in i:nrow(RLR.sim)){
    dists[j,i] <- aDist(RLR.sim[i,],RLR.sim[j,])
  }}
detach("package:robCompositions", unload=TRUE)
dista <- as.dist(dists)


RLR.RDA <- capscale(dista~as.factor(RLR.sim[,1]),comm=RLR.sim[,2:ncol(RLR.sim)])
#plot(RLR.RDA,t='n',xlim=c(-0.5,0.5),ylim=c(-1,1))
plot(t(apply(t(RLR.means.new)%*%RLR.RDA$CCA$v[,1:2],1,function(x){x*c(50,100)})),pch=1:6,ylim=c(-2,2.5))
points(clo(t(RLR.means))%*%RLR.RDA$CCA$v[,1:2]*c(50,100),pch=16)
legend(0.2,2.5,c('H.rubra','T.undulatus','Sea Urchins','Ascidians','Brown Algae','Red Algae'),pch=c(1:6,16,17))

plot(cumsum(sort(clo(rowSums(t(t(cbind(RLR.RDA$CCA$v)*c(RLR.RDA$CCA$eig))^2))),decreasing =T)),ylab='cumulative proportion',xlab='n')

which(cumsum(sort(clo(rowSums(t(t(cbind(RLR.RDA$CCA$v)*c(RLR.RDA$CCA$eig))^2))),decreasing =T))>0.95)[1]

sv = sort(clo(rowSums(t(t(cbind(RLR.RDA$CCA$v,RLR.RDA$CA$v))*c(RLR.RDA$CCA$eig,RLR.RDA$CA$eig))^2)),decreasing =T,index.return=T)
nv <- which(cumsum(sort(clo(rowSums(t(t(cbind(RLR.RDA$CCA$v)*c(RLR.RDA$CCA$eig))^2))),decreasing =T))>0.95)[1]
six <- sv$ix[1:nv]

n.fats=nv #nrow(RLR.means.new)


##### Real simulation --------------

n.preds=5
spreds <- clo(abs(mvrnorm(n.preds,RLR.means,diag(RLR.sd)^2)))

preds <- (data.matrix(alr(diets [,six])))

mprey <-  data.matrix(clo(t(RLR.means.new[six,])))
preym <-  alr(mprey)[,]

m.fats = (n.fats-1)
# get sums of squares for each prey matrix
RLR.sim[RLR.sim==0]=0.001

ni <- rep(NA,n.preys)
R <- array(,c(m.fats,m.fats,max(RLR.sim[,1])))
for (i in 1:(max(RLR.sim[,1]))){
  ni[i] = max(n.fats,nsamples[i])
  R[,,i]=diag(diag(cov(alr(RLR.sim[RLR.sim[,1]==i,2:(n.fats+1)]))*(as.numeric(ni[i]))),m.fats)#diag(1,m.fats)#
  #R[,,i]=
}


# set uninformative prior SS matrix for wishart prior alr transformed predator data
Rnot =diag(1,m.fats)    

mean_c = matrix(1,n.preys,n.fats)
tau_coeffs = matrix(1,n.preys,n.fats)

fc_mean = as.numeric(fc.mean.new)
fc_tau = as.numeric(fc.sd.new)
fc_tau[fc_tau==0] = mean(fc_tau)
fc_tau=1/fc_tau
  
initialss=list(list(
  fcc = fc_mean,
  fracs=matrix(mean_c,n.preys,n.fats),             
  #size_c_eff = matrix(0,n.preys,2),
  ps = rep(1/n.preys,n.preys),#matrix(1/n.preys,n.preds,n.preys),
  #pmean=rep(0,n.preys),
  #A = matrix(1/n.preys,n.preys,2),#
  #pprec=diag(1,n.preys), # precision on proportion covariance
  prey.means=preym,
  predprec = diag(1,m.fats),
  prey.precs = array(1,c(n.preys,m.fats,m.fats))
))



############# don't try this unless you chose to simualte relatively few FA (else it will take literally forever)

#bugsInits(initialss,1,'FattyInits.R')
#                 
#bugsData(c('R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','ni','preds','preym'),'SQdata.r')
datas=list('R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','ni','preds','preym')

vars = c('prop','fcc')

# compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...

Fatty_Lobster_res2 <- BRugsFit('Fatty_dir_Single_prop_fc.R', datas, inits=initialss, numChains = 1, vars,
                         nBurnin = 5000, nIter = 10000, nThin = 1, coda = T,
                         DIC = F, working.directory = getwd(), digits = 4, seed=NULL,
                         BRugsVerbose = getOption("BRugsVerbose"))


Fattyl = as.data.frame((Fatty_Lobster_res2)[[1]])
names(Fattyl)

prop.ix <- grep('prop',names(Fattyl))
hist(Fattyl[,prop.ix[1]],2)
plot(Fattyl[,prop.ix[1]])
post.props = matrix(colMeans(Fattyl[,prop.ix]),n.preys,1,byrow=T)
post.lq = matrix(apply(Fattyl[,prop.ix],2,function(x){quantile(x,0.05)}),n.preys,1,byrow=T)
post.uq = matrix(apply(Fattyl[,prop.ix],2,function(x){quantile(x,0.95)}),n.preys,1,byrow=T)

ggs_density(Fattyl[,prop.ix],colMeans(Fattyl[,prop.ix]))

fc.ix <- grep('fcc',names(Fattyl))
post.fx = matrix(colMeans(Fattyl[,fc.ix]),n.preys,1,byrow=T)

pprops <- as.data.frame(Fattyl[,prop.ix])
names(pprops) <- c('H.rubra','T.undulatus','Sea Urchins','Ascidians','Brown Algae','Red Algae')
plot(pprops)

multiplot(pp1,pp3,pp5,pp2,pp4,pp6,cols=2)

####### same with SI -----

# try with Lobster Reserve data -----------

RLR.SI <- read.csv("/home/philbobsqp/Work/Dropbox/Projects private files/Fatty Acid model/Lobster diet/Reserve_SI.csv",header=F)
RLR.SI <- RLR.SI[c(2,4),]
# SQuid means and SD
RLR.SI.means <- RLR.SI[,1]  
RLR.SI.sd <- RLR.SI[,2]*sqrt(RLR[1,1])

# Prey means and SD
RLR.SI.prey.means <- RLR.SI[,seq(3,ncol(RLR.SI),2)]  
RLR.SI.prey.sd <- t(apply(data.matrix(RLR.SI[,seq(4,ncol(RLR.SI),2)]),1,function(x){x*nsamples}))

isos=2
n_SI <- nsamples
##### Real simulation --------------

# simulate 130 predators from each size class
n.preds=5 

#preds <- t(data.matrix(alr(RLR.means)))
preds.SI <- mvrnorm(n.preds,RLR.SI.means,diag(RLR.SI.sd)^2)#t(data.matrix(RLR.SI.means))

mprey <-  data.matrix(clo(t(RLR.means.new)))
preym <-  alr(mprey)[,]
preym.SI <- t(RLR.SI.prey.means )


R_SI <- array(,c(isos,isos,n.preys))
for (i in 1:n.preys){  
  R_SI[,,i] <- (diag(1/RLR.SI.prey.sd[,i])^2)*nsamples[i]
}

## first some data and inits ----

Rnot_SI = diag(1,isos)

mean_cs = c(2.1,2.9)
tau_cs =1/c(0.4,0.4)^2 # 1.61

initials.SI=list(list(
  fc = fc_mean,
  fracs=matrix(mean_c,n.preys,n.fats),             
  #size_c_eff = matrix(0,n.preys,2),
  ps = rep(1/n.preys,n.preys),#matrix(1/n.preys,n.preds,n.preys),
  #pmean=rep(0,n.preys),
  #A = matrix(1/n.preys,n.preys,2),#
  #pprec=diag(1,n.preys), # precision on proportion covariance
  prey.means=preym,
  predprec = diag(1,m.fats),
  prey.precs = array(0.1,c(n.preys,m.fats,m.fats)),
  
  #isotope inits
  cs=matrix(mean_cs,n.preys,isos,byrow=T),
  prey.means_SI=preym.SI,
  predprec_SI = diag(0.01,isos),
  prey.precs_SI = array(1,c(n.preys,isos,isos))
  
))



############# don't try this unless you chose to simualte relatively few FA (else it will take literally forever)

bugsInits(initials.SI,1,'FattyInits_SI.R')
#                 
bugsData(c('R','R_SI','fc_mean','fc_tau','mean_c','tau_coeffs','mean_cs','tau_cs','Rnot','Rnot_SI','n.preys','n.preds','isos','n.fats','m.fats','n_SI','nsamples','preds','preds.SI','preym.SI','preym'),'SIdata.r')
datas.SI=list('R','R_SI','fc_mean','fc_tau','mean_c','tau_coeffs','mean_cs','tau_cs','Rnot','Rnot_SI','n.preys','n.preds','isos','n.fats','m.fats','n_SI','nsamples','preds','preds.SI','preym.SI','preym')

vars = c('prop','fc')

# compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...

Fatty_Lobster_SI_res <- BRugsFit('Fatty_and_SI_single prop.R', datas.SI, inits=initials.SI, numChains = 1, vars,
                                 nBurnin = 5000, nIter = 10000, nThin = 10, coda = T,
                                 DIC = F, working.directory = getwd(), digits = 8, seed=NULL,
                                 BRugsVerbose = getOption("BRugsVerbose"))


Fattyla = as.data.frame((Fatty_Lobster_SI_res)[[1]])
names(Fattyla)

prop.ixa <- grep('prop',names(Fattyla))
hist(Fattyla[,prop.ixa[1]],20)
plot(Fattyla[,prop.ixa[1]])
posts.props = matrix(colMeans(Fattyla[,prop.ixa]),n.preys,1,byrow=T)
post.lq = matrix(apply(Fattyla[,prop.ixa],2,function(x){quantile(x,0.05)}),n.preys,1,byrow=T)
post.uq = matrix(apply(Fattyla[,prop.ixa],2,function(x){quantile(x,0.95)}),n.preys,1,byrow=T)

pps1 <- ggs_density(cbind(Fattyla[,prop.ixa[c(1)]],Fattyl[,prop.ix[c(1)]]),c(mean(Fattyla[,prop.ixa[c(1)]]),mean(Fattyl[,prop.ix[c(1)]])))
pps2 <- ggs_density(cbind(Fattyla[,prop.ixa[c(2)]],Fattyl[,prop.ix[c(2)]]),c(mean(Fattyla[,prop.ixa[c(2)]]),mean(Fattyl[,prop.ix[c(2)]])))
pps3 <- ggs_density(cbind(Fattyla[,prop.ixa[c(3)]],Fattyl[,prop.ix[c(3)]]),c(mean(Fattyla[,prop.ixa[c(3)]]),mean(Fattyl[,prop.ix[c(3)]])))
pps4 <- ggs_density(cbind(Fattyla[,prop.ixa[c(4)]],Fattyl[,prop.ix[c(4)]]),c(mean(Fattyla[,prop.ixa[c(4)]]),mean(Fattyl[,prop.ix[c(4)]])))
pps5 <- ggs_density(cbind(Fattyla[,prop.ixa[c(5)]],Fattyl[,prop.ix[c(5)]]),c(mean(Fattyla[,prop.ixa[c(5)]]),mean(Fattyl[,prop.ix[c(5)]])))
pps6 <- ggs_density(cbind(Fattyla[,prop.ixa[c(6)]],Fattyl[,prop.ix[c(6)]]),c(mean(Fattyla[,prop.ixa[c(6)]]),mean(Fattyl[,prop.ix[c(6)]])))


multiplot(pps1,pps3,pps5,pps2,pps4,pps6,cols=2)
write.table(rbind(post.props,post.lq,post.uq),file='SQ_real_props.csv')

# try a denstrip plot


library(ggplot2)
library(denstrip)
library(lattice)

pre <- (cbind(Fattyl[,prop.ixa],Fattyla[,prop.ix]))[,c(1,7,2,8,3,9,4,10,5,11,6,12)]

pre <- as.data.frame(pre)

# apply and sort labels
labs = c('H.rubra','T.undulatus','Sea Urchins','Ascidians','Brown Algae','Red Algae')


#tiff('relative importance posteriors.tiff', pointsize = 2,res=300,width=1500,height=1500)

mp=melt(pre)

rpp = bwplot(variable~value,data=mp,xlim=c(0,1),xlab=list(label=expression("posterior diet proportions"),cex=1),
             ,panel = function(x, y) { 
               grid.segments(1,0,0,0)
               xlist <- split(x, factor(y))
               for (i in seq(along=xlist))
                 panel.denstrip(x=xlist[[i]], at=i,colmax=(i%%2)+2,mticks=mean(xlist[[i]]))
             },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T,labels=labs,at=seq(1.5,13.5,2))))
print(rpp)

# draw line at 0 across
trellis.focus("panel", 1, 1)
panel.abline(h=seq(2.5,10.5,2),col=1,lty=2)
trellis.unfocus()




############ do fished diets ------------

FLR <- read.csv("/home/philbobsqp/Work/Dropbox/Projects private files/Fatty Acid model/Lobster diet/Fished_FA.csv",header=F)

# Lobster means and SD
FLR.means <- clo(FLR[2:nrow(FLR),1])
FLR.sd <- FLR[2:nrow(FLR),2]*sqrt(FLR[1,1])/100

# Prey means and SD
FLR.prey.means <- apply(FLR[2:nrow(FLR),seq(3,ncol(FLR),2)],2,clo)
FLR.prey.sd <- t(apply(data.matrix(FLR[2:nrow(FLR),seq(4,ncol(FLR),2)]),1,function(x){x*sqrt(data.matrix(FLR[1,seq(4,ncol(FLR),2)]))}))/100

# pool data 

FLR.means.new <- cbind(FLR.prey.means[,1:2],apply(FLR.prey.means[,3:4],1,gmean),apply(FLR.prey.means[,5:7],1,gmean),apply(FLR.prey.means[,8:9],1,gmean),FLR.prey.means[,10])
# just use mean sd - numbers are comparable and it's just an illustration...
FLR.sd.new <- cbind(FLR.prey.sd[,1:2],rowMeans(FLR.prey.sd[,3:4]),rowMeans(FLR.prey.sd[,5:7]),rowMeans(FLR.prey.sd[,8:9]),FLR.prey.sd[,10])

# replace zeros in sds
for (i in 1:nrow(FLR.sd.new))
  FLR.sd.new[i,FLR.sd.new[i,]==0] <- 1e-8#min(FLR.sd.new[i,FLR.sd.new[i,]!=0])


#get fat contents
fc.mean.new = rep(1,n.preys)#c(1.4,1.3,4.8,1,0.5,0.3)
fc.sd.new = c(fc.mean.new /10)

FLR.means.new[FLR.means.new==0]=0.00001

nsample = as.numeric(FLR[1,seq(3,ncol(FLR-1),2)])
nsamples = c(nsample[1:2],sum(nsample[3:4]),sum(nsample[5:7]),sum(nsample[8:9]),nsample[10])

FLR.sim <- rep(NA,nrow(FLR.means.new)+1)
for (i in 1:ncol(FLR.means.new))
  FLR.sim <- rbind(FLR.sim,cbind(rep(i,50),abs(mvrnorm(50,FLR.means.new[,i],diag(FLR.sd.new[,i])^2))))
FLR.sim <- (FLR.sim[-1,])
FLR.sim[,2:ncol(FLR.sim)] <- clo(FLR.sim[,2:ncol(FLR.sim)])


FLR.RDA <- capscale(dista~as.factor(FLR.sim[,1]),comm=FLR.sim[,2:ncol(FLR.sim)])
#plot(FLR.RDA,t='n',xlim=c(-0.5,0.5),ylim=c(-1,1))
plot(t(apply(t(FLR.means.new)%*%FLR.RDA$CCA$v[,1:2],1,function(x){x*c(50,100)})),pch=1:6,ylim=c(-2,2.5))
points(clo(t(FLR.means))%*%FLR.RDA$CCA$v[,1:2]*c(50,100),pch=16)
legend(0.2,2.5,c('H.rubra','T.undulatus','Sea Urchins','Ascidians','Brown Algae','Red Algae'),pch=c(1:6,16,17))

plot(cumsum(sort(clo(rowSums(t(t(cbind(FLR.RDA$CCA$v)*c(FLR.RDA$CCA$eig))^2))),decreasing =T)),ylab='cumulative proportion',xlab='n')

sv = sort(clo(rowSums(t(t(cbind(FLR.RDA$CCA$v,FLR.RDA$CA$v))*c(FLR.RDA$CCA$eig,FLR.RDA$CA$eig))^2)),decreasing =T,index.return=T)
nv <- 12
six <- sv$ix[1:nv]

n.fats=nv #nrow(RLR.means.new)
n.preys=ncol(FLR.means.new)

##### Real simulation --------------

# simulate 130 predators from each size class
n.preds=10 
spreds <- clo(abs(mvrnorm(n.preds,FLR.means,diag(FLR.sd)^2)))

preds <- (data.matrix(alr(spreds[,six])))[,]

mprey <-  data.matrix(clo(t(FLR.means.new[six,])))
preym <-  alr(mprey)[,]

m.fats = (n.fats-1)
# get sums of squares for each prey matrix
FLR.sim[FLR.sim==0]=0.001

ni <- rep(NA,n.preys)
R <- array(,c(m.fats,m.fats,max(FLR.sim[,1])))
for (i in 1:(max(FLR.sim[,1]))){
  ni[i] = max(n.fats,nsamples[i])
  R[,,i]=diag(diag(cov(alr(FLR.sim[FLR.sim[,1]==i,2:(n.fats+1)]))*(as.numeric(ni[i]))),m.fats)#diag(1,m.fats)#
  #R[,,i]=
}
# set uninformative prior SS matrix for wishart prior alr transformed predator data
Rnot =diag(1,m.fats)    

mean_c = matrix(1,n.preys,n.fats)
tau_coeffs = matrix(1,n.preys,n.fats)

fc_mean = as.numeric(fc.mean.new)
fc_tau = as.numeric(fc.sd.new^2)
fc_tau[fc_tau==0] = mean(fc_tau)
fc_tau=1/fc_tau


initialss=list(list(
  fc = fc_mean,
  fracs=matrix(mean_c,n.preys,n.fats),             
  #size_c_eff = matrix(0,n.preys,2),
  ps = rep(1/n.preys,n.preys),#matrix(1/n.preys,n.preds,n.preys),
  #pmean=rep(0,n.preys),
  #A = matrix(1/n.preys,n.preys,2),#
  #pprec=diag(1,n.preys), # precision on proportion covariance
  prey.means=preym,
  predprec = diag(1,m.fats),
  prey.precs = array(0.1,c(n.preys,m.fats,m.fats))
))



############# don't try this unless you chose to simualte relatively few FA (else it will take literally forever)

#bugsInits(initialss,1,'FattyInits.R')
#                 
#bugsData(c('R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','nsamples','preds','preym'),'SQdata.r')
#datas=list('R','fc_mean','fc_tau','mean_c','tau_coeffs','Rnot','n.preys','n.preds','n.fats','m.fats','nsamples','preds','preym')

vars = c('prop','fcc')

# compilation time and return time once OpenBUGS has finished can be very long. Patience is of the essence...

Fatty_Lobster_Fishw <- BRugsFit('Fatty_dir_Single_prop_fc.R', datas, inits=initialss, numChains = 1, vars,
                          nBurnin = 10000, nIter = 50000, nThin = 50, coda = T,
                          DIC = F, working.directory = getwd(), digits = 4, seed=NULL,
                          BRugsVerbose = getOption("BRugsVerbose"))


Fattyls = as.data.frame((Fatty_Lobster_Fishw)[[1]])
names(Fattyls)

prop.ix <- grep('prop',names(Fattyls))
hist(Fattyls[,prop.ix[1]],20)
plot(Fattyls[,prop.ix[1]])
postas.propas = matrix(colMeans(Fattyls[,prop.ix]),n.preys,1,byrow=T)
post.lq = matrix(apply(Fattyls[,prop.ix],2,function(x){quantile(x,0.05)}),n.preys,1,byrow=T)
post.uq = matrix(apply(Fattyls[,prop.ix],2,function(x){quantile(x,0.95)}),n.preys,1,byrow=T)

fc.ix <- grep('fcc',names(Fattyls))
post.fx = matrix(colMeans(Fattyls[,fc.ix]),n.preys,1,byrow=T)


pp1 <- ggs_density(Fattylsss[,prop.ix[c(1)]],mean(Fattylsss[,prop.ix[c(1)]]))
pp2 <- ggs_density(Fattylsss[,prop.ix[c(2)]],mean(Fattylsss[,prop.ix[c(2)]]))
pp3 <- ggs_density(Fattylsss[,prop.ix[c(3)]],mean(Fattylsss[,prop.ix[c(3)]]))
pp4 <- ggs_density(Fattylsss[,prop.ix[c(4)]],mean(Fattylsss[,prop.ix[c(4)]]))
pp5 <- ggs_density(Fattylsss[,prop.ix[c(5)]],mean(Fattylsss[,prop.ix[c(5)]]))
pp6 <- ggs_density(Fattylsss[,prop.ix[c(6)]],mean(Fattylsss[,prop.ix[c(6)]]))


multiplot(pp1,pp3,pp5,pp2,pp4,pp6,cols=2)
write.table(rbind(post.props,post.lq,post.uq),file='SQ_real_props.csv')




############ SI in fished -----

FLR.SI <- read.csv("/home/philbobsqp/Work/Dropbox/Projects private files/Fatty Acid model/Lobster diet/Fished_SI.csv",header=F)
FLR.SI <- FLR.SI[c(2,4),]
# SQuid means and SD
FLR.SI.means <- FLR.SI[,1]  
FLR.SI.sd <- FLR.SI[,2]*sqrt(FLR[1,1])

# Prey means and SD
FLR.SI.prey.means <- FLR.SI[,seq(3,ncol(FLR.SI),2)]  
FLR.SI.prey.sd <- t(apply(data.matrix(FLR.SI[,seq(4,ncol(FLR.SI),2)]),1,function(x){x*nsamples}))

isos=2
n_SI <- nsamples
##### Real simulation --------------

# simulate 130 predators from each size class
n.preds=10 

#preds <- #t(data.matrix(alr(FLR.means)))
preds.SI <- mvrnorm(n.preds,FLR.SI.means,diag(FLR.SI.sd)^2)#t(data.matrix(FLR.SI.means))

mprey <-  data.matrix(clo(t(FLR.means.new)))
preym <-  alr(mprey[,six])[,]
preym.SI <- t(FLR.SI.prey.means )


R_SI <- array(,c(isos,isos,n.preys))
for (i in 1:n.preys){  
  R_SI[,,i] <- (diag(1/FLR.SI.prey.sd[,i])^2)*nsamples[i]
}

## first some data and inits ----

Rnot_SI = diag(0.01,isos)

mean_cs = c(2.1,2.9)
tau_cs =1/c(0.4,0.4)^2 # 1.61

initials.SI=list(list(
  fc = fc_mean,
  fracs=matrix(mean_c,n.preys,n.fats),             
  #size_c_eff = matrix(0,n.preys,2),
  ps = rep(1/n.preys,n.preys),#matrix(1/n.preys,n.preds,n.preys),
  #pmean=rep(0,n.preys),
  #A = matrix(1/n.preys,n.preys,2),#
  #pprec=diag(1,n.preys), # precision on proportion covariance
  prey.means=preym,
  predprec = diag(1,m.fats),
  prey.precs = array(0.1,c(n.preys,m.fats,m.fats)),
  
  #isotope inits
  cs=mean_cs,
  prey.means_SI=preym.SI,
  predprec_SI = diag(0.01,isos),
  prey.precs_SI = array(1,c(n.preys,isos,isos))
  
))



############# don't try this unless you chose to simualte relatively few FA (else it will take literally forever)

 bugsInits(initials.SI,1,'FattyInits_SI.R')
#                 
 bugsData(c('R','R_SI','fc_mean','fc_tau','mean_c','tau_coeffs','mean_cs','tau_cs','Rnot','Rnot_SI','n.preys','n.preds','isos','n.fats','m.fats','n_SI','nsamples','preds','preds.SI','preym.SI','preym'),'SIdata.r')
datas.SI=list('R','R_SI','fc_mean','fc_tau','mean_c','tau_coeffs','mean_cs','tau_cs','Rnot','Rnot_SI','n.preys','n.preds','isos','n.fats','m.fats','n_SI','nsamples','preds','preds.SI','preym.SI','preym')

vars = c('prop')

# compilation time and return time oncj,e OpenBUGS has finished can be very long. Patience is of the essence...

Fatty_Lobster_SI_fished <- BRugsFit('Fatty_and_SI_single prop.R', datas.SI, inits=initials.SI, numChains = 1, vars,
                                 nBurnin = 10000, nIter = 50000, nThin = 50, coda = T,
                                 DIC = F, working.directory = getwd(), digits = 8, seed=NULL,
                                 BRugsVerbose = getOption("BRugsVerbose"))


Fattylas = as.data.frame((Fatty_Lobster_SI_fished)[[1]])
names(Fattylas)

prop.ixa <- grep('prop',names(Fattylas))
ggs_density(Fattylas[,prop.ixa],colMeans(Fattylas[,prop.ixa]))

fc.ix <- grep('fc',names(Fattylass))
post.fx = matrix(colMeans(Fattyla[,fc.ix]),n.preys,1,byrow=T)
ggs_density(Fattyla[,fc.ix],colMeans(Fattyla[,fc.ix]))


pprops <- as.data.frame(clr(Fattyla[,prop.ix]))
names(pprops) <- c('H.rubra','T.undulatus','Sea Urchins','Ascidians','Brown Algae','Red Algae')
plot(pprops)

Fattylass = as.data.frame((Fatty_Lobster_SI_fished2)[[1]])
fc.ix <- grep('cs',names(Fattylass))
colMeans(Fattylass[,fc.ix])

hist(Fattylas[,prop.ixa[1]],20)
plot(Fattylas[,prop.ixa[1]])
post.propas = matrix(colMeans(Fattylas[,prop.ixa]),n.preys,1,byrow=T)
post.lq = matrix(apply(Fattylas[,prop.ixa],2,function(x){quantile(x,0.05)}),n.preys,1,byrow=T)
post.uq = matrix(apply(Fattylas[,prop.ixa],2,function(x){quantile(x,0.95)}),n.preys,1,byrow=T)

ppss1 <- ggs_density(cbind(Fattylas[,prop.ixa[c(1)]],Fattyls[,prop.ix[c(1)]]),c(mean(Fattylas[,prop.ixa[c(1)]]),mean(Fattyls[,prop.ix[c(1)]])))
ppss2 <- ggs_density(cbind(Fattylas[,prop.ixa[c(2)]],Fattyls[,prop.ix[c(2)]]),c(mean(Fattylas[,prop.ixa[c(2)]]),mean(Fattyls[,prop.ix[c(2)]])))
ppss3 <- ggs_density(cbind(Fattylas[,prop.ixa[c(3)]],Fattyls[,prop.ix[c(3)]]),c(mean(Fattylas[,prop.ixa[c(3)]]),mean(Fattyls[,prop.ix[c(3)]])))
ppss4 <- ggs_density(cbind(Fattylas[,prop.ixa[c(4)]],Fattyls[,prop.ix[c(4)]]),c(mean(Fattylas[,prop.ixa[c(4)]]),mean(Fattyls[,prop.ix[c(4)]])))
ppss5 <- ggs_density(cbind(Fattylas[,prop.ixa[c(5)]],Fattyls[,prop.ix[c(5)]]),c(mean(Fattylas[,prop.ixa[c(5)]]),mean(Fattyls[,prop.ix[c(5)]])))
ppss6 <- ggs_density(cbind(Fattylas[,prop.ixa[c(6)]],Fattyls[,prop.ix[c(6)]]),c(mean(Fattylas[,prop.ixa[c(6)]]),mean(Fattyls[,prop.ix[c(6)]])))


multiplot(ppss1,ppss3,ppss5,ppss2,ppss4,ppss6,cols=2)
write.table(rbind(post.props,post.lq,post.uq),file='SQ_real_props.csv')

# try a denstrip plot


library(ggplot2)
library(denstrip)
library(lattice)

pre <- (cbind(Fattyls[,prop.ix],Fattylas[,prop.ixa]))[,c(1,7,2,8,3,9,4,10,5,11,6,12)]

pre <- as.data.frame(pre)

# apply and sort labels
labs = c('H.rubra','T.undulatus','Sea Urchins','Ascidians','Brown Algae','Red Algae')


#tiff('relative importance posteriors.tiff', pointsize = 2,res=300,width=1500,height=1500)

mp=melt(pre)

rpp = bwplot(variable~value,data=mp,xlim=c(0,1),xlab=list(label=expression("posterior diet proportions"),cex=1),
             ,panel = function(x, y) { 
               grid.segments(1,0,0,0)
               xlist <- split(x, factor(y))
               for (i in seq(along=xlist))
                 panel.denstrip(x=xlist[[i]], at=i,colmax=(i%%2)+2,mticks=mean(xlist[[i]]))
             },par.settings = list(axis.line = list(col=NA)),scales=list(col=1,cex=1,x=list(col=1,at=seq(0,1,0.2)),y=list(draw=T,labels=labs,at=seq(1.5,13.5,2))))
print(rpp)

# draw line at 0 across
trellis.focus("panel", 1, 1)
panel.abline(h=seq(2.5,10.5,2),col=1,lty=2)
trellis.unfocus()


