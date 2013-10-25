# simulate fatty acid proportions in 3 prey species ----
require(FASTIN)

# Lobster Reserve data -----------

kappas <- function(mat,vars){
        ks <- vector(,length(vars))
        for (i in 1:length(vars)){
            ks[i] <- kappa(mat[,vars[1:i]],exact=T)
        }
        return(ks)
    }
    
clo <- function(x){
      
        if (is.null(dim(x))){xc <- x/sum(x)} else {xc <- t(apply(x,1,function(y){y/sum(y)}))}
        return(xc)
    }

alr <- function(x){
        x<-clo(x)
        if (is.null(dim(x))){xc <- log(x[1:(length(x)-1)]/x[length(x)])} else { t(apply(x,1,function(y){log(y[1:(length(y)-1)]/y[length(y)])}))}
    }

adist <- function(mat){

        dims <- dim(mat)
        dists <- matrix(,dims[1],dims[1])
        for (i in 1:(dims[1]-1)){
            for (j in (i+1):dims[1]){
                dists[j,i] <- robCompositions::aDist(mat[i,],mat[j,])
            }}
        dista <- as.dist(dists)
              return(dista)
    }

setwd("../Lobster/")

RLR <- read.csv("Reserve_FA.csv",header=F)

# Lobster means and SD
RLR.means <- clo(RLR[2:nrow(RLR),1])
RLR.sd <- RLR[2:nrow(RLR),2]*sqrt(RLR[1,1])/100

# Prey means and SD
RLR.prey.means <- apply(RLR[2:nrow(RLR),seq(3,ncol(RLR),2)],2,clo)
RLR.prey.sd <- t(apply(data.matrix(RLR[2:nrow(RLR),seq(4,ncol(RLR),2)]),1,function(x){x*sqrt(data.matrix(RLR[1,seq(4,ncol(RLR),2)]))}))/100

###########################################################
###### Data grooming: combining sources, ##################
###### and a few (minor) assumptions...  ##################
###########################################################

RLR.means.new <- t(cbind(RLR.prey.means[,1:2],apply(RLR.prey.means[,3:4],1,function(x){exp(mean(log(x)))}),apply(RLR.prey.means[,5:7],1,function(x){exp(mean(log(x)))}),apply(RLR.prey.means[,8:9],1,function(x){exp(mean(log(x)))}),RLR.prey.means[,10]))
# just use mean sd - numbers are comparable and it's just an illustration...
RLR.sd.new <- t(cbind(RLR.prey.sd[,1:2],rowMeans(RLR.prey.sd[,3:4]),rowMeans(RLR.prey.sd[,5:7]),rowMeans(RLR.prey.sd[,8:9]),RLR.prey.sd[,10]))


for (i in 1:ncol(RLR.means.new))
 RLR.means.new[RLR.means.new[,i]==0,i] <- min(RLR.means.new[RLR.means.new[,i]!=0,i])/10
# replace zeros in sds
for (i in 1:ncol(RLR.sd.new))
  RLR.sd.new[RLR.sd.new[,i]==0,i] <- min(RLR.sd.new[RLR.sd.new[,i]!=0,i])

n.preys=ncol(RLR.means.new)
n.fats=ncol(RLR.means.new)

#get fat contents -set to one as lit values are jsut as unreliable
fc.mean.new = rep(1,n.preys)# c(1.4,1.3,4.8,1,0.5,0.3)
fc.sd.new = c(fc.mean.new /1)

nsample = as.numeric(RLR[1,seq(3,ncol(RLR-1),2)])
nsamples = c(nsample[1:2],sum(nsample[3:4]),sum(nsample[5:7]),sum(nsample[8:9]),nsample[10])

dista <- adist(RLR.means.new)

RLR.RDA <- vegan::capscale(dista~as.factor(1:n.preys),comm=RLR.means.new)
#plot(RLR.RDA,t='n',xlim=c(-0.5,0.5),ylim=c(-1,1))
plot(RLR.means.new%*%RLR.RDA$CCA$v[,1:2],pch=1:6)
points(clo(RLR.means)%*%RLR.RDA$CCA$v[,1:2],pch=16)
legend('topleft',c('H.rubra','T.undulatus','Sea Urchins','Ascidians','Brown Algae','Red Algae','Reserve Lobster'),pch=c(1:6,16,17))


######################################
######## variable selection ##########
######################################

n.fats=ncol(RLR.means.new)
par(mfcol=c(2,1))
sv = sort(clo(rowSums(t(t(RLR.RDA$CCA$v)*RLR.RDA$CCA$eig)^2)),decreasing =T,index.return=T)
plot(cumsum(sort(clo(rowSums(t(t(RLR.RDA$CCA$v)*RLR.RDA$CCA$eig)^2)),decreasing =T)),axes=F,xlab='',ylab='Cumulative source separation',ylim=c(0,1))
axis(2)
axis(1,at=1:n.fats,labels=F)
text(1:n.fats, par("usr")[3] - 0.2, srt = 45, adj = 1,
     labels = colnames(RLR.means.new)[sv$ix], xpd = TRUE)

ks <- kappas(RLR.means.new,sv$ix)
plot(ks,axes=F,ylab='Prey matrix condition number',xlab='',ylim=c(0,max(ks)))
axis(2)
axis(1,at=1:n.fats,labels=F)
text(1:n.fats, par("usr")[3] - max(ks)/5 , srt = 45, adj = 1,
     labels = colnames(RLR.means.new)[sv$ix], xpd = TRUE)
   
cumsum(sort(clo(rowSums(t(t(RLR.RDA$CCA$v)*RLR.RDA$CCA$eig)^2)),decreasing =T))

 nv <- select.list(title='please choose the fatty acids to use (at least 3)',choices = colnames(RLR.means.new)[sv$ix],graphics=T,multiple=T)
      
      #nv <- readline(prompt = "please enter number of variables for analysis \n")
six <- match(nv,colnames(RLR.means.new))

n.fats=length(six)


##### Real simulation --------------

n.preds=1
#spreds <- clo(abs(mvrnorm(n.preds,RLR.means,diag(RLR.sd)^2)))
spreds <- RLR.means
preds <- (data.matrix(alr(spreds[six])))

mprey <-  data.matrix(clo(RLR.means.new[,six]))
preym <-  alr(mprey)[,]

m.fats = (n.fats-1)

# get sums of squares for each prey matrix
#number of simulated prey samples per species/group
ns <- 300
RLR.sim <- array(,c(ns,n.fats,n.preys))

for (i in 1:n.preys) RLR.sim[,,i] <- mvrnorm(ns,RLR.means.new[i,six],diag(RLR.sd.new[i,six]^2))

for (i in (1:n.preys))
RLR.sim[,,i] <- apply(RLR.sim[,,i],2,function(x){x[x<=0] = min(x[x>0])/10;return(x)})

R <- array(,c(m.fats,m.fats,n.preys))
ni<-rep(NA,n.preys)
for (i in 1:n.preys){
    ni[i] <- max(n.fats,nsamples[i]-1)
    R[,,i]=cov(alr(RLR.sim[,,i]))*ni[i]
}

# set uninformative prior SS matrix for wishart prior alr transformed predator data
Rnot =diag(1,m.fats)    

# # set SS matrix for wishart prior on proportions
even=0.05
2*(1-pnorm(log(95)/2,0,sqrt(1/even))) # chance that 95*p1<p2

mean_c = matrix(1,n.preys,n.fats)
tau_coeffs = matrix(1,n.preys,n.fats)

fc_mean = as.numeric(fc.mean.new)
fc_tau = as.numeric(fc.sd.new)
fc_tau[fc_tau==0] = mean(fc_tau)
fc_tau=1/fc_tau
  
#############################################################
########### make data object and run analysis ###############
#############################################################

# Fatty Acid data (stable isotopes are intigrated in the same way

datas.FA <- list(n.fats = n.fats,m.fats=m.fats,fc_mean=fc_mean,fc_tau =fc_tau,R=R,Rnot = Rnot,preym=preym,preds = preds,ni=ni,mean_c=mean_c,tau_c = tau_coeffs)

RLR.data <- list(datas.FA=datas.FA,n.preys=n.preys,n.preds=n.preds,even=even,prey.ix=c('H.rubra','T.undulatus','Sea Urchins','Ascidians','Brown Algae','Red Algae'))

guiSet('datas',RLR.data)

# MCMC defaults to population proportions only from FA data,
RLR.out <- run_MCMC(nIter=100000,nBurnin=10000,nChains=3,nThin=100,datas = RLR.data)

plot(RLR.out)

summary(RLR.out)


###############################
####### same with SI -----
###############################

# try with Lobster Reserve data -----------

RLR.SI <- read.csv("Reserve_SI.csv",header=F)
RLR.SI <- RLR.SI[c(2,4),]
# SQuid means and SD
RLR.SI.means <- RLR.SI[,1]  
RLR.SI.sd <- RLR.SI[,2]*sqrt(RLR[1,1])

# Prey means and SD
RLR.SI.prey.means <- t(RLR.SI[,seq(3,ncol(RLR.SI),2)])
RLR.SI.prey.sd <- apply(data.matrix(RLR.SI[,seq(4,ncol(RLR.SI),2)]),1,function(x){x*nsamples})

isos=2
n_SI <- nsamples
##### Real simulation --------------

# simulate 130 predators from each size class
n.preds=1

preds.SI <- RLR.SI.means#mvrnorm(n.preds,RLR.SI.means,diag(RLR.SI.sd)^2)#t(data.matrix(RLR.SI.means))

preym.SI <- RLR.SI.prey.means


R_SI <- array(,c(isos,isos,n.preys))
ni.SI<-rep(NA,n.preys)
for (i in 1:n.preys){
    ni.SI[i] <- max(isos,nsamples[i]-1)
  R_SI[,,i] <- (diag(1/RLR.SI.prey.sd[i,])^2)*ni.SI[i]
}

## first some data and inits ----

Rnot_SI = diag(0.1,isos)

mean_cs = c(2.1,2.9)
mean_cs <- matrix(unlist(mean_cs), n.preys, isos,byrow=T)
tau_cs =1/c(0.4,0.4)^2 # 1.61
tau_cs <- matrix(unlist(tau_cs), n.preys, isos,byrow=T)

datas.SI <- list(isos = isos,R.SI=R_SI,Rnot.SI = Rnot_SI,preym.SI=preym.SI,preds.SI = preds.SI,ni.SI=ni.SI,mean_cs=mean_cs,tau_cs = tau_cs)

RLR.data.full <- list(datas.FA=datas.FA,datas.SI=datas.SI,n.preys=n.preys,n.preds=n.preds,even=even,prey.ix=c('H.rubra','T.undulatus','Sea Urchins','Ascidians','Brown Algae','Red Algae'))

guiSet('datas',RLR.data.full)

# MCMC defaults to population proportions only from FA data, so need to put SI data in explicitly
RLR.out.SI <- run_MCMC(nIter=100000,nBurnin=10000,nChains=3,nThin=10,datas = RLR.data.full,Data.Type='Stable.Isotopes')

plot(RLR.out.SI)

# Do the combined analysis
RLR.out.full <- run_MCMC(nIter=100000,nBurnin=10000,nChains=3,nThin=10,datas = RLR.data.full,Data.Type='Combined.Analysis')

plot(RLR.out.full)


##########################
## show all together #####
##########################

reserve.list <- list('Reserve Full'=RLR.out.full,'Reserve SI'=RLR.out.SI,'Reserve FA'=RLR.out)

multiplot(reserve.list,density=F)
multiplot(fished.list,density=T)







###########################################
############ do fished diets ------------##
###########################################

FLR <- read.csv("Fished_FA.csv",header=F)

# Lobster means and SD
FLR.means <- clo(FLR[2:nrow(FLR),1])
FLR.sd <- FLR[2:nrow(FLR),2]*sqrt(FLR[1,1])/100

# Prey means and SD
FLR.prey.means <- apply(FLR[2:nrow(FLR),seq(3,ncol(FLR),2)],2,clo)
FLR.prey.sd <- t(apply(data.matrix(FLR[2:nrow(FLR),seq(4,ncol(FLR),2)]),1,function(x){x*sqrt(data.matrix(FLR[1,seq(4,ncol(FLR),2)]))}))/100

###########################################################
###### Data grooming: combining sources, ##################
###### and a few (minor) assumptions...  ##################
###########################################################

FLR.means.new <- t(cbind(FLR.prey.means[,1:2],apply(FLR.prey.means[,3:4],1,function(x){exp(mean(log(x)))}),apply(FLR.prey.means[,5:7],1,function(x){exp(mean(log(x)))}),apply(FLR.prey.means[,8:9],1,function(x){exp(mean(log(x)))}),FLR.prey.means[,10]))
# just use mean sd - numbers are comparable and it's just an illustration...
FLR.sd.new <- t(cbind(FLR.prey.sd[,1:2],rowMeans(FLR.prey.sd[,3:4]),rowMeans(FLR.prey.sd[,5:7]),rowMeans(FLR.prey.sd[,8:9]),FLR.prey.sd[,10]))


for (i in 1:ncol(FLR.means.new))
 FLR.means.new[FLR.means.new[,i]==0,i] <- min(FLR.means.new[FLR.means.new[,i]!=0,i])/10
# replace zeros in sds
for (i in 1:ncol(FLR.sd.new))
  FLR.sd.new[FLR.sd.new[,i]==0,i] <- min(FLR.sd.new[FLR.sd.new[,i]!=0,i])

n.preys=nrow(FLR.means.new)
n.fats=ncol(FLR.means.new)

#get fat contents -set to one as lit values are jsut as unreliable
fc.mean.new = rep(1,n.preys)# c(1.4,1.3,4.8,1,0.5,0.3)
fc.sd.new = c(fc.mean.new /1)

nsample = as.numeric(FLR[1,seq(3,ncol(FLR-1),2)])
nsamples = c(nsample[1:2],sum(nsample[3:4]),sum(nsample[5:7]),sum(nsample[8:9]),nsample[10])

dista <- adist(FLR.means.new)

FLR.RDA <- vegan::capscale(dista~as.factor(1:n.preys),comm=FLR.means.new)
#plot(FLR.RDA,t='n',xlim=c(-0.5,0.5),ylim=c(-1,1))
plot(FLR.means.new%*%FLR.RDA$CCA$v[,1:2],pch=1:6)
points(clo(FLR.means)%*%FLR.RDA$CCA$v[,1:2],pch=16)
legend('topleft',c('H.rubra','T.undulatus','Sea Urchins','Ascidians','Brown Algae','Red Algae','Reserve Lobster'),pch=c(1:6,16,17))


######################################
######## variable selection ##########
######################################

n.fats=ncol(FLR.means.new)
par(mfcol=c(2,1))
sv = sort(clo(rowSums(t(t(FLR.RDA$CCA$v)*FLR.RDA$CCA$eig)^2)),decreasing =T,index.return=T)
plot(cumsum(sort(clo(rowSums(t(t(FLR.RDA$CCA$v)*FLR.RDA$CCA$eig)^2)),decreasing =T)),axes=F,xlab='',ylab='Cumulative source separation',ylim=c(0,1))
axis(2)
axis(1,at=1:n.fats,labels=F)
text(1:n.fats, par("usr")[3] - 0.2, srt = 45, adj = 1,
     labels = colnames(FLR.means.new)[sv$ix], xpd = TRUE)

ks <- kappas(FLR.means.new,sv$ix)
plot(ks,axes=F,ylab='Prey matrix condition number',xlab='',ylim=c(0,max(ks)))
axis(2)
axis(1,at=1:n.fats,labels=F)
text(1:n.fats, par("usr")[3] - max(ks)/5 , srt = 45, adj = 1,
     labels = colnames(FLR.means.new)[sv$ix], xpd = TRUE)
   

 nv <- select.list(title='please choose the fatty acids to use (at least 3)',choices = colnames(FLR.means.new)[sv$ix],graphics=T,multiple=T)
      
      #nv <- readline(prompt = "please enter number of variables for analysis \n")
six <- match(nv,colnames(FLR.means.new))

n.fats=length(six)


##### Real simulation --------------

n.preds=1
#spreds <- clo(abs(mvrnorm(n.preds,FLR.means,diag(FLR.sd)^2)))
spreds <- FLR.means
preds <- (data.matrix(alr(spreds[six])))

mprey <-  data.matrix(clo(FLR.means.new[,six]))
preym <-  alr(mprey)[,]

m.fats = (n.fats-1)

# get sums of squares for each prey matrix
#number of simulated prey samples per species/group
ns <- 30
FLR.sim <- array(,c(ns,n.fats,n.preys))

for (i in 1:n.preys) FLR.sim[,,i] <- mvrnorm(ns,FLR.means.new[i,six],diag(FLR.sd.new[i,six]^2))

for (i in (1:n.preys))
FLR.sim[,,i] <- apply(FLR.sim[,,i],2,function(x){x[x<=0] = min(x[x>0])/10;return(x)})

R <- array(,c(m.fats,m.fats,n.preys))
ni<-rep(NA,n.preys)
for (i in 1:n.preys){
    ni[i] <- max(n.fats,nsamples[i]-1)
    R[,,i]=cov(alr(FLR.sim[,,i]))*ni[i]
}

# set uninformative prior SS matrix for wishart prior alr transformed predator data
Rnot =diag(1,m.fats)    

# # set SS matrix for wishart prior on proportions
even=0.05
2*(1-pnorm(log(95)/2,0,sqrt(1/even))) # chance that 95*p1<p2

mean_c = matrix(1,n.preys,n.fats)
tau_coeffs = matrix(1,n.preys,n.fats)

fc_mean = as.numeric(fc.mean.new)
fc_tau = as.numeric(fc.sd.new)
fc_tau[fc_tau==0] = mean(fc_tau)
fc_tau=1/fc_tau
  
#############################################################
########### make data object and run analysis ###############
#############################################################

# Fatty Acid data (stable isotopes are intigrated in the same way

datas.FA <- list(n.fats = n.fats,m.fats=m.fats,fc_mean=fc_mean,fc_tau =fc_tau,R=R,Rnot = Rnot,preym=preym,preds = preds,ni=ni,mean_c=mean_c,tau_c = tau_coeffs)

FLR.data <- list(datas.FA=datas.FA,n.preys=n.preys,n.preds=n.preds,even=even,prey.ix=c('H.rubra','T.undulatus','Sea Urchins','Ascidians','Brown Algae','Red Algae'))

guiSet('datas',FLR.data)

# MCMC defaults to population proportions only from FA data,
FLR.out <- run_MCMC(nIter=100000,nBurnin=10000,nChains=3,nThin=100,datas = FLR.data)

plot(FLR.out)

summary(FLR.out)


###############################
####### same with SI -----
###############################

# try with Lobster Reserve data -----------

FLR.SI <- read.csv("Fished_SI.csv",header=F)
FLR.SI <- FLR.SI[c(2,4),]
# SQuid means and SD
FLR.SI.means <- FLR.SI[,1]  
FLR.SI.sd <- FLR.SI[,2]*sqrt(FLR[1,1])

# Prey means and SD
FLR.SI.prey.means <- t(FLR.SI[,seq(3,ncol(FLR.SI),2)])
FLR.SI.prey.sd <- apply(data.matrix(FLR.SI[,seq(4,ncol(FLR.SI),2)]),1,function(x){x*nsamples})

isos=2
n_SI <- nsamples
##### Real simulation --------------

# simulate 130 predators from each size class
n.preds=1

preds.SI <- FLR.SI.means#mvrnorm(n.preds,FLR.SI.means,diag(FLR.SI.sd)^2)#t(data.matrix(FLR.SI.means))

preym.SI <- FLR.SI.prey.means


R_SI <- array(,c(isos,isos,n.preys))
ni.SI<-rep(NA,n.preys)
for (i in 1:n.preys){
    ni.SI[i] <- max(isos,nsamples[i]-1)
  R_SI[,,i] <- (diag(1/FLR.SI.prey.sd[i,])^2)*ni.SI[i]
}

## first some data and inits ----

Rnot_SI = diag(1,isos)

mean_cs = c(2.1,2.9)
mean_cs <- matrix(unlist(mean_cs), n.preys, isos,byrow=T)
tau_cs =1/c(0.4,0.4)^2 # 1.61
tau_cs <- matrix(unlist(tau_cs), n.preys, isos,byrow=T)

datas.SI <- list(isos = isos,R.SI=R_SI,Rnot.SI = Rnot_SI,preym.SI=preym.SI,preds.SI = preds.SI,ni.SI=ni.SI,mean_cs=mean_cs,tau_cs = tau_cs)

FLR.data.full <- list(datas.FA=datas.FA,datas.SI=datas.SI,n.preys=n.preys,n.preds=n.preds,even=even,prey.ix=c('H.rubra','T.undulatus','Sea Urchins','Ascidians','Brown Algae','Red Algae'))

guiSet('datas',FLR.data.full)

# MCMC defaults to population proportions only from FA data, so need to put SI data in explicitly
FLR.out.SI <- run_MCMC(nIter=100000,nBurnin=10000,nChains=3,nThin=100,datas = FLR.data.full,Data.Type='Stable.Isotopes')

plot(FLR.out.SI)

# Do the combined analysis
FLR.out.full <- run_MCMC(nIter=100000,nBurnin=10000,nChains=3,nThin=100,datas = FLR.data.full,Data.Type='Combined.Analysis')

plot(FLR.out.full)

##########################
## show all together #####
##########################

fished.list <- list('Fished Full'=FLR.out.full,'Fished SI'=FLR.out.SI,'Fished FA'=FLR.out)

multiplot(fished.list,density=F)
multiplot(fished.list,density=T)


##########################################################
###### Fished vs reserve from full dataset (FA+SI) #######
##########################################################

Lobster.list <- list('Fished'=FLR.out.full,'Reserve'=RLR.out.full)

multiplot(Lobster.list,density=F)
multiplot(Lobster.list,density=T)

#####
