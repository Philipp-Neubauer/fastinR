# simulate fatty acid proportions in 3 prey species ----
require(FASTIN)

# try with Squid data -----------

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

#setwd("../Squid/")

SQ <- read.csv("Squid data.csv",header=T,sep=',')
PR <- read.csv("Squid prey data.csv",header=T)

# SQuid means and SD from std error - divided by 100 to get proportion rather than percent
SQ.means <- t(apply(SQ[2:nrow(SQ),seq(1,ncol(SQ),2)],2,clo) )
SQ.sd <- (apply(data.matrix(SQ[2:nrow(SQ),seq(2,ncol(SQ),2)]),1,function(x){x*sqrt(data.matrix(SQ[1,seq(2,ncol(SQ),2)]))}))
SQ.sd <- SQ.sd/100

n.preds = sum(as.numeric(SQ[1,c(1,3)]))


###########################################################
###### Data grooming: combining sources, ##################
###### and a few (minor) assumptions...  ##################
###########################################################

#get fat contents

fc.means <- as.numeric(PR[2,seq(1,ncol(PR),2)]/100)
fc.sd = as.numeric(PR[2,seq(2,ncol(PR),2)]/100)

# combine sources
fc.mean.new = c(exp(mean(log(fc.means[1:4]))),exp(mean(log(fc.means[5:8]))),exp(mean(log(fc.means[c(9,11:12,14)]))),fc.means[c(10,13)])

fc.sd.new = c(mean(fc.sd[1:4]),mean(fc.sd[5:8]),mean(fc.sd[c(9,11:12,14)]),t(fc.sd[c(10,13)]))

# Prey means and SD
PR.means <- apply(PR[3:nrow(PR),seq(1,ncol(PR),2)],2,clo)

# replace zeros in means with min/10 (pretending we had more than 1 sample...) 
for (i in 1:nrow(PR.means))
  PR.means[i,PR.means[i,]==0] <- min(PR.means[i,PR.means[i,]!=0])/10

# pool data for crustatcians, myctophid fish
# pool means using the geometric mean
PR.means.new <- cbind(apply(PR.means[,1:4],1,function(x){exp(weighted.mean(log(x),w=clo(fc.means[1:4])))}),apply(PR.means[,5:8],1,function(x){exp(weighted.mean(log(x),w=clo(fc.means[5:8])))}),apply(PR.means[,c(9,11:12,14)],1,function(x){exp(weighted.mean(log(x),w=clo(fc.means[c(9,11:12,14)])))}),PR.means[,c(10,13)])

PR.means.new <- t(apply(PR.means.new,2,clo))

#SDs

PR.sds <- t(apply(data.matrix(PR[3:nrow(PR),seq(2,ncol(PR),2)]),1,function(x){x*sqrt(data.matrix(PR[1,seq(2,ncol(PR),2)]))}))
PR.sd <- PR.sds/100

# use simple sd instead of geometric sd since we do not have the original data
PR.sd.new <- cbind(apply(PR.sd[,1:4],1,function(x){sqrt(sum(clo(fc.means[1:4])^2*x^2))}),apply(PR.sd[,5:8],1,function(x){sqrt(sum(clo(fc.means[5:8])^2*x^2))}),apply(PR.sd[,c(9,11:12,14)],1,function(x){sqrt(sum(clo(fc.means[c(9,11:12,14)])^2*x^2))}),PR.sd[,c(10,13)])

# replace zeros in sds with min sd (pretending we had more than 1 sample...) 
for (i in 1:nrow(PR.sd.new))
  PR.sd.new[i,PR.sd.new[i,]==0] <- min(PR.sd.new[i,PR.sd.new[i,]!=0]) 

PR.sd.new <- t(PR.sd.new)

nsample = as.numeric(PR[1,seq(1,ncol(PR-1),2)])
nsamples = c(sum(nsample[1:4]),sum(nsample[5:8]),sum(nsample[c(9,11:13)]),nsample[10],nsample[14])
# remove n and fat content
#PR<- PR[-c(1,2),]

n.fats=ncol(PR.means.new)

#################################################
#### Plot data in multivariate (CAP) space ######
#################################################

dista <- adist(PR.means.new)

PR.RDA <- vegan::capscale(dista~as.factor(1:5),comm=PR.means.new)
#plot(PR.RDA,t='n',xlim=c(-0.5,0.5),ylim=c(-1,1))
plot(t(apply((PR.means.new)%*%PR.RDA$CCA$v[,1:2],1,function(x){x})),pch=1:6)
points((SQ.means[1,])%*%PR.RDA$CCA$v[,1:2],pch=16)
points((SQ.means[2,])%*%PR.RDA$CCA$v[,1:2],pch=17)
legend('bottomleft',c('Crustaceans','Myctophid Fish','Other Fish','Salilota australis','Loligo gahi','Small Squid','Large Squid'),pch=c(1:5,16,17))

par(mfcol=c(2,1))
sv = sort(clo(rowSums(t(t(PR.RDA$CCA$v)*PR.RDA$CCA$eig)^2)),decreasing =T,index.return=T)
plot(cumsum(sort(clo(rowSums(t(t(PR.RDA$CCA$v)*PR.RDA$CCA$eig)^2)),decreasing =T)),axes=F,xlab='',ylab='Cumulative source separation',ylim=c(0,1))
axis(2)
axis(1,at=1:n.fats,labels=F)
text(1:n.fats, par("usr")[3] - 0.2, srt = 45, adj = 1,
     labels = colnames(PR.means.new)[sv$ix], xpd = TRUE)

ks <- kappas(PR.means.new,sv$ix)
plot(ks,axes=F,ylab='Prey matrix condition number',xlab='',ylim=c(0,max(ks)))
axis(2)
axis(1,at=1:n.fats,labels=F)
text(1:n.fats, par("usr")[3] - max(ks)/5 , srt = 45, adj = 1,
     labels = colnames(PR.means.new)[sv$ix], xpd = TRUE)
   

 nv <- select.list(title='please choose the fatty acids to use (at least 3)',choices = colnames(PR.means.new)[sv$ix],graphics=T,multiple=T)
      
      #nv <- readline(prompt = "please enter number of variables for analysis \n")
six <- match(nv,colnames(PR.means.new))

n.fats=length(six)
n.preys=nrow(PR.means.new)

#################################
###### simulation study #########
#################################

# number of predators to simulate
n.preds=4

# since we do not have individual prey samples from the study, we simulate them, to get an idea of the correlation structure in the composition introduced by the sum constraint itself. More complex correlations tructure is unfortuantelly lost...

#number of simulated prey samples per species/group
ns <- 30
PR.sim <- array(,c(ns,n.fats,n.preys))

for (i in 1:n.preys) PR.sim[,,i] <- mvrnorm(ns,PR.means.new[i,six],diag(PR.sd.new[i,six]^2))

for (i in (1:n.preys))
PR.sim[,,i] <- apply(PR.sim[,,i],2,function(x){x[x<=0] = min(x[x>0])/10;return(x)})


v1 = 0.1; v2 = 0.3
Q = rlnorm(n.preys,rlnorm(n.preys,0,v1),v2)
# individual proportions
props=MCMCpack::rdirichlet(n.preds,Q)
# show simulated proportions
props

mprey <-  data.matrix(clo(PR.means.new))
preda <-  clo((props) %*%(as.numeric(fc.mean.new)*mprey))
# normalize and transform
preds <- unclass(data.matrix(alr(preda[,six])))

# transform prey stats
preym <-  unclass(alr(mprey[,six]))

m.fats = (n.fats-1)


#####################################################
############ Priors #################################
#####################################################

# get sums of squares for each prey matrix

R <- array(,c(m.fats,m.fats,n.preys))
ni<-rep(NA,n.preys)
for (i in 1:n.preys){
    ni[i] <- max(n.fats+1,n.preys-1)
    R[,,i]=cov(alr(PR.sim[,,i]))*ni[i]
}

# set uninformative prior SS matrix for wishart prior alr transformed predator data
Rnot =diag(0.01,m.fats)    
# # set SS matrix for wishart prior on proportions
even=0.05
2*(1-pnorm(log(95)/2,0,sqrt(1/even))) # chance that 95*p1<p2

mean_c = matrix(1,n.preys,n.fats)
tau_coeffs = matrix(10000,n.preys,n.fats)

fc_mean = as.numeric(fc.mean.new)
fc_tau = 1/as.numeric(fc.sd.new^2)

#############################################################
########### make data object and run analysis ###############
#############################################################

# Fatty Acid data (stable isotopes are intigrated in the same way

datas.FA <- list(n.fats = n.fats,m.fats=m.fats,fc_mean=fc_mean,fc_tau =fc_tau,R=R,Rnot = Rnot,preym=preym,preds = preds,ni=ni,mean_c=mean_c,tau_c = tau_coeffs)

Squid.data <- list(datas.FA=datas.FA,n.preys=n.preys,n.preds=n.preds,even=even,prey.ix=c('Crustaceans','Myctophid Fish','Other Fish','Salilota australis','Loligo gahi'))

guiSet('datas',Squid.data)

# MCMC defaults to population proportions only from FA data,
outs <- run_MCMC(nIter=50000,nBurnin=10000,nChains=3,nThin=50,datas = Squid.data)

#compare props to outs
 colMeans(outs[[1]])[1:n.preys]-colMeans(props)

#plot outputs
plot(outs)

summary(outs)

preya.names <- unique(guiGetSafe('datas')$prey.ix)

outs <- matrix(unlist(x),(mcpar(x[[1]])[2]-(mcpar(x[[1]])[1]-mcpar(x[[1]])[3]))/mcpar(x[[1]])[3],length(preya.names))

source('ggdensity.R')
ggs_density(outs,colMeans(props))
ggsave('pop_prop_squid.pdf')

# now try individual proportions from FA data,
outs.ind <- run_MCMC(nIter=100000,nBurnin=10000,nChains=1,nThin=10,datas = Squid.data,Analysis.Type='Individual.proportions')

#compare props to outs
colMeans(outs.ind[[1]])[1:n.preys]-colMeans(props)

#plot outputs
plot(outs.ind)

pop.ix <- grep('pop',names(outs[[1]]))

ggs_density(outs[[1]][,pop.ix],colMeans(props))





#############################################
##### look at tradeoffs           ###########
#############################################

# these simualtions produce nsims independent simpulations and estiamtions to
# look at correlations in errors made when estimating diets

nsims=50
delta.p <- matrix(,nsims,n.preys)

for (k in 1:nsims){
    n.preds=1
    
    v1 = 0.5; v2 = 0.7
    Q = rlnorm(n.preys,rlnorm(n.preys,0,v1),v2)
                                        # individual proportions
    props= MCMCpack::rdirichlet(n.preds,Q)
                                       # re-calculate simulated predators       
    preda <-  clo((props) %*%(as.numeric(fc.mean.new)*mprey))
                                        # normalize and transform
    preds <- unclass(data.matrix(alr(preda[,six])))
    
                                        # new data object
    datas.FA <- list(n.fats = n.fats,m.fats=m.fats,fc_mean=fc_mean,fc_tau =fc_tau,R=R,Rnot = Rnot,preym=preym,preds = preds,ni=ni,mean_c=mean_c,tau_c = tau_coeffs)
    
    Squid.data <- list(datas.FA=datas.FA,n.preys=n.preys,n.preds=n.preds,even=even,prey.ix=c('Crustaceans','Myctophid Fish','Other Fish','Salilota australis','Loligo gahi'))
    
                                        # MCMC defaults to population proportions only from FA data,
    outs <- run_MCMC(nIter=10000,nBurnin=1000,nChains=1,nThin=10,datas = Squid.data,plott=F)
    
    delta.p[k,] = colMeans(outs[[1]])[1:n.preys]-colMeans(props)
}

delta.p <- as.data.frame(delta.p)
names(delta.p)<-c('Crustaceans','Myctophid Fish','Other Fish','Salilota australis','Loligo gahi')

pdf("boxplot_squid_errors.pdf")
boxplot(delta.p)
dev.off()

pdf("pairwise_squid_errors.pdf")
plot(delta.p)
dev.off()



###########################################
##### 'Real' analysis  ########
###########################################

########## small squid ############

n.preds =1

small <- SQ.means[1,]#clo(abs(mvrnorm(n.preds/2,SQ.means[1,],diag(SQ.sd[1,])^2)))

                                        # normalize and transform
preds <- unclass(data.matrix(alr(small[six])))

                                        # new data object
datas.FA <- list(n.fats = n.fats,m.fats=m.fats,fc_mean=fc_mean,fc_tau =fc_tau,R=R,Rnot = Rnot*100,preym=preym,preds = preds,ni=ni,mean_c=mean_c,tau_c = tau_coeffs)

Squid.data <- list(datas.FA=datas.FA,n.preys=n.preys,n.preds=n.preds,even=even,prey.ix=c('Crustaceans','Myctophid Fish','Other Fish','Salilota australis','Loligo gahi'))

guiSet('datas',Squid.data)

                                        # MCMC defaults to population proportions only from FA data,
outs.small <- run_MCMC(nIter=100000,nBurnin=5000,nChains=3,nThin=100,datas = Squid.data)

plot(outs.small)

summary(outs.small)


###### Large squid ###########

large <- SQ.means[2,]#clo(abs(mvrnorm(n.preds/2,SQ.means[2,],diag(SQ.sd[2,])^2)))

preds <- unclass(data.matrix(alr(large[six])))

                                        # new data object
datas.FA <- list(n.fats = n.fats,m.fats=m.fats,fc_mean=fc_mean,fc_tau =fc_tau,R=R,Rnot = Rnot*100,preym=preym,preds = preds,ni=ni,mean_c=mean_c,tau_c = tau_coeffs)

Squid.data <- list(datas.FA=datas.FA,n.preys=n.preys,n.preds=n.preds,even=even,prey.ix=c('Crustaceans','Myctophid Fish','Other Fish','Salilota australis','Loligo gahi'))

guiSet('datas',Squid.data)

                                        # MCMC defaults to population proportions only from FA data,
outs.large <- run_MCMC(nIter=100000,nBurnin=5000,nChains=3,nThin=100,datas = Squid.data)

plot(outs.large)

summary(outs.large)


Squid.list <- list('Small Squid'=outs.small,'Large Squid' = outs.large)
multiplot(Squid.list,density=F)
multiplot(Squid.list,density=T)
    
