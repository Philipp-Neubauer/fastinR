prey.table <- read.csv('prey_SI.csv',header=T)
pred.table <- read.csv('predator_SI.csv',header=T,stringsAsFactors=F)

prey.type <- as.numeric(prey.table[,1] == 'Grass shrimp')+1
n.preys <- 2 
n.preys.samps <- length(prey.type)

prey <- prey.table[,2:3]

# feed type index
feed.type <- pred.table[,1]
feed.type[pred.table[,1] == 'F' | pred.table[,1] == 'SF'] <- 1
feed.type[pred.table[,1] == 'C' | pred.table[,1] == 'SC'] <- 2
# use only fish and crustations diets, rest will be assessed
idx <- which(feed.type %in% c(1,2))
feed.type <- as.numeric(feed.type[idx])

pred <- pred.table[idx,2:3]
n.preds.samps <- length(idx)

# Set up priors ------

# prior for \Delta N follows from Hussey et al 2014
# beta.not - 95% interval contains 2SD, if symmetric then
bnot.prior.N <- 5.92
sd.bnot.N <- 5.92-4.55
bnot.tau.prior.N <- 1/(sd.bnot.N^2)

beta.prior.N <- -0.27
sd.beta.N <- -0.27--0.41
beta.tau.prior.N <- 1/(sd.beta.N^2)

# prior for \Delta C is adapted from Caut et al. 2009
# 95% interval contains 2SD, if symmetric then
bnot.prior.C <- -2.85
sd.bnot.C <- 0.5
bnot.tau.prior.C <- 1/(sd.bnot.C^2)

beta.prior.C <- -0.21
sd.beta.C <- 0.05
beta.tau.prior.C <- 1/(sd.beta.C^2)

# final priors
bnot.prior <- c(bnot.prior.C,bnot.prior.N)
bnot.tau.prior <- c(bnot.tau.prior.C,bnot.tau.prior.N)
beta.prior <- c(beta.prior.C,beta.prior.N)
beta.tau.prior <- c(beta.tau.prior.C,beta.tau.prior.N)


require(dplyr)
# predator mean priors
prior.mu <- data.frame(feed = prey.type,prey) %.% group_by(feed)  %.% summarise(mu.C = mean(X_13C),mu.N = mean(X_15N)) %.% arrange(feed)
prior.mu <- as.matrix(prior.mu[,2:3])


# fix sigmas?
sigma = data.frame(feed = prey.type,prey) %.% group_by(feed)  %.% summarise(mu.C = 1/var(X_13C),mu.N = 1/var(X_15N)) %.% arrange(feed)
sigma <- as.matrix(sigma[,2:3])

sigma.pred = data.frame(feed = feed.type,pred) %.% group_by(feed)  %.% summarise(mu.C = 1/var(X_13C),mu.N = 1/var(X_15N)) %.% arrange(feed)
sigma.pred <- as.matrix(sigma.pred[,2:3])

input <- list(prey=prey,
              pred=pred,
              n.preys.samps=n.preys.samps,
              n.preys=n.preys,
              n.preds.samps=n.preds.samps,
              prey.type=prey.type,
              feed.type=feed.type,
              prior.mu=prior.mu,
              #sigma=sigma,
              #sigma.pred=sigma.pred,
              bnot.prior=bnot.prior,
              beta.prior=beta.prior,
              bnot.tau.prior=bnot.tau.prior,
              beta.tau.prior=beta.tau.prior)

require(rjags)

DM <- jags.model('../../Discrim.model.SI.R',n.chains=3,inits = list(mu=prior.mu,beta.not=bnot.prior,beta.reg=beta.prior),data=input)

update(DM,100000)

samps <- coda.samples(DM,c('pred.discr','beta.reg','beta.not','mu','sigma'),n.iter=1e5,thin=100)

x11()
plot(samps,ask=T)
summary(samps)

# get estiamted discrimintation from all chains:
ix <- grep('pred.discr',colnames(samps[[1]]))

r.samps <- do.call('rbind',samps)[,ix]
dim(r.samps)

#into format for analysis
discr.means <- matrix(apply(r.samps,2,mean),2,2)
discr.var <- matrix(apply(r.samps,2,var),2,2)

# write to file
colnames(discr.means) <- colnames(prey.table)[2:3]
rownames(discr.means) <- c('fish','shrimp')
write.csv(discr.means,file='discr.means')

# write to file
colnames(discr.var) <- colnames(prey.table)[2:3]
rownames(discr.var) <- c('fish','shrimp')
write.csv(discr.var,file='discr.var')

####
### do the same for fatty acids -----
####

prey.table.FA <- t(read.csv('prey_FA.csv',header=T,stringsAsFactors=F,row.names=1))
pred.table.FA <- t(read.csv('predator_FA.csv',header=T,stringsAsFactors=F,row.names=1))

prey.type <- rep(1,nrow(prey.table.FA))
prey.type[grep('Grass.Shrimp',rownames(prey.table.FA))] <- 2

n.preys <- 2 
n.preys.samps <- length(prey.type)

for (i in 1:ncol(prey.table.FA))
  prey.table.FA[prey.table.FA[,i]==0,i] = min(prey.table[prey.table.FA[,i]>0,i])
prey=alr(prey.table.FA)

for (i in 1:ncol(pred.table.FA))
  pred.table.FA[pred.table.FA[,i]==0,i] = min(pred.table.FA[pred.table.FA[,i]>0,i])
pred=alr(pred.table.FA)

# feed type index
feed.type <- pred.table.FA[,1]
feed.type[grep('F',rownames(pred.table.FA))] <- 1
feed.type[grep('C',rownames(pred.table.FA))] <- 2
# use only fish and crustations diets, rest will be assessed
idx <- which(feed.type %in% c(1,2))
feed.type <- as.integer(feed.type[idx])

pred <- pred[idx,]
n.preds.samps <- length(idx)

# Set up FA priors ------
# not sure here, conflicting evidence for many FA's, better 

n.fats <- ncol(prey.table)
m.fats <- n.fats-1
S <- diag(0.1,m.fats)
R <- diag(0.1,n.fats)

require(dplyr)
# predator mean priors
pt <- aggregate(prey.table.FA,list(prey.type),gmean)
prior.mu <- data.matrix(alr(pt[,2:length(pt)]))
p=rep(1/n.fats,n.fats)
zeros = rep(0,n.fats)

# input <- list(prey=(prey),
#               n.fats=n.fats,
#               m.fats=m.fats,
#               pred=data.matrix(pred),
#               n.preys.samps=n.preys.samps,
#               n.preys=n.preys,
#               n.preds.samps=n.preds.samps,
#               prey.type=prey.type,
#               feed.type=feed.type,
#               prior.mu=prior.mu,
#               p=p
#               )

require(rjags)

DM <- jags.model('../../Discrim.model.FA.R',n.chains=2)

update(DM,1000)

samps <- coda.samples(DM,c('beta.reg'),n.iter=1e4,thin=10)

x11()
plot(samps,ask=T)

summary(samps)

# get estiamted discrimintation from all chains:

r.samps <- do.call('rbind',samps)
dim(r.samps)

fish.cc.samples <- r.samps[,seq(1,2*n.fats,2)]
shrimp.cc.samples <- r.samps[,seq(2,2*n.fats,2)]

fish.cc <- colMeans(fish.cc.samples)
shrimp.cc <- colMeans(shrimp.cc.samples)
#combine
ccs <- rbind(fish.cc,shrimp.cc)
colnames(ccs) <- colnames(prey.table)
write.csv(ccs,file='cc_FA.csv')

crosscorr.plot(samps)
# using independent ccs seems warranted here
fish.cc.var <- apply(fish.cc.samples,2,var)
shrimp.cc.var <- apply(shrimp.cc.samples,2,var)
# combine
ccs.var <- rbind(fish.cc.var,shrimp.cc.var)
colnames(ccs.var) <- colnames(prey.table)
write.csv(ccs.var,file='cc_FA_var.csv')


cor(fish.cc.samples)
