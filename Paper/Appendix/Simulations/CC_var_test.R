test.cc.dep <- function(k=0.1,n.fats=8,n.preys=3,nsamples=20,sep=3,n.preds=6){

    require(fastinR)

    # loop over ks

    cc_test <- array(,c(length(k),3))
    a=0
    for (cvar in k){
        a=a+1
                                        # simulations for CC var dependence --------------

    As = matrix(rlnorm(n.preys*n.fats,rnorm(n.preys*n.fats,sep)),n.preys,n.fats,byrow=T)
    
    preys <- array(,c(nsamples,n.fats,n.preys))
    for (i in 1:n.preys){
      preys[,,i] <- MCMCpack::rdirichlet(nsamples,As[i,])
    }
 
    preys.ix <- rep(1:n.preys,each=nsamples)
    preya<- apply(preys,2,I)
   
    #randomly draw fat content
    fc_mean <- (sample(40,n.preys))/10+1
    
    fc.var = log(0.4 + fc_mean^2) - 2*log(fc_mean)
    fc_mean = log(fc_mean)-fc.var/2
        
    # set variance of fat content to 0.4 - this is loosely based on literature
    fc_tau <- 1/fc.var
        
# draw mean conversion coeffs - cvar is the variance of mean conversion coeffs around 1; 0.2 seems reasonable  (c.f., Rosen & Tollit histogram)

mean_cs=matrix(rlnorm(n.fats*n.preys,0,cvar),n.preys,n.fats) 
mean_c=mean_cs # just so I can modify mean_cs later without loosing the actual mean
var_c= matrix((0.03)^2,n.preys,n.fats)  # set the variance to some constant - pretty close, at least for FA whose CC isn't too far from 1 (Iverson 2004 graphs)
     
 # covert to gamma parameters
  rate <- mean_c/var_c
  shape <- mean_c^2/var_c
  
  mean_c <- shape
  var_c <- rate
        

v1 = 0.5; v2 = 0.7
Q = rlnorm(n.preys,rlnorm(n.preys,0,v1),v2)
# individual proportions
props=MCMCpack::rdirichlet(n.preds,Q)


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

 R[,,i]=(cov(alr(preys[,,i]))*(nsamples-1))
 
}

# set uninformative prior SS matrix for wishart prior alr transformed predator data
Rnot =diag(0.01,m.fats)

#############################################################
########### make data object and run analysis ###############
#############################################################

# Fatty Acid data (stable isotopes are intigrated in the same way

    datas.FA <- list(n.fats = n.fats,m.fats=m.fats,fc_mean=fc_mean,fc_tau =fc_tau,R=R,Rnot = NULL,preym=preym,preds = preds,preds.FA=preda,preys=apply(preys,2,I),ni=rep(nsamples,n.preys),mean_c=mean_c,tau_c = var_c)

        test.data <- list(datas.FA=datas.FA,n.preys=n.preys,n.preds=n.preds,even=0.05,prey.ix=preys.ix)
  
# MCMC defaults to population proportions only from FA data,
    outs <- run_MCMC(nIter=10000,nBurnin=1000,nChains=1,nThin=10,datas = test.data,plott=F)   

    cc_test[a,1] <- cvar
    cc_test[a,2] <- robCompositions::aDist(colMeans(outs[[1]]),colMeans(props))

########## set mean_c to 1
mean_cs=matrix(1,n.preys,n.fats) 
mean_c=mean_cs # 
var_c= matrix(1,n.preys,n.fats)  #
     
 # covert to gamma parameters
  rate <- mean_c/var_c
  shape <- mean_c^2/var_c
  
  mean_c <- shape
  var_c <- rate
   
    
    datas.FA <- list(n.fats = n.fats,m.fats=m.fats,fc_mean=fc_mean,fc_tau =fc_tau,R=R,Rnot = NULL,preym=preym,preds = preds,preds.FA=preda,preys=apply(preys,2,I),ni=rep(nsamples,n.preys),mean_c=mean_c,tau_c = var_c)

        test.data <- list(datas.FA=datas.FA,n.preys=n.preys,n.preds=n.preds,even=0.05,prey.ix=preys.ix)
  
    
# MCMC defaults to population proportions only from FA data,
    outs <- run_MCMC(nIter=10000,nBurnin=1000,nChains=1,nThin=10,datas = test.data,plott=F)   

    cc_test[a,3] <- robCompositions::aDist(colMeans(outs[[1]]),colMeans(props))
  }
    return(cc_test)
}

