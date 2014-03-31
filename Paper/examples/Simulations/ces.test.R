ces.tests <- function(seps=seq(0.5,5,0.5),n.fats=12,n.preys=5,nsamples=20,n.preds=6){

    require(FASTIN)

     clo <- function(x){
        if (is.null(dim(x))){xc <- x/sum(x)} else {xc <- t(apply(x,1,function(y){y/sum(y)}))}
        return(xc)
    }

    alr <- function(x){
        x<-clo(x)
        if (is.null(dim(x))){xc <- log(x[1:(length(x)-1)]/x[length(x)])} else { t(apply(x,1,function(y){log(y[1:(length(y)-1)]/y[length(y)])}))}
    }

    adist <- function(mat,R){

        dims <- dim(mat)
        dists <- matrix(,dims[1],dims[1])
        for (i in 1:(dims[1]-1)){
            for (j in (i+1):dims[1]){
                dists[j,i] <- mahalanobis(alr(mat[i,]),alr(mat[j,]),R[,,j])
            }}
        dista <- as.dist(dists)
              return(dista)
    }

    kappas <- function(mat,vars){
        ks <- vector(,length(vars))
        for (i in 1:length(vars)){
            ks[i] <- kappa(mat[,vars[1:i]],exact=T)
        }
        return(ks)
    }
    
    # loop over ks

    cc_test <- array(,c(length(seps),4))
   

    a=0
    for (sep in seps){
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
    
    # set variance of fat content to 0.4 - this is loosely based on literature
    fc_tau <- rep(1/0.4,n.preys)
    
# draw mean conversion coeffs - cvar is the variance of mean conversion coeffs around 1; 0.2 seems reasonable  (c.f., Rosen & Tollit histogram)

    mean_cs=matrix(rlnorm(n.fats*n.preys,0,0.2),n.preys,n.fats) # rtnorm is a truncated normal, with lower truncation set to 0
    mean_c=mean_cs # just so I can modify mean_cs later without loosing the actual mean
    var_coeffs = matrix((0.03)^2,n.preys,n.fats)  # set the variance to some constant - pretty close, at least for FA whose CC isn't too far from 1 (Iverson 2004 graphs)
    tau_coeffs = 1/var_coeffs

# draw individual coeffs for each predator (e.g., may depend on individual diet status - fasting vs low fat diet versus hight fat diet) from a distribution about mean coeffs
#cs <- t(matrix(rlnorm(n.preys*n.fats,log(mean_c),sqrt(var_coeffs)),n.fats,n.preys))
#hist(cs)
# these are the proportion that each predator eats of each prey
# population 'mean', Dirichlet alpha notation, hard to interpret but computationally efficient.
# v(1,2) determines the variance of the alpha coeffs - the higher this variance, the lower the among predator variation in diet proportions.
    v1 = 0.5; v2 = 0.3
    Q = rlnorm(n.preys,rlnorm(n.preys,0,v1),v2)
# individual proportions

    
    props=rep(1/n.preys,n.preys)
    maxent <- -sum(props*log(props))
    
    props=MCMCpack::rdirichlet(n.preds,Q)

    this.ent=-sum(colMeans(props)*log(colMeans(props)))

# the mixture in matrix algebra (just a matrix product)
# prey mean is the geometric mean of each column
    mprey <-  clo(t(apply(preys,3,function(x){apply(x,2,function(y){exp(mean(log(y)))})})))
    preda <-  clo(props%*%(fc_mean*mprey*mean_c))
 
                                        # normalize and transform
        preds <- unclass(data.matrix(alr(preda[,])))

                                        # transform prey stats
        preym <-  alr(mprey[,])
        
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

        datas.FA <- list(n.fats = n.fats,m.fats=m.fats,fc_mean=fc_mean,fc_tau =fc_tau,R=R,Rnot = Rnot,preym=preym,preds = preds,ni=rep(nsamples,n.preys),mean_c=mean_c,tau_c = tau_coeffs)

        test.data <- list(datas.FA=datas.FA,n.preys=n.preys,n.preds=n.preds,even=0.05)
    
# MCMC defaults to population proportions only from FA data,
        outs <- run_MCMC(nIter=10000,nBurnin=1000,nChains=1,nThin=10,datas = test.data,plott=F)   

        cc_test[a,1] <- this.ent/maxent
        cc_test[a,2] <- kappa(preya,exact=T)
        cc_test[a,3] <- mean(adist(mprey,R))
        cc_test[a,4] <- robCompositions::aDist(colMeans(outs[[1]]),colMeans(props))

  }
    return(cc_test)
}

