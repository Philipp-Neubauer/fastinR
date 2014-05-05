model{
  for (i in 1:2){
    for (j in 1:n.preys.samps){
    
      # prey samples
      prey[j,i] ~ dnorm(mu[prey.type[j],i],tau[prey.type[j],i])
      
    }
    for (k in 1:n.preys){
      #prey mean and precision
      mu[k,i] ~ dnorm(prior.mu[k,i],0.001)
      tau[k,i] ~ dgamma(0.1,0.1)
      
      # predator variance and discrimination variance
      tau.pred[k,i] ~ dgamma(0.1,0.1)
      tau.pred.mu[k,i] ~ dgamma(0.1,0.1)
      
      pred.mu[k,i] <- mu[k,i] + beta.not[i] + beta.reg[i]*mu[k,i]
      pred.discr[k,i] <- beta.not[i] + beta.reg[i]*mu[k,i]
    }
    
    for (n in 1:n.preds.samps)
    {
      pred[n,i] ~ dnorm(pred.mu[feed.type[n],i],tau.pred[feed.type[n],i])
    }
    # regression priors
    beta.not[i] ~ dnorm(bnot.prior[i],bnot.tau.prior[i])
    beta.reg[i] ~ dnorm(beta.prior[i],beta.tau.prior[i])
  }
  
}