model{
  for (i in 1:2){
    for (j in 1:n.preys.samps){
    
      prey[j,i] ~ dnorm(mu[prey.type[j],i],sigma[prey.type[j],i])
      
    }
    for (k in 1:n.preys){
      mu[k,i] ~ dnorm(prior.mu[k,i],0.001)
      sigma[k,i] ~ dgamma(0.1,0.1)
      
      
      sigma.pred[k,i] ~ dgamma(0.1,0.1)
      sigma.pred.mu[k,i] ~ dgamma(0.1,0.1)
      
      pred.mu[k,i] ~ dnorm(mu.p[k,i],sigma.pred.mu[k,i])
      mu.p[k,i] <- mu[k,i] + beta.not[i] + beta.reg[i]*mu[k,i]
      pred.discr[k,i] <- beta.not[i] + beta.reg[i]*mu[k,i]
    }
    
    for (n in 1:n.preds.samps)
    {
      pred[n,i] ~ dnorm(pred.mu[feed.type[n],i],sigma.pred[feed.type[n],i])
    }
    
    beta.not[i] ~ dnorm(bnot.prior[i],bnot.tau.prior[i])
    beta.reg[i] ~ dnorm(beta.prior[i],beta.tau.prior[i])
  }
  
}