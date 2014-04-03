model{
  
    for (j in 1:n.preys.samps){
      prey[j,1:m.fats] ~ dmnorm(mu[prey.type[j],1:m.fats],sigma[prey.type[j],1:m.fats,1:m.fats])
    }
    
    for (k in 1:n.preys){
      mu[k,1:m.fats] ~ dmnorm(prior.mu[k,1:m.fats],S[1:m.fats,1:m.fats])
      sigma[k,1:m.fats,1:m.fats] ~ dwish(S[1:m.fats,1:m.fats],m.fats)
      
      sigma.pred[k,1:m.fats,1:m.fats] ~ dwish(S[1:m.fats,1:m.fats],m.fats)
      sigma.pred.mu[k,1:m.fats,1:m.fats] ~ dwish(S[1:m.fats,1:m.fats],m.fats)
      
      # is also the number of feeds
      
      pred.mu[k,1:m.fats] ~ dmnorm(mu.p[k,1:m.fats],sigma.pred.mu[k,1:m.fats,1:m.fats])
      
      #closure & log-ratio transform
      for (f in 1:m.fats)
      {
        l.mu[k,f]  <- mu.p.org[k,f]/mu.p.org[k,n.fats]
        mu.p[k,f] <- log(l.mu[k,f] )
        mu.o[k,f] <- exp(mu[k,f])
        mu.p.org[k,f] <- beta.reg[k,f]*mu.org[k,f]
        # uniform dirichlet
        ps[k,f] ~ dgamma(1/n.fats,1)I(0.01,)
      }
      ps[k,n.fats] ~ dgamma(1/n.fats,1)I(0.01,)
      
      mu.p.org[k,n.fats] <- beta.reg[k,n.fats]*mu.org[k,n.fats]
      
      #closure
      mu.o[k,n.fats] <- 1
      mu.org[k,1:n.fats] <- mu.o[k,1:n.fats]/sum(mu.o[k,1:n.fats])
      # fractionation....
      beta.reg[k,1:n.fats] <- ps[k,1:n.fats]/sum(ps[k,1:n.fats] )
      
    }
    
    for (n in 1:n.preds.samps)
    {
      pred[n,1:m.fats] ~ dmnorm(pred.mu[feed.type[n],1:m.fats],sigma.pred[feed.type[n],1:m.fats,1:m.fats])
    }
    
  
}