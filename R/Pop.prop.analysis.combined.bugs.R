Pop.prop.analysis.combined.bugs <- function(){
  for (j in 1:n.preys){
    
    prop[j] <- ps[j]/sum(ps[1:n.preys] )
    ps[j] ~ dgamma(0.1,1)%_%I(0.01,)
    
    # draw fat content from prior
    fc[j] ~ dnorm(fc_mean[j],fc_tau[j])%_%I(0.001,)    
    
    # take out of joint posterior for mixing model = make analysis conditional on fc
    #fcc[j] <- cut(fc[j])
    
    # model for sources with Jeffrey's (uninformative, improper) prior
    prey.precs[j,1:m.fats,1:m.fats] ~ dwish(R[,,j],ni[j])
    prey.means[j,1:m.fats] ~  dmnorm(preym[j,],Rs[j,,])    
    
    # same for SI
    prey.precs_SI[j,1:isos,1:isos] ~ dwish(R_SI[,,j],ni.SI[j])
    prey.means_SI[j,1:isos] ~  dmnorm(preym.SI[j,],Rs_SI[j,,])
    
    for (i in 1:isos)
    { cs[j,i] ~ dnorm(mean_cs[i],tau_cs[i])
      # take out of joint posterior for mixing model = make analysis conditional on prey mean distribution
      prey.mc_SI[j,i] <-cut(prey.means_SI[j,i])
      prey_SI[j,i] <- cs[j,i]+prey.mc_SI[j,i]
      
      for (k in 1:isos)
      {
        Rs_SI[j,i,k] <- prey.precs_SI[j,i,k]*ni.SI[j]
      }
    }
    
    for (f in 1:m.fats){
      # take out of joint posterior for mixing model = make analysis conditional on prey mean distribution
      prey.mc[j,f] <-cut(prey.means[j,f]) 
            
      for (g in 1:m.fats){ 	
        Rs[j,f,g] <- prey.precs[j,f,g]*ni[j]
      }
      
      # backtransform to apply fractionation and fat content
      preydd[j,f] <- exp(prey.mc[j,f])
      preyd[j,f] <- preydd[j,f]/(sum(preydd[j,1:m.fats])+1)
      
    }
    
    preyd[j,n.fats] <- 1/(sum(preydd[j,1:m.fats])+1)
    
    for (f in 1:n.fats){      
      # draw fractionation coeffs for prey j and fatty acid f
      fracs[j,f] ~ dnorm(mean_c[j,f],tau_coeffs[j,f])%_%I(0,)
      frac[j,f] <- cut(fracs[j,f])
      preyf[j,f] <- (preyd[j,f]*frac[j,f]*fc[j])
    }
    
  }
  
  for(i in 1:n.preds) {     
       
    # predator likelihood
    preds[i,1:m.fats] ~ dmnorm(mu[],predprec[,])
    
    preds.SI[i,1:isos] ~ dmnorm(mu_SI[],predprec_SI[,])
      
    }
 
  # mixing for predator SI signature for likelihood
  for (k in 1:isos){ 
    mu_SI[k] <- inprod(prop[1:n.preys],prey_SI[1:n.preys,k])
    
  }
  
  # mixing and alr transformation of predator fat signature for likelihood
  for (f in 1:n.fats){
    
    predm[f] <- inprod(prop[1:n.preys],preyf[1:n.preys,f])
    predf[f] <- predm[f]/sum(predm[1:m.fats])
  }
  for (f in 1:m.fats){
    mu[f] <- log(predf[f]/predf[n.fats])
    
  }  
  
  # priors for proportions and predator covariance
  
  predprec[1:m.fats,1:m.fats] ~ dwish(Rnot[,],m.fats)
  predprec_SI[1:isos,1:isos] ~ dwish(Rnot_SI[,],isos)
  
}
