Analysis.with.Cov.combined.bugs <- function(){
  for (j in 1:n.preys){
  
 # draw fat content from prior
    fc[j] ~ dnorm(fc_mean[j],fc_tau[j])%_%I(0.001,)    
    
    # take out of joint posterior for mixing model = make analysis conditional on fc
    #fcc[j] <- cut(fc[j])
    
    # model for sources with Jeffrey's (uninformative, improper) prior
    prey.precs[j,1:m.fats,1:m.fats] ~ dwish(R[,,j],ni[j])
    prey.means[j,1:m.fats] ~  dmnorm(preym[j,],Rs[j,,])    

  # model for sources with Jeffrey's (uninformative, improper) prior
   
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
    
       # draw clr(proportions)
   pnorm[i,1:n.preys] ~ dmnorm(mus[i,],pprec[,])
    for(j in 1:n.preys) {
      
      # backtransform clr proportions by clr^-1
      p_unn[i,j] <- exp(pnorm[i,j])
            
      #closure
      prop[i,j] <- (p_unn[i,j])/sum(p_unn[i,1:n.preys])
      
      mus[i,j] <-inprod(beta[j,1:n.covs],Covs[i,1:n.covs])
    }
    
   # predator likelihood
    
    preds.SI[i,1:isos] ~ dmnorm(mu_SI[i,],predprec_SI[,])
      
    
 
  # mixing for predator SI signature for likelihood
  for (k in 1:isos){ 
    mu_SI[i,k] <- inprod(prop[i,1:n.preys],prey_SI[1:n.preys,k])
    
  }
 # predator likelihood
    preds[i,1:m.fats] ~ dmnorm(mu[i,],predprec[,])
    
    # mixing and alr transformation of predator signature for likelihood
    for (f in 1:n.fats){ 
      predm[i,f] <- inprod(prop[i,1:n.preys],preyf[1:n.preys,f])
      predf[i,f] <- predm[i,f]/sum(predm[i,1:m.fats])
    }
    for (f in 1:m.fats){
      mu[i,f] <- log(predf[i,f]/predf[i,n.fats])
     }        

  }   
        
# regression coeffs
for (k in 1:n.covs){
  beta[1:n.preys,k] ~ dmnorm(zeros[],S[,])
}

  # priors for proportions and predator covariance
predprec[1:m.fats,1:m.fats] ~ dwish(Rnot[,],m.fats)
   predprec_SI[1:isos,1:isos] ~ dwish(Rnot_SI[,],isos)
  
  for(j in 1:n.preys) {  
    pm_unn[j] <- exp(beta[j,1])       
    #closure
    pop.prop[j] <- (pm_unn[j])/sum(pm_unn[1:n.preys])
  }
    
  pprec[1:n.preys,1:n.preys] ~ dwish(SS[,],n.preys)
}
