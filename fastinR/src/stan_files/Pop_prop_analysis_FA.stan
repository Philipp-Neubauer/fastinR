functions {
  
   vector clo(vector pp)
    return pp/sum(pp);  
    
  vector alr(vector pp, int dim)
    return log(pp[1:(dim-1)]/pp[dim]);  
  
  vector inv_alr(vector pp, int dim){
    vector[dim] ppo;
    
    ppo[1:(dim-1)] = exp(pp);
    ppo[dim] = 1;
    return clo(ppo);  
  }
}

data {
  int<lower = 1> n_prey_samps;
  int<lower = 1> n_preys;
  int<lower = 1> n_preds;
  int<lower = 1> n_fats;
  int<lower = 1> m_fats;
  int prey_ix[n_prey_samps];
  vector[n_preys] fc_mean;
  vector[n_preys] fc_tau;
  vector[m_fats] preys[n_prey_samps];
  vector[m_fats] preym[n_preys];
  vector[n_fats] mean_c[n_preys];
  vector[n_fats] tau_coeffs[n_preys];
  vector[m_fats] preds[n_preds];
}
parameters{
  vector[m_fats] prey_means[n_preys];
  vector<lower=0>[n_fats] cs[n_preys];
  real<lower=0> fc[n_preys];
  vector[m_fats] cons_prey[n_preys];
  
  cholesky_factor_corr[m_fats] corr_prey[n_preys];
  cholesky_factor_corr[m_fats] corr_pred;
  cholesky_factor_corr[m_fats] corr_mean;
  vector<lower = 0>[m_fats] tau_prey[n_preys];
  vector<lower = 0>[m_fats] tau_mean;
  vector<lower = 0>[m_fats] tau_pred;
  vector<lower = 0>[n_preys] props;
}
transformed parameters{
  matrix[m_fats,m_fats] prey_precs[n_preys];
  matrix[m_fats,m_fats] pred_prec;
  simplex[n_preys] prop;
  matrix[n_fats,n_preys] prey;
  vector[n_fats] c_prey[n_preys];
  vector[m_fats] mu;
  
  prop = clo(props);
  
  for (j in 1:n_preys){
    
    prey_precs[j] = diag_pre_multiply(tau_prey[j], corr_prey[j]);
    c_prey[j,1:(n_fats-1)] = exp(cons_prey[j]);
    c_prey[j,n_fats] = 1;
    prey[,j] = fc[j] * cs[j] .* c_prey[j];
      // mixing for predator FA signature for likelihood
  }
  mu = alr(prey*prop, n_fats);
  
  pred_prec = diag_pre_multiply(tau_pred, corr_pred);
  
}
model{

//prey samples
  for(ps in 1:n_prey_samps)  preys[ps] ~ multi_normal_cholesky(prey_means[prey_ix[ps]], prey_precs[prey_ix[ps]]);

// predator likelihood
    
  for (p in 1:n_preds) preds[p] ~ multi_normal_cholesky(mu,pred_prec);

  props ~ gamma(1.5,1);
  fc ~ lognormal(fc_mean,sqrt(1.0./fc_tau));   

  for (j in 1:n_preys){
    
    tau_prey[j] ~ normal(0, 10);
    corr_prey[j] ~ lkj_corr_cholesky(3);
    
    prey_means[j] ~ multi_normal_cholesky(preym[j], diag_pre_multiply(tau_mean, corr_mean));
    cons_prey[j] ~ multi_normal_cholesky(prey_means[j],prey_precs[j]);

    cs[j] ~ gamma(mean_c[j],tau_coeffs[j]);
     
  }

      
  // priors for predator covariance
   tau_pred ~ normal(0, 10);
   corr_pred ~ lkj_corr_cholesky(3);
  
   tau_mean ~ normal(0, 10);
   corr_mean ~ lkj_corr_cholesky(3);
 
  
}
