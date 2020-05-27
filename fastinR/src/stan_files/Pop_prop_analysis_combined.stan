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
  int<lower = 1> n_prey_samps_SI;
  
  int<lower = 1> n_preys;
  int<lower = 1> n_preds;
  int<lower = 1> n_fats;
  int<lower = 1> m_fats;
  int<lower = 1> isos;
  
  int prey_ix[n_prey_samps];
  vector[n_preys] fc_mean;
  vector[n_preys] fc_tau;
  int<lower=0, upper=1> fc_data;
  vector[m_fats] preys[n_prey_samps];
  vector[m_fats] preym[n_preys];
  vector[n_fats] mean_c[n_preys];
  vector[n_fats] tau_coeffs[n_preys];
  vector[m_fats] preds[n_preds];
  
  int prey_ix_SI[n_prey_samps_SI];
  vector[isos] preys_SI[n_prey_samps_SI];
  vector[isos] preym_SI[n_preys];
  vector[isos] mean_cs[n_preys];
  vector[isos] sigma_cs[n_preys];
  vector[isos] preds_SI[n_preds];
}
transformed data{
  real fc_base[n_preys];
  
  for (p in 1:n_preys) fc_base[p] = 1.0;
}
parameters{
  vector[m_fats] prey_means[n_preys];
  vector<lower=0>[n_fats] cs[n_preys];
  real<lower=0, upper=10> fcs[fc_data ? n_preys : 0];
  vector[m_fats] cons_prey[n_preys];
  
  cholesky_factor_corr[m_fats] corr_prey[n_preys];
  cholesky_factor_corr[m_fats] corr_pred;
  cholesky_factor_corr[m_fats] corr_mean;
  vector<lower = 0>[m_fats] tau_prey[n_preys];
  vector<lower = 0>[m_fats] tau_mean;
  vector<lower = 0>[m_fats] tau_pred;
  
  vector[isos] prey_means_SI[n_preys];
  vector[isos] cc[n_preys];
  vector[isos] cons_prey_SI[n_preys];
  
  cholesky_factor_corr[isos] corr_prey_SI[n_preys];
  cholesky_factor_corr[isos] corr_pred_SI;
  cholesky_factor_corr[isos] corr_mean_SI;
  vector<lower = 0>[isos] tau_prey_SI[n_preys];
  vector<lower = 0>[isos] tau_mean_SI;
  vector<lower = 0>[isos] tau_pred_SI;
  
  
  vector<lower = 0>[n_preys] props;
}
transformed parameters{
  matrix[m_fats,m_fats] prey_precs[n_preys];
  matrix[m_fats,m_fats] pred_prec;
  matrix[n_fats,n_preys] prey;
  vector[n_fats] c_prey[n_preys];
  vector[m_fats] mu;
  
  matrix[isos,isos] prey_precs_SI[n_preys];
  matrix[isos,isos] pred_prec_SI;
  matrix[isos,n_preys] prey_SI;
  vector[isos] mu_SI;
  real fc[n_preys]; 
  simplex[n_preys] prop;
  
  
  prop = clo(props);
  
  if(fc_data){
    fc = fcs;
  } else {
    fc = fc_base;
  }
  
  for (j in 1:n_preys){
    
    prey_precs[j] = diag_pre_multiply(tau_prey[j], corr_prey[j]);
    c_prey[j,1:(n_fats-1)] = exp(cons_prey[j]);
    c_prey[j,n_fats] = 1;
    prey[,j] = fc[j] * cs[j] .* c_prey[j];
    // mixing for predator FA signature for likelihood
    
    prey_precs_SI[j] = diag_pre_multiply(tau_prey_SI[j], corr_prey_SI[j]);
    prey_SI[,j] = cc[j]+cons_prey_SI[j];
    // mixing for predator SI signature for likelihood
  }
  //for (i in 1:isos){
    mu_SI = prey_SI*prop;
    pred_prec_SI = diag_pre_multiply(tau_pred_SI, corr_pred_SI);
    
    mu = alr(prey*prop, n_fats);
    pred_prec = diag_pre_multiply(tau_pred, corr_pred);
    
}
model{
  
  //prey samples
  for(ps in 1:n_prey_samps)  preys[ps] ~ multi_normal_cholesky(prey_means[prey_ix[ps]], prey_precs[prey_ix[ps]]);
  for(ps in 1:n_prey_samps_SI)  preys_SI[ps] ~ multi_normal_cholesky(prey_means_SI[prey_ix_SI[ps]], prey_precs_SI[prey_ix_SI[ps]]);
  
  // predator likelihood
  for (p in 1:n_preds) {
    preds[p] ~ multi_normal_cholesky(mu,pred_prec);
    preds_SI[p] ~ multi_normal_cholesky(mu_SI,pred_prec_SI);
  }
  
  props ~ gamma(1.5,1);
  if(fc_data) fcs ~ lognormal(fc_mean,sqrt(1.0./fc_tau));  
  
  for (j in 1:n_preys){
    
    tau_prey[j] ~ normal(0, 10);
    corr_prey[j] ~ lkj_corr_cholesky(3);
    
    prey_means[j] ~ multi_normal_cholesky(preym[j], diag_pre_multiply(tau_mean, corr_mean));
    cons_prey[j] ~ multi_normal_cholesky(prey_means[j],prey_precs[j]);
    
    cs[j] ~ gamma(mean_c[j],tau_coeffs[j]);
    
    
    tau_prey_SI[j] ~ normal(0, 10);
    corr_prey_SI[j] ~ lkj_corr_cholesky(3);
    
    prey_means_SI[j] ~ multi_normal_cholesky(preym_SI[j], diag_pre_multiply(tau_mean_SI, corr_mean_SI));
    cons_prey_SI[j] ~ multi_normal_cholesky(prey_means_SI[j],prey_precs_SI[j]);
    
    cc[j] ~ normal(mean_cs[j],sigma_cs[j]);
    
  }
  
  
  // priors for predator covariance
  tau_pred ~ normal(0, 10);
  corr_pred ~ lkj_corr_cholesky(3);
  tau_pred_SI ~ normal(0, 10);
  corr_pred_SI ~ lkj_corr_cholesky(3);
  
  tau_mean ~ normal(0, 10);
  corr_mean ~ lkj_corr_cholesky(3);
  tau_mean_SI ~ normal(0, 10);
  corr_mean_SI ~ lkj_corr_cholesky(3);
  
}
