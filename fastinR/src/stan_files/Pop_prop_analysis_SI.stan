data {
  int<lower = 1> n_prey_samps_SI;
  int<lower = 1> n_preys;
  int<lower = 1> n_preds;
  int<lower = 1> isos;
  int prey_ix_SI[n_prey_samps_SI];
  vector[isos] preys_SI[n_prey_samps_SI];
  vector[isos] preym_SI[n_preys];
  vector[isos] mean_cs[n_preys];
  vector[isos] sigma_cs[n_preys];
  vector[isos] preds_SI[n_preds];
}
parameters{
  vector[isos] prey_means_SI[n_preys];
  vector[isos] cc[n_preys];
  vector[isos] cons_prey_SI[n_preys];
  
  cholesky_factor_corr[isos] corr_prey[n_preys];
  cholesky_factor_corr[isos] corr_pred;
  cholesky_factor_corr[isos] corr_mean;
  vector<lower = 0>[isos] tau_prey[n_preys];
  vector<lower = 0>[isos] tau_mean;
  vector<lower = 0>[isos] tau_pred;
  real<lower = 0> props[n_preys];
}
transformed parameters{
  matrix[isos,isos] prey_precs_SI[n_preys];
  matrix[isos,isos] pred_prec_SI;
  vector[n_preys] prop;
  matrix[isos,n_preys] prey_SI;
  vector[isos] mu_SI;
  
  for (j in 1:n_preys){
    prop[j] = props[j]/sum(props[1:n_preys] );
    prey_precs_SI[j] = diag_pre_multiply(tau_prey[j], corr_prey[j]);
    prey_SI[,j] = cc[j]+cons_prey_SI[j];
      // mixing for predator SI signature for likelihood
  }
  //for (i in 1:isos){
  mu_SI = prey_SI*prop;
  
  pred_prec_SI = diag_pre_multiply(tau_pred, corr_pred);
  
}
model{

//prey samples
  for(ps in 1:n_prey_samps_SI)  preys_SI[ps] ~ multi_normal_cholesky(prey_means_SI[prey_ix_SI[ps]], prey_precs_SI[prey_ix_SI[ps]]);

// predator likelihood
    
  for (p in 1:n_preds) preds_SI[p] ~ multi_normal_cholesky(mu_SI,pred_prec_SI);

  props ~ gamma(1.5,1);    

  for (j in 1:n_preys){
    
    tau_prey[j] ~ normal(0, 10);
    corr_prey[j] ~ lkj_corr_cholesky(3);
    
    prey_means_SI[j] ~ multi_normal_cholesky(preym_SI[j], diag_pre_multiply(tau_mean, corr_mean));
    cons_prey_SI[j] ~ multi_normal_cholesky(prey_means_SI[j],prey_precs_SI[j]);

          
    cc[j] ~ normal(mean_cs[j],sigma_cs[j]);
    
    
  }

      
  // priors for predator covariance
  
   
   tau_pred ~ normal(0, 10);
   corr_pred ~ lkj_corr_cholesky(3);
  
   tau_mean ~ normal(0, 10);
   corr_mean ~ lkj_corr_cholesky(3);
 
  
}
