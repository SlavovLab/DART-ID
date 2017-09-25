// -*- mode: C -*-

data {

  int<lower=1> num_experiments;
  int<lower=1> num_peptides;
  int<lower=1> num_total_observations;
  int<lower=1> num_pep_exp_pairs;

  int<lower=1> muij_map[num_total_observations];
  int<lower=1> muij_to_exp[num_pep_exp_pairs];
  int<lower=1> muij_to_pep[num_pep_exp_pairs];
  
  int<lower=1, upper=num_experiments> experiment_id[num_total_observations];
  int<lower=1, upper=num_peptides> peptide_id[num_total_observations];
  real<lower=0> retention_times[num_total_observations];

  
  
}
parameters {

  // canonical retention time for given peptide
  real<lower=0> mu[num_peptides];

  real beta_0[num_experiments];
  real<lower=0> beta_1[num_experiments];
  real beta_2[num_experiments];
  real<lower=0> split_point[num_experiments];

  real<lower=0> sigma_global;
  real<lower=0> sigma[num_peptides];

  //real<lower=0> sigma_shrinkage;
}
transformed parameters {

  ## alignment is based on a segemented linear regression

  ## muij is the mean retention time for peptide i in experiment j 
  real muij[num_pep_exp_pairs];
  for (i in 1:num_pep_exp_pairs) {
    if(mu[muij_to_pep[i]] < split_point[muij_to_exp[i]]) {
      muij[i] = beta_0[muij_to_exp[i]] + beta_1[muij_to_exp[i]] * mu[muij_to_pep[i]];
    } else if( mu[muij_to_pep[i]] >=  split_point[muij_to_exp[i]] ) {
      muij[i] = beta_0[muij_to_exp[i]] + beta_1[muij_to_exp[i]] * split_point[muij_to_exp[i]] + beta_2[muij_to_exp[i]] * (mu[muij_to_pep[i]] - split_point[muij_to_exp[i]]);
      //muij[i] = beta_0[muij_to_exp[i]] + beta_1[muij_to_exp[i]] * mu[muij_to_pep[i]];
    
    }
  }
}
model {
  mu ~ lognormal(mean(log(retention_times)), sd(log(retention_times)));
  split_point ~ lognormal(mean(log(retention_times)), sd(log(retention_times)));
  sigma_global ~ lognormal(0, 10);
  //sigma_shrinkage ~ lognormal(0, 10);
  
  beta_0 ~ normal(0, 1);
  beta_1 ~ lognormal(0, 0.5);
  beta_2 ~ lognormal(0, 0.5);
  
  for(p in 1:num_peptides) {
    sigma[p] ~ lognormal(log(sigma_global), 0.2);
  }
  
  for (i in 1:num_total_observations) {
    retention_times[i] ~ normal(muij[muij_map[i]], sigma[peptide_id[i]]);
  }
}
