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
  real<lower=0, upper=1> pep[num_total_observations];

  real<lower=0> max_retention_time;
}
parameters {

  // canonical retention time for given peptide
  real<lower=0, upper=max(retention_times)> mu[num_peptides];


  real<lower=0> beta_1[num_experiments];
  //real<lower=0> beta_2[num_experiments];
  real beta_0[num_experiments];
  //real<lower=0> split_point[num_experiments];

  
  real<lower=0> sigma_slope_global;
  real<lower=0> sigma_slope[num_experiments];
  real<lower=0> sigma_intercept[num_experiments];

}
transformed parameters {

  real<lower=0> sigma_ij[num_pep_exp_pairs];
  
  // alignment is based on a segemented linear regression

  // muij is the mean retention time for peptide i in experiment j 
  real<lower=0, upper=max_retention_time> muij[num_pep_exp_pairs];
  for (i  in 1:num_pep_exp_pairs) {

    muij[i] = max_retention_time*inv_logit(beta_1[muij_to_exp[i]] * (mu[muij_to_pep[i]] - beta_0[muij_to_exp[i]]));

  }

  // standard deviation grows linearly with time
  for(i in 1:num_pep_exp_pairs) { 
    sigma_ij[i] = sigma_intercept[muij_to_exp[i]] + 
      sigma_slope[muij_to_exp[i]] / 100 * mu[muij_to_pep[i]];
  }

}
model {
  mu ~ lognormal(mean(log(retention_times)), sd(log(retention_times)));

  sigma_slope_global ~ lognormal(0.1, 0.5);
  //sigma_slope ~ lognormal(log(sigma_slope_global), 0.5);
  sigma_slope ~ lognormal(log(sigma_slope_global), 1);
  sigma_intercept ~ lognormal(0, 2);
  
  beta_0 ~ normal(max_retention_time/2, 60);
  beta_1 ~ exponential(1); //lognormal(-3, 1);

  for (i in 1:num_total_observations) {

    //mixture between laplace and a uniform prior (distribution over all peptides)
    real comp1 = log(1-pep[i]) + double_exponential_lpdf(retention_times[i] | muij[muij_map[i]], sigma_ij[peptide_id[i]]);
    real comp2 = log(pep[i]) + uniform_lpdf(retention_times[i] | 0, max_retention_time);

    real lse = log_sum_exp(comp1, comp2); 

    target += lse;

  }


}
