---
layout: default
title: Model Customization
nav_order: 8
permalink: docs/model-customization
---

# Alignment Model Customization

The DART-ID alignment model is written in [STAN](https://mc-stan.org/), a statistical modelling language. STAN allows the user to specify inputs, set priors, and assign a target likelihood to minimize. The full model is available [on GitHub](https://raw.githubusercontent.com/SlavovLab/DART-ID/master/dart_id/models/fit_RT3e.stan) and reproduced below:

```c
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
  //real<lower=0> rt_means[num_experiments];
  //real<lower=0> rt_stds[num_experiments];
  real<lower=0> rt_mean;
  real<lower=0> rt_std;
  real<lower=0, upper=1> pep[num_total_observations];

  real<lower=0> max_retention_time;
}
parameters {
  // canonical retention time for given peptide
  real<lower=0, upper=max_retention_time> mu[num_peptides];

  // segmented linear regression parameters for the 
  // retention time of each pep-exp pair
  //
  // beta_0 = y-intercept
  // beta_1 = slope of first segment
  // beta_2 = slope of second segment
  // split_point = x-value of the split b/n the two segments
  real<lower=0> beta_1[num_experiments];
  real<lower=0> beta_2[num_experiments];
  real<lower= -1.0 * min(beta_1) * min(mu) > beta_0[num_experiments];
  real<lower=0> split_point[num_experiments];

  // linear regression parameters for the standard deviation/spread 
  // for each pep-exp pair
  //
  // sigma_slope_global = average sigma slope for the set of pep-exp pairs.
  //                      "moving prior" for sigma_slopes of pairs moving forwards
  // sigma_slope = slope of line
  // sigma_intercept = y-intercept of line - i.e., the smallest sigma_ij possible in that experiment
  real<lower=0> sigma_slope_global;
  real<lower=0> sigma_slope[num_experiments];
  real<lower=0.01> sigma_intercept[num_experiments];

}
transformed parameters {
  // alignment is based on a segemented linear regression
  
  // sigma_ij is the spread of the laplace distribution of 
  // peptide i in experiment j
  real<lower=0.01> sigma_ij[num_pep_exp_pairs];
  
  // muij is the mean retention time for peptide i in experiment j 
  real<lower=0, upper=max_retention_time> muij[num_pep_exp_pairs];
  
  for (i in 1:num_pep_exp_pairs) {
    if(mu[muij_to_pep[i]] < split_point[muij_to_exp[i]]) {
      // first segment of the linear regression (before split_point)
      muij[i] = beta_0[muij_to_exp[i]] + 
        beta_1[muij_to_exp[i]] * mu[muij_to_pep[i]];
    } else if( mu[muij_to_pep[i]] >=  split_point[muij_to_exp[i]]) {
      // second segment of the linear regression (after split_point)
      muij[i] = beta_0[muij_to_exp[i]] + 
        beta_1[muij_to_exp[i]] * split_point[muij_to_exp[i]] + 
        beta_2[muij_to_exp[i]] * (mu[muij_to_pep[i]] - split_point[muij_to_exp[i]]);
    }
    
    // simple linear regression
    //muij[i] = beta_0[muij_to_exp[i]] + beta_1[muij_to_exp[i]] * mu[muij_to_pep[i]];

    // make sure that the retention time is within the physical limits
    if( muij[i] > max_retention_time ) {
      muij[i] = max_retention_time * 0.99;
    }
  }

  // standard deviation grows linearly with time - simple linear regression
  for(i in 1:num_pep_exp_pairs) { 
    sigma_ij[i] = sigma_intercept[muij_to_exp[i]] + 
      sigma_slope[muij_to_exp[i]] / 100 * mu[muij_to_pep[i]];
  }

}
model {
  mu ~ normal(rt_mean, rt_std);

  sigma_slope_global ~ lognormal(0.1, 0.5);
  sigma_slope ~ lognormal(log(sigma_slope_global), 1);
  sigma_intercept ~ lognormal(0, 2);
  
  beta_0 ~ normal(0, 10);
  beta_1 ~ lognormal(0, 0.5);
  beta_2 ~ lognormal(0, 0.5);

  //for(i in 1:num_experiments){
  split_point ~ uniform(0, max_retention_time);
  //}
  
  for(i in 1:num_total_observations) {
    target += log_mix(pep[i],
      normal_lpdf(retention_times[i] | rt_mean, rt_std),
      double_exponential_lpdf(retention_times[i] | muij[muij_map[i]], sigma_ij[muij_map[i]]));
  }
}
```

## Data

The first `data` block specifies the inputs to the model. Most variables are used for mapping values back and forth. If you need to pass in more data to the model to, for example, change the prior distributions, you need to specify them here and pass them into the model from the code in `align.py` as well.

## Parameters

This section is for defining the parameters: the reference retetion time `mu` and the experiment fits. 

`beta_0`, `beta_1`, `beta_2`, and `split_point` are the necessary parameters to define the two-segment linear model. If you wish to change the model, you will have to define your new parameters here.

`sigma_slope_global`, `sigma_slope`, `sigma_intercept` are used to model per-experiment RT variances as well as RT-dependent shifts in variances. Possible extensions would be to define a fixed per-experiment RT variance, or a global RT variance, or even a per-peptide RT variance.

## Transformed Parameters

This section derives the peptide-experiment predicted retention time `muij` and variance `sigma_ij` from the parameters in the previous block. This section defines the logic behind the two-segment linear model for the reference RT and the linear model for the variance.

## Model

This section defines the priors for the parameters. New experimental designs, such as longer gradients, a mix of gradient lengths, or a restriction in the data (RT in seconds instead of minutes) may necessitate a change in the prior distributions.

The final piece is the "objective function", which is represented here as `target` being the log likelihood of the model over all observations.

## Implementing the Model

We recommend that you test and validate any customized model to ensure that it is returning sensible results. 

In order to interface the model with the DART-ID python code, you must compile the `.stan` model code with [cmdstan](https://mc-stan.org/users/interfaces/cmdstan). Link the new model to a model definition by adding the definition in `models.py`, and ensure that the new model binary is placed in the right location. (Will write up more detailed guide soon...)

