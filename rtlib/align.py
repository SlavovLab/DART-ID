#!/usr/bin/env python3
# coding: utf-8

import argparse
import logging
import numpy as np
import os
import pandas as pd
import pickle
import pkg_resources
import pystan
import re
import time

from rtlib.converter import process_files
from rtlib.helper import *
from hashlib import md5
from scipy.stats import norm, lognorm

pd.options.mode.chained_assignment = None
logger = logging.getLogger('root')

dirname = os.path.dirname(__file__)
print(dirname)

def align(dfa, config):

  # take subset of confident observations to use for alignment
  dff = dfa[-(dfa['exclude'])]
  dff = dff.reset_index(drop=True)

  logger.info('{} / {} ({:.2%}) confident, alignable observations (PSMs) after filtering.'.format(dff.shape[0], dfa.shape[0], dff.shape[0] / dfa.shape[0]))

  # refactorize peptide id into stan_peptide_id, 
  # to preserve continuity when feeding data into STAN
  dff['stan_peptide_id'] = dff['sequence'].map({ind: val for val, ind in enumerate(dff['sequence'].unique())})

  exp_names = np.sort(dff['raw_file'].unique())
  num_experiments = len(exp_names)
  num_observations = dff.shape[0]
  num_peptides = dff['stan_peptide_id'].max() + 1
  
  logger.info('Building peptide-experiment pairs...')

  # build a unique peptide-experiment ID
  pep_exp_all = dff['stan_peptide_id'].map(str) + ' - ' + dff['exp_id'].map(str)
  pep_exp_pairs = pep_exp_all.unique()
  num_pep_exp_pairs = len(pep_exp_pairs)

  # muij_map - maps pair ID back to dataframe index
  muij_map = pep_exp_all.map({ind: val for val, ind in enumerate(pep_exp_pairs)})
  # maps to experiment and peptide ID
  splt = pd.Series(pep_exp_pairs).str.split(' - ')
  muij_to_pep = splt.str.get(0).map(int)
  muij_to_exp = splt.str.get(1).map(int)

  logger.info('{} peptide-experiment pairs.'.format(num_pep_exp_pairs))

  # build dictionary of data to feed to STAN
  # revert all data to primitive types to avoid problems later
  # STAN code is all 1-indexed, so add 1 to all indexed forms of data
  stan_data = {
    'num_experiments': num_experiments,
    'num_peptides': num_peptides,
    'num_pep_exp_pairs': num_pep_exp_pairs,
    'num_total_observations': num_observations,
    'muij_map': muij_map+1,
    'muij_to_pep': (muij_to_pep+1).tolist(),
    'muij_to_exp': (muij_to_exp+1).tolist(),
    'experiment_id': (dff['exp_id']+1).tolist(),
    'peptide_id': (dff['stan_peptide_id']+1).tolist(),
    'retention_times': dff['retention_time'].tolist(),
    'mean_log_rt': np.mean(np.log(dff['retention_time'])),
    'sd_log_rt': np.std(np.log(dff['retention_time'])),
    'rt_mean': np.mean(dff['retention_time']),
    'rt_std': np.std(dff['retention_time']),
    'pep': dff['pep'].tolist(),
    'max_retention_time': dff['retention_time'].max()
  }

  logger.info('Initializing fit priors for {} peptides...'.format(num_peptides))

  # get the PEP filter, if it exists
  filter_pep = next((f for f in config['filters'] if f['name'] == 'pep'))
  if filter_pep is not None: filter_pep = filter_pep['value']
  else:                      filter_pep = np.max(df['pep'])

  # get the average retention time for a peptide, weighting by PEP
  def get_mu(x):
      weights = ((1 - x['pep'].values) - (1 - filter_pep)) / filter_pep
      return np.sum(x['retention_time'].values * weights) / np.sum(weights)
  
  # apply the get_mu function on all peptides and add some distortion
  mu_init = dff.groupby('stan_peptide_id')[['pep', 'retention_time']].apply(get_mu).values + np.random.normal(0, config['rt_distortion'], num_peptides)

  # negative or very low retention times not allowed. floor at 5 minutes
  mu_init[mu_init <= config['mu_min']] = config['mu_min']
  # canonical retention time shouldn't be bigger than largest real RT
  mu_max = dff['retention_time'].max()
  mu_init[mu_init > mu_max] = mu_max

  # take retention times and distort
  if config['rt_distortion'] > 0:
    logger.info('Distorting RTs by {} minutes for initial value generation'.format(config['rt_distortion']))

  rt_distorted = dff['retention_time'] + np.random.normal(0, config['rt_distortion'], len(dff['retention_time']))
  # make sure distorted retention times stay within bounds of real ones
  rt_distorted[rt_distorted > dff['retention_time'].max()] = dff['retention_time'].max()
  rt_distorted[rt_distorted < dff['retention_time'].min()] = dff['retention_time'].min()

  # initialize priors for the segmented linear regression
  # first element of vector is beta_0, or the intercept
  # second element is beta_1 and beta_2, the slopes of the two segments
  beta_init = np.array((
    np.repeat(10, num_experiments), 
    np.repeat(1, num_experiments)), dtype=float)

  logger.info('Optimizing priors with linear approximation for {} iterations.'.format(config['prior_iters']))

  mu_pred = np.zeros(num_peptides)
  # temporary data frame to quickly map over in the loop
  dft = pd.DataFrame(dict(
    stan_peptide_id=dff['stan_peptide_id'], 
    exp_id=dff['exp_id'], 
    pep=dff['pep'], 
    retention_time=mu_init[dff['stan_peptide_id']]))

  for i in range(0, config['prior_iters']):
    # for each experiment, fit a simple linear regression
    # between the distorted RTs and the initial canonical retention times
    for j in range(0, num_experiments):
        idx     = (dff['exp_id'] == j)
        rt_cur  = rt_distorted[idx]
        mu_cur  = mu_init[dff['stan_peptide_id'][idx]]
        pep_cur = dff['pep'][idx]

        # for this experiment, run a linear regression (1 degree polyfit)
        # of the mus and the distorted RTs. store the linear regression params
        m, c = np.polyfit(mu_cur, rt_cur, 1, w=(1 - pep_cur))
        beta_init[(0,1), j] = [c, m]

    # calculate new set of canonical RTs based on linear regression params
    mu_pred = (rt_distorted - beta_init[0][dff['exp_id']]) / beta_init[1][dff['exp_id']] 
    # make sure new canonical RTs are within same range as distorted RTs
    mu_pred[mu_pred <= 0] = rt_distorted.min()
    mu_pred[mu_pred >= rt_distorted.max()] = rt_distorted.max()
    dft['retention_time'] = np.copy(mu_pred)

    mu_prev = np.copy(mu_init)

    # new set of priors for canonical RTs based on weighted combination of
    # this set of predicted canonical RTs
    mu_init = dft.groupby('stan_peptide_id')[['pep', 'retention_time']].apply(get_mu).values

    logger.info('Iter {} | Avg. canonical RT shift: {:.5f}'.format(i + 1, pow(np.sum(mu_prev - mu_init), 2) / len(mu_init)))

  # grab linear regression params
  # set beta_2 (slope of second segment) to the same as the slope of the 1st segment
  beta_0 = beta_init[0]
  beta_1 = beta_init[1]
  beta_2 = np.copy(beta_1)

  # apply lower bound of (-1.0 * min(beta_1) * min(muInit)) to beta_0
  # where (-1.0 * min(beta_1) * min(muInit)) is the lowest possible intercept
  # given the lowest possible mu and lowest possible beta_1
  beta_0[beta_0 <= (-1 * beta_1.min() * mu_init.min())] = (-1 * beta_1.min() * mu_init.min()) + 1e-3
  # apply lower bound of 0 to beta_1 and beta_2
  beta_1[beta_1 <= 0] = 1e-3
  beta_2[beta_2 <= 0] = 1e-3

  # apply upper bound to prior canonical RTs
  mu_init[mu_init >= dff['retention_time'].max()] = 0.95 * dff['retention_time'].max()

  # init sigmas (spread) for each peptide
  sigma_init = np.zeros(num_peptides)
  muijs = beta_init[0][dff['exp_id']] + (beta_init[1][dff['exp_id']] * mu_init[dff['stan_peptide_id']])
  for i in range(0, num_peptides):
    sigma_init[i] = np.std(muijs[dff['stan_peptide_id']==i])

  # add sigma statistics to data list, to build sampling dist off of
  stan_data['sigma_mean'] = np.mean(sigma_init)
  stan_data['sigma_std'] = np.std(sigma_init)

  # create prior list for STAN
  init_list = {
    'mu': mu_init,
    'sigma': sigma_init,
    'beta_0': beta_0,
    'beta_1': beta_1,
    'beta_2': beta_2,
    'sigma_slope': np.repeat(0.1, num_experiments),
    'sigma_intercept': np.repeat(0.1, num_experiments),
    'split_point': np.repeat(np.median(mu_init), num_experiments)
  }

  # run STAN, store optimization parameters
  sm = StanModel_cache()

  logger.info('Running STAN for {} iterations and {} attempts in case of failure'.format(config['stan_iters'], config['stan_attempts']))

  # sometimes STAN will error out due to bad RNG or bad priors
  # set a limit on how many times we will try this stan configuration 
  # before running it all the way back again
  counter = 1
  op = None
  while op is None and counter <= config['stan_attempts']:
    try:
      # run STAN and time it
      logger.info('Starting STAN Model | Attempt #{} ...'.format(counter))
      start = time.time()
      op = sm.optimizing(data=stan_data, init=init_list, iter=config['stan_iters'], verbose=config['verbose'])
      logger.info('STAN Model Finished. Run time: {:.3f} seconds.'.format(time.time() - start))
    except RuntimeError as e:
      logger.error(str(e))
      counter = counter + 1

  # if loop terminates without any optimization parameters, then STAN failed
  if op is None:
    raise Exception('Maximum number of tries exceeded for STAN. Please re-run process or choose different parameters.')

  # exp_params - contains regression parameters for each experiment
  # peptide_params - contains canonical RT (mu) for each peptide sequence
  # pair_params - contains muij and sigmaij for each peptide-experiment pair
  #               not exactly necessary as these can be extrapolated from 
  #               exp_params and peptide_params
  # 
  exp_params = pd.DataFrame({ key: op[key] for key in ['beta_0', 'beta_1', 'beta_2', 'split_point', 'sigma_slope', 'sigma_intercept']})
  peptide_params = pd.DataFrame({ key: op[key] for key in ['mu']})
  pair_params = pd.DataFrame({ key: op[key] for key in ['muij', 'sigma_ij']})

  # add exp_id to exp_params
  exp_params['exp_id'] = np.sort(dff['exp_id'].unique())

  # put the maps in the pair parameters file as well
  pair_params = pair_params.assign(muij_to_pep=muij_to_pep)
  pair_params = pair_params.assign(muij_to_exp=muij_to_exp)

  if config['save_params']:
    # write parameters to file, so operations can be done on alignment data without
    # the entire alignment to run again
    logger.info('Writing STAN parameters to file...')
    exp_params.to_csv(os.path.join(config['output'], 'exp_params.txt'), sep='\t', index=False)
    pair_params.to_csv(os.path.join(config['output'], 'pair_params.txt'), sep='\t', index=True, index_label='pair_id')
    peptide_params.to_csv(os.path.join(config['output'], 'peptide_params.txt'), sep='\t', index=True, index_label='peptide_id')

  # write parameters to dict for further operations
  params = {}
  params['exp'] = exp_params
  params['pair'] = pair_params
  params['peptide'] = peptide_params

  return params

# taken from: https://pystan.readthedocs.io/en/latest/avoiding_recompilation.html#automatically-reusing-models
def StanModel_cache():

  model_name = 'FitRT3D'
  model_code = pkg_resources.resource_string('rtlib', '/'.join(('fits', 'fit_RT3d.stan')))

  # convert from bytes to a string
  model_code = model_code.decode('utf-8')

  # Use just as you would `stan`
  code_hash = md5(model_code.encode('ascii')).hexdigest()
  cache_fn = os.path.join(dirname, 'fits/cached-model-{}-{}.pkl'.format(model_name, code_hash))

  try:
      # load cached model from file
      sm = pickle.load(open(cache_fn, 'rb'))
  except:
      # compile model
      logger.info('Compiling STAN Model {} ...'.format(model_name))
      sm = pystan.StanModel(model_name=model_name, model_code=model_code)
      with open(cache_fn, 'wb') as f:
          # save model to file
          pickle.dump(sm, f)
  else:
      logger.info('Using cached StanModel: {}_{}'.format(model_name, code_hash))

  return sm


def main():
  # load command-line args
  parser = argparse.ArgumentParser()  
  add_global_args(parser)
  args = parser.parse_args()

  # load config file
  # this function also creates the output folder
  config = read_config_file(args)

  # initialize logger
  init_logger(config['verbose'], os.path.join(config['output'], 'align.log'))

  logger.info('Converting files and filtering PSMs')
  df, df_original = process_files(config)
  logger.info('Finished converting files and filtering PSMs.')

  logger.info('Beginning alignment procedure')
  align(df, config)
  logger.info('Alignment procedure finished')

if __name__ == '__main__':
  main()


