#!/usr/bin/env python3
# coding: utf-8

import argparse
import json
import logging
import numpy as np
import os
import pandas as pd
import pickle
import pkg_resources
import platform
#import pystan
import time
import sys

#sys.path.append('/Users/albert/git/pystan')
import pystan


from dart_id.converter import process_files
from dart_id.exceptions import *
from dart_id.models import models, get_model_from_config
from dart_id.helper import *
from distutils.sysconfig import get_config_var
from distutils.version import LooseVersion
from hashlib import md5
from scipy.stats import norm, lognorm

pd.options.mode.chained_assignment = None
logger = logging.getLogger('root')

dirname = os.path.dirname(__file__)

def align(dfa, config):

  # take subset of confident observations to use for alignment
  dff = dfa[-(dfa['exclude'])].reset_index(drop=True)

  #logger.info('{} / {} ({:.2%}) confident, alignable observations (PSMs) after filtering.'.format(dff.shape[0], dfa.shape[0], dff.shape[0] / dfa.shape[0]))

  # refactorize peptide id into stan_peptide_id, 
  # to preserve continuity when feeding data into STAN
  dff['stan_peptide_id'] = dff['sequence'].map({
    ind: val for val, ind in enumerate(dff['sequence'].unique())})

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
    'num_peptides': num_peptides.tolist(),
    'num_pep_exp_pairs': num_pep_exp_pairs,
    'num_total_observations': num_observations,
    'muij_map': (muij_map+1).tolist(),
    'muij_to_pep': (muij_to_pep+1).tolist(),
    'muij_to_exp': (muij_to_exp+1).tolist(),
    'experiment_id': (dff['exp_id']+1).tolist(),
    'peptide_id': (dff['stan_peptide_id']+1).tolist(),
    'retention_times': dff['retention_time'].tolist(),
    'mean_log_rt': np.mean(np.log(dff['retention_time'])).tolist(),
    'sd_log_rt': np.std(np.log(dff['retention_time'])).tolist(),
    'rt_mean': np.mean(dff['retention_time']).tolist(),
    'rt_std': np.std(dff['retention_time']).tolist(),
    'pep': dff['pep'].tolist(),
    'max_retention_time': dff['retention_time'].max().tolist()
  }

  logger.info('Initializing fit priors for {} peptides...'.format(num_peptides))

  model = get_model_from_config(config)

  logger.info('Aligning RTs using the \"{}\" model'.format(model['name']))

  # generate initial values with the given function in models.py
  init_list = model['init_func'](dff, config)

  # init sigmas (spread) for each peptide
  #sigma_init = np.zeros(num_peptides)
  #muijs = init_list['beta_0'][dff['exp_id']] + \
  #  (init_list['beta_1'][dff['exp_id']] * init_list['mu'][dff['stan_peptide_id']])
  #for i in range(0, num_peptides):
  #  sigma_init[i] = np.std(muijs[dff['stan_peptide_id']==i])

  # add sigma statistics to data list, to build sampling dist off of
  #stan_data['sigma_mean'] = np.mean(sigma_init)
  #stan_data['sigma_std'] = np.std(sigma_init)
  
  # allow user to save STAN input data to separate files, so they can run STAN standalone
  # (with cmdstan). Useful for debugging alignment errors.
  if 'save_stan_input' in config and config['save_stan_input']:
    logger.info('Saving STAN input data, for standalone STAN runs, in JSON format.')

    with open(os.path.join(config['output'], 'stan_input.json'), 'w') as f:
      logger.info('Saving input data to stan_input.json')
      f.write(json.dumps(stan_data))
    
    with open(os.path.join(config['output'], 'stan_init_list.json'), 'w') as f:
      logger.info('Saving initial values to stan_init_list.json')
      f.write(json.dumps(init_list))

  # run STAN, store optimization parameters
  sm = StanModel_cache(model)

  """
  logger.info('Running STAN for {} iterations and {} attempt(s) in case of failure'.format(config['stan_iters'], config['stan_attempts']))

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

      op = sm.optimizing(data=stan_data, init=init_list, as_vector=False,
        iter=config['stan_iters'], verbose=config['verbose'],
        save_iterations=True, tol_param=1e-8)
      print(op)

      logger.info('STAN Model Finished. Run time: {:.3f} seconds.'.format(
        time.time() - start))

    except RuntimeError as e:
      logger.error(str(e))
      counter = counter + 1

  # if loop terminates without any optimization parameters, then STAN failed
  if op is None:
    raise STANError('Maximum number of tries exceeded for STAN. Please re-run process or choose different parameters.')
  """

  # run STAN and time it
  logger.info('Starting STAN Model...')
  start = time.time()

  op = sm.optimizing(data=stan_data, init=init_list, as_vector=False,
    iter=config['stan_iters'], verbose=config['verbose'])
  print(op)

  if op is None:
    raise STANError('STAN did not return any parameters. Please check that pystan was installed correctly.')

  logger.info('STAN Model Finished. Run time: {:.3f} seconds.'.format(
    time.time() - start))

  # exp_params - contains regression parameters for each experiment
  # peptide_params - contains canonical RT (mu) for each peptide sequence
  # pair_params - contains muij and sigmaij for each peptide-experiment pair
  #               not exactly necessary as these can be extrapolated from 
  #               exp_params and peptide_params
  
  pars = op['par']
  exp_params = pd.DataFrame({ 
    key: pars[key] for key in model['exp_keys']})
  peptide_params = pd.DataFrame({ key: pars[key] for key in model['peptide_keys']})
  pair_params = pd.DataFrame({ key: pars[key] for key in model['pair_keys']})

  # add exp_id to exp_params
  exp_params['exp_id'] = np.sort(dff['exp_id'].unique())

  # put the maps in the pair parameters file as well
  pair_params = pair_params.assign(muij_to_pep=muij_to_pep)
  pair_params = pair_params.assign(muij_to_exp=muij_to_exp)

  # add initial values to parameters files, and prepend those column names with 'init_'
  exp_params = pd.concat([exp_params,
    pd.DataFrame({ 'init_'+key: init_list[key] for key in model['exp_keys']})],
    axis=1)
  peptide_params = pd.concat([peptide_params,
    pd.DataFrame({ 'init_'+key: init_list[key] for key in model['peptide_keys']})],
    axis=1)
  #pair_params = pd.concat([pair_params,
  #  pd.DataFrame({ 'init_'+key: init_list[key] for key in model['pair_keys']})],
  #  axis=1, ignore_index=True)

  if config['save_params']:
    # write parameters to file, so operations can be done on alignment data without
    # the entire alignment to run again
    logger.info('Writing STAN parameters to file...')
    exp_params.to_csv(os.path.join(config['output'], 'exp_params.txt'), 
      sep='\t', index=False)
    pair_params.to_csv(os.path.join(config['output'], 'pair_params.txt'), 
      sep='\t', index=True, index_label='pair_id')
    peptide_params.to_csv(os.path.join(config['output'], 'peptide_params.txt'), 
      sep='\t', index=True, index_label='peptide_id')

  # write parameters to dict for further operations
  params = {}
  params['exp'] = exp_params
  params['pair'] = pair_params
  params['peptide'] = peptide_params

  return params

# taken from: https://pystan.readthedocs.io/en/latest/avoiding_recompilation.html
# #automatically-reusing-models
def StanModel_cache(model):

  model_name = model['model_name']
  model_code = pkg_resources.resource_string('dart_id', '/'.join((
    'models', model['stan_file'])))
   
  # convert from bytes to a string
  model_code = model_code.decode('utf-8')

  # Use just as you would `stan`
  code_hash = md5(model_code.encode('ascii')).hexdigest()
  cache_fn = os.path.join(dirname, 'models/cached-model-{}-{}.pkl'.format(
    model_name, code_hash))
  
  try:
      # load cached model from file
      sm = pickle.load(open(cache_fn, 'rb'))
  except:
      # compile model
      logger.info('Compiling STAN Model {} ...'.format(model_name))

      # https://github.com/pandas-dev/pandas/commit/459ebb2ad309938e9e68ae79e3bdb312efac0ca2
      # For mac, ensure extensions are built for macos 10.9 when compiling on a
      # 10.9 system or above, overriding distuitls behaviour which is to target
      # the version that python was built for. This may be overridden by setting
      # MACOSX_DEPLOYMENT_TARGET before calling setup.py
      if sys.platform == 'darwin':
        if 'MACOSX_DEPLOYMENT_TARGET' not in os.environ:
          current_system = LooseVersion(platform.mac_ver()[0])
          python_target = LooseVersion(
            get_config_var('MACOSX_DEPLOYMENT_TARGET'))
          if python_target < '10.9' and current_system >= '10.9':
            os.environ['MACOSX_DEPLOYMENT_TARGET'] = '10.9'

      sm = pystan.StanModel(model_name=model_name, model_code=model_code)
      with open(cache_fn, 'wb') as f:
          # save model to file
          pickle.dump(sm, f)
  else:
      logger.info('Using cached StanModel: {}_{}'.format(model_name, code_hash))
  
  return sm

"""
  logger.info('Compiling STAN Model {} ...'.format(model_name))
  sm = pystan.StanModel(model_name=model_name, model_code=model_code)
  with open(cache_fn, 'wb') as f:
      # save model to file
      pickle.dump(sm, f)
"""

def main():
  # load command-line args
  parser = argparse.ArgumentParser()  
  add_global_args(parser)
  args = parser.parse_args()

  # load config file
  # this function also creates the output folder
  config = read_config_file(args)

  # initialize logger
  init_logger(config['verbose'], os.path.join(config['output'], 'align.log'), config['log_file'])

  logger.info('Converting files and filtering PSMs')
  df, df_original = process_files(config)
  logger.info('Finished converting files and filtering PSMs.')

  logger.info('Beginning alignment procedure')
  align(df, config)
  logger.info('Alignment procedure finished')

if __name__ == '__main__':
  main()


