#!/usr/bin/env python3
# coding: utf-8

import argparse
import json
import logging
import numpy as np
import os
import pandas as pd
import pkg_resources
import platform
import re
try:
  import subprocess32 as subprocess
except ImportError:
  import subprocess
import sys
import time

from dart_id.converter import process_files
from dart_id.exceptions import STANError
from dart_id.models import models, get_model_from_config
from dart_id.helper import init_logger, add_global_args, read_config_file, convert_numpy_scalar

pd.options.mode.chained_assignment = None
logger = logging.getLogger('root')

def get_os():
  x = platform.system()
  if x == 'Linux': 
    # get the linux distro
    p = platform.platform()
    if re.search(r'centos-7', p) is not None or re.search(r'rhel', p) is not None:
      return 'rhel'
    elif re.search(r'debian', p) is not None:
      return 'debian'
    else:
      logger.warning('Linux distribution {} not recognized. Defaulting to generic linux tag...'.format(p))
      return 'linux_generic'

  elif x == 'Darwin': return 'mac'
  elif x == 'Windows': return 'windows'
  else: raise Exception('Operating system {} not identified'.format(x))

def get_exec_name(model):
  x = platform.system()
  if x == 'Windows':
    return '{}.exe'.format(model)
  else: return model

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
  #
  # convert any numpy objects to base python types/lists, so they can
  # be printed in non-binary form when this dict is converted to JSON
  stan_data = {
    'num_experiments': convert_numpy_scalar(num_experiments),
    'num_peptides': convert_numpy_scalar(num_peptides),
    'num_pep_exp_pairs': convert_numpy_scalar(num_pep_exp_pairs),
    'num_total_observations': convert_numpy_scalar(num_observations),
    'muij_map': (muij_map+1).tolist(),
    'muij_to_pep': (muij_to_pep+1).tolist(),
    'muij_to_exp': (muij_to_exp+1).tolist(),
    'experiment_id': (dff['exp_id']+1).tolist(),
    'peptide_id': (dff['stan_peptide_id']+1).tolist(),
    'retention_times': dff['retention_time'].tolist(),
    'mean_log_rt': convert_numpy_scalar(np.mean(np.log(dff['retention_time']))),
    'sd_log_rt': convert_numpy_scalar(np.std(np.log(dff['retention_time']))),
    'rt_mean': convert_numpy_scalar(np.mean(dff['retention_time'])),
    'rt_std': convert_numpy_scalar(np.std(dff['retention_time'])),
    'pep': dff['pep'].tolist(),
    'max_retention_time': convert_numpy_scalar(dff['retention_time'].max())
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
  
  
  logger.info('Preparing for STAN Alignment')

  logger.info('Saving data as JSON to pass to STAN executable...')

  stan_input_path = os.path.join(config['output'], 'stan_input.json')
  with open(stan_input_path, 'w') as f:
    logger.info('Saving input data to stan_input.json')
    f.write(json.dumps(stan_data))
    f.close()
  
  stan_init_list_path = os.path.join(config['output'], 'stan_init_list.json')
  with open(stan_init_list_path, 'w') as f:
    logger.info('Saving initial values to stan_init_list.json')
    f.write(json.dumps(init_list))
    f.close()

  exec_folder = get_os()
  logger.info('Detected operating system: {} ({} | {})'.format(exec_folder, platform.system(), platform.platform()))
  exec_name = get_exec_name(model['model_binary'])
  
  # find the model folder in models/
  model_exists = pkg_resources.resource_exists('dart_id', '/'.join(['models', model['model_binary'], exec_folder, exec_name]))
  if not model_exists:
    raise Exception('Compiled model \"{}\" does not exist for operating system \"{}\". Please contact the authors for the executable or build it locally on your machine with cmdstan.'.format(model['model_binary'], exec_folder))

  model_path = pkg_resources.resource_filename('dart_id', '/'.join(('models', model['model_binary'], exec_folder, exec_name)))
  logger.info('Found executable for model \"{}\" in: {}'.format(model['model_name'], model_path))


  stan_output_path = os.path.join(config['output'], 'stan_output.csv')

  cmd = [
    model_path,
    'optimize', 'iter={}'.format(config['stan_iters']),
    'algorithm=lbfgs',
      'init_alpha={}'.format(config['init_alpha']),
      'tol_obj={}'.format(config['tol_obj']),
      'tol_rel_obj={}'.format(config['tol_rel_obj']),
      'tol_grad={}'.format(config['tol_grad']),
      'tol_rel_grad={}'.format(config['tol_rel_grad']),
      'tol_param={}'.format(config['tol_param']),
      'history_size={}'.format(config['history_size']),
    'init={}'.format(stan_init_list_path),
    'data', 'file={}'.format(stan_input_path),
    'output', 'file={}'.format(stan_output_path)
  ]

  # run STAN and time it
  logger.info('Running STAN with command: {}'.format(' '.join(cmd)))
  start = time.time()

  #stan_log_path = os.path.join(config['output'], 'stan_alignment.log')
  #output_fd = open(stan_log_path, 'w')
  proc = subprocess.Popen(cmd)

  while True:
    try:
      proc.wait(1)
      break
    except subprocess.TimeoutExpired:
      continue
    except KeyboardInterrupt:
      break
  
  logger.info('STAN Model Finished. Run time: {:.3f} seconds.'.format(time.time() - start))

  #output_fd.close()

  #stdout = ''
  #with open(stan_log_path, 'r') as fd:
  #  stdout = fd.read()

  if proc.returncode != 0:
    logger.warning('Stan model exited with error {}'.format(proc.returncode))
  
  # parse the output file
  with open(stan_output_path, 'r') as f:
    stan_out = f.readlines()
    # we only need the parameters from the last iteration,
    # so just read the headers and the last line
    headers = stan_out[-2]
    values = stan_out[-1]
    f.close()

  # transform into a dict
  stan_results = dict()
  headers = headers.split(',')
  values = values.split(',')
  
  start = 0
  cur_header = headers[0].split('.')[0]
  for i, header in enumerate(headers):
    header = header.split('.')[0]

    if header != cur_header or i == (len(headers) - 1):
      # where to end the slice
      end = i
      if i == (len(headers) - 1): end = i + 1
      # store value slice for the previous header
      stan_results[cur_header] = np.array(values[start:end])
      # reset counter
      start = i
      # reset previous header
      cur_header = header
    else: continue

  if len(stan_results) == 0:
    raise STANError('STAN did not return any parameters, or the returned parameters could not be read. Please check your output folder for the file \"stan_output.csv\" and verify that it exists and holds the parameter data.')

  # convert strings to ints or floats
  for param in stan_results:
    is_float = '.' in stan_results[param][0]
    if is_float:
      stan_results[param] = stan_results[param].astype(float)
    else:
      stan_results[param] = stan_results[param].astype(int)

  # exp_params - contains regression parameters for each experiment
  # peptide_params - contains canonical RT (mu) for each peptide sequence
  # pair_params - contains muij and sigmaij for each peptide-experiment pair
  #               not exactly necessary as these can be extrapolated from 
  #               exp_params and peptide_params
  
  exp_params = pd.DataFrame({ 
    key: stan_results[key] for key in model['exp_keys']})
  peptide_params = pd.DataFrame({ key: stan_results[key] for key in model['peptide_keys']})
  pair_params = pd.DataFrame({ key: stan_results[key] for key in model['pair_keys']})

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


