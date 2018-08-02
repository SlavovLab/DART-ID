#!/usr/bin/env python3
# coding: utf-8

import argparse
import logging
import numpy as np
import os
import pandas as pd
import time

from dart_id.align import align
from dart_id.converter import process_files
from dart_id.exceptions import *
from dart_id.figures import figures
from dart_id.models import models, get_model_from_config
from dart_id.helper import *
from scipy.stats import norm, lognorm, laplace

logger = logging.getLogger('root')

def update(dfa, params, config):
  dfa = dfa.reset_index(drop=True)

  #logger.info('{} / {} ({:.2%}) confident, alignable observations (PSMs) after filtering.'.format(dff.shape[0], dfa.shape[0], dff.shape[0] / dfa.shape[0]))

  # refactorize peptide id into stan_peptide_id, 
  # to preserve continuity when feeding data into STAN
  dfa['stan_peptide_id'] = dfa['sequence'].map({
    ind: val for val, ind in enumerate(dfa['sequence'].unique())})

  num_experiments = dfa['exp_id'].max() + 1
  num_peptides = dfa['peptide_id'].max() + 1
  exp_names = np.sort(dfa['raw_file'].unique())
  pep_id_list = dfa['peptide_id'].unique()

  # validate parameters file. make sure it is from the same filters
  # or else the program will crash in the code below
  # check num_experiments, num_peptides
  if params['exp'].shape[0] != num_experiments or \
     params['peptide'].shape[0] != (dfa['stan_peptide_id'].max() + 1):
    raise ConfigFileError('Parameters files have different data than the input data provided. Ensure that both the input list and filters used to generate the alignment parameters and those provided to the current update are the __exact__ same.')
  
  model = get_model_from_config(config)

  # mu from the STAN alignment
  dfa['mu'] = params['peptide']['mu'].values[dfa['stan_peptide_id']] 
  
  # concatenate transformation parameters
  exp_params = pd.DataFrame({ key: params['exp'][key][dfa['exp_id']] \
    for key in model['exp_keys']}).reset_index(drop=True)
  dfa = pd.concat([dfa, exp_params], axis=1)

  # predict mus with RTs, and RTs with aligned mus
  dfa['mu_pred'] = model['rt_to_ref'](dfa, dfa['mu'], params)
  dfa['muij']    = model['ref_to_rt'](dfa, dfa['mu'], params)
  dfa['sigmaij'] = model['sigmaij_func'](dfa, params)

  # PEP ceiling at 1, otherwise will result in 
  # incorrect negative densities when plugging into Bayes' theorem
  dfa['pep'][dfa['pep'] > 1.0] = 1.0

  # output table
  df_new = pd.DataFrame()

  bootstrap_method = 'none'
  if 'bootstrap_method' in config:
    bootstrap_method = config['bootstrap_method']
    logger.info('Using \"{}\" bootstrap method'.format(bootstrap_method))
  else:
    logger.info('Bootstrap method not defined, using point estimates to update confidence instead.')

  k = 20 # default
  if 'bootstrap_iters' in config:
    k = config['bootstrap_iters']
    if bootstrap_method != 'none':
      logger.info('Using {} bootstrap iterations'.format(k))

  for i, e in enumerate(np.sort(dfa['exp_id'].unique())):
    exp_name = exp_names[i]
    logger.info('Updating PEPs: Exp ({} / {}) - {}'.format(i+1, num_experiments, exp_name))

    exp = dfa[dfa['exp_id'] == e]
    exp = exp.reset_index(drop=True)

    exp_peptides = exp['stan_peptide_id'].unique()

    # vector of P(RT|delta=1) for this experiment.
    rt_plus = pd.Series(np.zeros(exp.shape[0]))

    if bootstrap_method != 'none':

      # to avoid using this experiment's own data to update the confidence
      # of its own observations, recalculate the reference RTs (mu) without the
      # data from this experiment, by: 
      # 1) non-parametric bootstrapping over the median of the predicted mus.
      # OR
      # 2) parametric bootstrapping, using the RT distribution parameters
      
      # get predicted mus of peptides in this experiment, excluding predicted mus
      # transformed from RTs observed in this experiment
      dfe = dfa.loc[((dfa['stan_peptide_id'].isin(exp_peptides)) & (dfa['exp_id'] != e)), \
        ['stan_peptide_id', 'pep', 'mu_pred']]

      # extract predicted mus and PEPs    
      mu_preds = dfe.groupby('stan_peptide_id')['mu_pred']\
        .apply(lambda x: x.values).values.tolist()
      peps = dfe.groupby('stan_peptide_id')['pep']\
        .apply(lambda x: x.values).values.tolist()
      
      # number of observations per peptide sequence
      obs_per_seq = [len(peptide) for peptide in mu_preds]
      num_peptides = len(mu_preds)

      # the number of observations per peptide -- used in loop
      num_obs = 0
      # matrix of n by k estimated mus from the bootstrapping
      # will iterate over in the loop after the immediate one
      mu_k = np.zeros((num_peptides, k))

      if bootstrap_method == 'parametric':
        # parametric bootstrap
        for i in range(0, num_peptides):
          num_obs = obs_per_seq[i]

          # generate all random numbers in one go, and mold into matrix
          # where each column is a bootstrap iteration
          # then, take the median of each column (bootstrap iter)
          mu_k[i] = np.apply_along_axis(np.median, 0, 
            laplace.rvs(\
              loc=np.median(mu_preds[i]), 
              scale=np.std(mu_preds[i]), 
              size=k*num_obs)\
            .reshape(num_obs, k))

      #elif bootstrap_method == 'parametric_mixture':
        # generate random samples from the mixture model instead of just the
        # correct RT distribution
        # to do this, need to calculate then apply an FDR (sum of PEPs / # observations)
        # for each peptide.
        #for i in range(0, num_peptides):

      elif bootstrap_method == 'non-parametric':
        # non-parametric bootstrap
        # instead of generating random indices for the sampling for each
        # iteration, and for each peptide, we'll generate a batch of random numbers
        # now and pull from them later.
        # the counter will keep track of which portion of the pool we're using
        counter = 0
        rand_pool = np.random.rand(np.sum(obs_per_seq) * k)

        for i in range(0, num_peptides): # for each peptide sequence
          num_obs = obs_per_seq[i]
          for j in range(0, k): # for each iteration:
            # re-estimate mu from the resampled mu_preds
            # TODO: choice also of mean, weighted mean
            mu_k[i][j] = np.median(mu_preds[i][\
              (rand_pool[counter:(counter+num_obs)] * num_obs).astype(int, copy=False)])

            counter = counter + num_obs

      else:
        raise ConfigFileError('Invalid bootstrap method. Please choose \"parametric\" or \"non-parametric.\"')
      
      # map of stan_peptide_id onto 1:num_peptides
      pep_inds = {ind: var for var, ind in enumerate(exp_peptides)}

      # for each bootstrap iteration:
      for j in range(0, k):
        # evaluate the transformed RTs (predicted mus) on distributions
        # with the bootstrapped, estimated mus as the means.
        rt_plus = rt_plus + laplace.pdf(exp['retention_time'], \
          loc=model['ref_to_rt'](\
            exp, mu_k[:,j][exp['stan_peptide_id'].map(pep_inds)], params), \
          scale=exp['sigmaij'])

      # divide total likelihood by # of iterations to normalize to area of 1
      rt_plus = rt_plus / k

    else:
      # not using bootstrap, but using adjusted mu as a point estimate
      # for updating the confidence
      rt_plus = model['rt_plus_func'](exp)

    #                                         P(RT|delta=0)*P(delta=0)
    # PEP.new = P(delta=0|RT) =   ---------------------------------------------------
    #                             P(RT|delta=0)*P(delta=0) + P(RT|delta=1)*P(delta=1)
    #                         
    # delta=1 = Correct ID (true positive)
    # delta=0 = Incorrect (false positive)
    
    # P(RT|delta=0) = probability of peptides RT, given that PSM is incorrect
    #           estimate empirical density of RTs over the experiment
    
    rt_minus = model['rt_minus_func'](exp)

    # P(delta=0) = probability that PSM is incorrect (PEP)
    # P(delta=1) = probability that PSM is correct (1-PEP)
    
    # P(RT|delta=1) = probability that given the correct ID, the RT falls in the
    #           normal distribution of RTs for that peptide, for that experiment
      
    # delta=1 = Correct ID (true positive)
    # delta=0 = Incorrect (false positive)
    # 
    pep_new = (rt_minus * exp['pep']) / \
      ((rt_minus * exp['pep']) + (rt_plus * (1.0 - exp['pep'])))

    # for PSMs for which we have alignment/update data
    exp_new = pd.DataFrame({
        'rt_minus':          rt_minus.tolist(),
        'rt_plus':           rt_plus.tolist(),
        'mu':                exp['mu'].values.tolist(),
        'muij':              exp['muij'].values.tolist(),
        'sigmaij':           exp['sigmaij'].values.tolist(),
        'pep_new':           pep_new.tolist(),

        'id':                exp['id'].values,
        'exp_id':            exp['exp_id'].values,
        'peptide_id':        exp['peptide_id'].values,
        'stan_peptide_id':   exp['stan_peptide_id'].values,
        'input_id':          exp['input_id'].values,
        'exclude':           exp['exclude'].values
    })
    # append to master DataFrame and continue
    df_new = df_new.append(exp_new)

  # reorder by ID and reset the index
  df_new = df_new.sort_values('id')
  df_new = df_new.reset_index(drop=True)

  return df_new

def write_output(df, out_path, config):
  # remove diagnostic columns, unless they are specified to be kept
  if 'add_diagnostic_cols' not in config or config['add_diagnostic_cols'] == False:
    df = df.drop(['remove', 'exclude', 'mu', 'muij', 
      'rt_minus', 'rt_plus', 'sigmaij', 
      'input_id', 'exp_id', 'peptide_id', 'stan_peptide_id'], axis=1)

  # filter by 1% FDR?
  if 'psm_fdr_threshold' in config and type(config['psm_fdr_threshold']) == float:
    if config['psm_fdr_threshold'] <= 0:
      logger.warning('FDR threshold equal to or below 0. Please provide a value between 0 and 1. Ignoring...')
    elif config['psm_fdr_threshold'] >= 1:
      logger.warning('FDR threshold equal to or greater than 1. Please provide a value between 0 and 1. Ignoring...')
    else:
      to_remove = (df['q-value'] > config['psm_fdr_threshold'])
      logger.info('{}/{} ({:.2%}) PSMs removed at a threshold of {:.2%} FDR.'.format(np.sum(to_remove), df.shape[0], np.sum(to_remove) / df.shape[0], config['psm_fdr_threshold']*100))
      df = df[~to_remove].reset_index(drop=True)
    
  df.to_csv(out_path, sep='\t', index=False)

def main():
  start = time.time()

  # load command-line args
  parser = argparse.ArgumentParser()  
  add_global_args(parser)
  args = parser.parse_args()

  # load config file
  # this function also creates the output folder
  config = read_config_file(args)

  # initialize logger
  init_logger(config['verbose'], os.path.join(config['output'], 'dart.log'), config['log_file'])

  logger.info('Converting files and filtering PSMs')
  df, df_original = process_files(config)
  logger.info('Finished converting files and filtering PSMs.')

  # load params, either from a defined folder or from running the alignment
  params = {}
  if type(config['params_folder']) is str:
    params = load_params_from_file(config['params_folder'])
  else:
    logger.info('Beginning alignment procedure')
    params = align(df, config)
    logger.info('Alignment procedure finished')

  # now we have the params, run the update
  logger.info('Updating PEPs with alignment data...')
  df_new = update(df, params, config)

  # save the sparse combined input file?
  #df_new.to_csv(os.path.join(args.output, 'df_converted.txt'), sep='\t', index=False)

  # add new columns to original DF, and remove the duplicate ID column
  logger.info('Concatenating results to original data...')

  df_adjusted = pd.concat([df_original.loc[~df_original['remove']].reset_index(drop=True), df_new.drop(['id', 'input_id'], axis=1).reset_index(drop=True)], axis=1)

  # add rows of PSMs originally removed from analysis
  if np.sum(df_original['remove']) > 0:
    logger.info('Reattaching {} PSMs excluded from initial filters'.format(df_original['remove'].sum()))
    # store a copy of the columns and their order for later
    df_cols = df_adjusted.columns
    # concatenate data frames
    df_adjusted = pd.concat([df_adjusted, df_original.loc[df_original['remove']]], axis=0, ignore_index=True)
    # pd.concat reindexes the order of the columns, 
    # so just order it back to what it used to be
    df_adjusted = df_adjusted.reindex(df_cols, axis=1)

  # sort by ID, and reset index
  df_adjusted = df_adjusted.sort_values(['id'])
  df_adjusted = df_adjusted.reset_index(drop=True)

  # add pep_updated column - which is pep_new, with the NaNs filled in
  # with the old PEPs.
  df_adjusted['pep_updated'] = df_adjusted['pep_new']
  df_adjusted['pep_updated'][pd.isnull(df_adjusted['pep_new'])] = \
    df_adjusted[config['col_names']['pep']][pd.isnull(df_adjusted['pep_new'])]

  # add q-value (FDR) column
  # rank-sorted, cumulative sum of PEPs is expected number of false positives
  # q-value is just that vector divided by # of observations, to get FDR
  df_adjusted['q-value'] = \
    ( \
      np.cumsum(df_adjusted['pep_updated'][np.argsort(df_adjusted['pep_updated'])]) / \
      np.arange(1, df_adjusted.shape[0]+1) \
    )[np.argsort(np.argsort(df_adjusted['pep_updated']))]


  # print figures?
  if config['print_figures']:
    figures(df_adjusted, config, params)

  # overwrite PEP?
  # if true, then store old PEP in "Spectra PEP" column,
  # and put the updated PEP in "PEP" column.
  # then drop the pep_new and pep_updated columns
  if config['overwrite_pep']:
    logger.info('Overwriting PEP column with new PEP. Saving old PEP in \"Spectra PEP\" column.')
    df_adjusted['Spectra PEP'] = df_adjusted[config['col_names']['pep']]
    df_adjusted[config['col_names']['pep']] = df_adjusted['pep_updated']
    df_adjusted = df_adjusted.drop(['pep_new', 'pep_updated'], axis=1)

  # tell the user whether or not to expect diagnostic columns
  if config['add_diagnostic_cols']:
    logger.info('Adding diagnostic columns to output')

  # if no output is specified, throw an error
  if config['save_combined_output'] == False and config['save_separate_output'] == False:
    raise ConfigFileError('No output format specified. Either set \"save_combined_output\" to true, or set \"save_separate_output\" to true.')

  # write to file
  if config['save_combined_output']:
    # if combining input files, then write to one combined file
    out_path = os.path.join(config['output'], config['combined_output_name'])
    logger.info('Combining input file(s) and writing adjusted data file to {} ...'.format(out_path))
    write_output(df_adjusted, out_path, config)
  
  if config['save_separate_output']:
    # if keeping input files separate, then use 'input_id' to retain the
    # order in which the input files were passed in
    logger.info('Saving output to separate files...')
    for i, f in enumerate(config['input']):

      # get output extension
      # default to the same extension as the input
      # if one in the config file exists, use that instead
      out_ext = os.path.splitext(os.path.basename(f))[1]
      if config['output_ext'] is not None:
        out_ext = config['output_ext']

      # contruct output path based on which input file it was
      out_path = os.path.join(config['output'], 
        os.path.splitext(os.path.basename(f))[0] + \
        config['output_suffix'] + '_' + str(i) + out_ext)

      # if saving back to the original input folder,
      # then base the output file name on the input file name instead.
      # no need to number these.
      if config['save_in_input_folder']:
        out_path = os.path.join(os.path.dirname(f),
          os.path.splitext(os.path.basename(f))[0] + \
          config['output_suffix'] + out_ext)

      logger.info('Saving input file {} to {}'.format(i, out_path))
      df_a = df_adjusted.loc[df_adjusted['input_id'] == i]
      # save to file
      # other data formats might have a different separator, or have an index column
      write_output(df_a, out_path, config)

  logger.info('Done! Process took {:.3f} seconds'.format(time.time() - start))

if __name__ == '__main__':
  main()




