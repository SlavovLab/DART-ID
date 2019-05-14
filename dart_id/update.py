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
from dart_id.exceptions import ConfigFileError
from dart_id.fido.BayesianNetwork import run_internal
from dart_id.figures import figures
from dart_id.models import models, get_model_from_config
from dart_id.helper import add_global_args, read_config_file, init_logger, load_params_from_file
from scipy.stats import norm, lognorm, laplace, bernoulli, uniform

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
  # scaled sigma is the same ratio of muij / mu applied to sigmaij
  dfa['sigma_pred'] = dfa['sigmaij'] * dfa['mu_pred'] / dfa['muij']

  # get parameters for the null distributions for each experiment
  null_dists = dfa.groupby('exp_id')['retention_time'].agg([np.mean, np.std])
  #null_dists = np.array([norm(loc=null_dists.loc[i, 'mean'], scale=null_dists.loc[i, 'std']) for i in range(0, num_experiments)])
  # first column is mean, second is std
  null_dists = np.array([null_dists['mean'].values, null_dists['std'].values]).T

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

  logger.info('Updating PEPs...')
  for i, e in enumerate(np.sort(dfa['exp_id'].unique())):

    exp_name = exp_names[i]

    exp = dfa[dfa['exp_id'] == e]
    exp = exp.reset_index(drop=True)

    exp_peptides = exp['stan_peptide_id'].unique()

    logger.info('Exp ({} / {}) - {} - ({} Peptides, {} PSMs)'.format(i+1, num_experiments, exp_name, len(exp_peptides), exp.shape[0]))

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
        ['stan_peptide_id', 'pep', 'mu_pred', 'mu', 'sigma_pred', 'exp_id']]

      # extract relevant values for each peptide
      mu_preds    = dfe.groupby('stan_peptide_id')['mu_pred'].apply(lambda x: x.values).values.tolist()
      mus         = dfe.groupby('stan_peptide_id')['mu'].apply(lambda x: x.values).values.tolist()
      sigma_preds = dfe.groupby('stan_peptide_id')['sigma_pred'].apply(lambda x: x.values).values.tolist()
      peps        = dfe.groupby('stan_peptide_id')['pep'].apply(lambda x: x.values).values.tolist()
      exp_ids     = dfe.groupby('stan_peptide_id')['exp_id'].apply(lambda x: x.values).values.tolist()
      
      # number of observations per peptide sequence
      obs_per_seq = [len(peptide) for peptide in mu_preds]
      num_peptides = len(mu_preds)

      # the number of observations per peptide -- used in loop
      num_obs = 0
      # matrix of n by k estimated mus from the bootstrapping
      # will iterate over in the loop after the immediate one
      mu_k = np.zeros((num_peptides, k))

      if bootstrap_method == 'parametric' or bootstrap_method == 'parametric_mixture':

        t_laplace_samples = 0
        t_coin_flips = 0
        t_null_samples = 0
        t_loop_indexing = 0
        t_medians = 0

        # create pool of coin flips, instead of sampling for every peptide
        # the pool is uniformly distributed from 0 to 1, and "successful" coin flip
        # is determined by whether or not the sample from the pool is less than the
        # measured PEP
        _time = time.time()
        coin_flip_pool = 0
        if bootstrap_method == 'parametric_mixture':
          coin_flip_pool = uniform.rvs(size=(np.sum(obs_per_seq) * k))

        t_coin_flips += (time.time() - _time)
        coin_counter = 0

        # create a pool of laplace samples, to pull from for each peptide
        _time = time.time()
        sample_pool = laplace.rvs(size=(np.sum(obs_per_seq) * k))
        t_laplace_samples += (time.time() - _time)
        # keep track of where we are in the pool with a counter
        sample_counter = 0

        # parametric bootstrap
        for i in range(0, num_peptides):
          num_obs = obs_per_seq[i]
          
          _time = time.time()
          # sample num_obs synthetic RTs for k bootstrap iterations
          # do the sampling in a big pool, then shape to matrix where
          # rows correspond to bootstrap iters and columns correspond to sample observations
          #samples = laplace.rvs(size=(k * num_obs)).reshape(k, num_obs)
          
          # draw samples from sample pool, reshape into matrix
          _time = time.time()
          samples = sample_pool[sample_counter:(sample_counter+(k*num_obs))]
          samples = samples.reshape(k, num_obs)

          # increment sample counter
          sample_counter += (k*num_obs)

          #mu_med = np.median()
          # shift and scale sampled RTs by mu and sigma_pred, respectively
          samples = (samples * sigma_preds[i]) + mu_preds[i]
          #samples = (samples * sigma_preds[i]) + mu_med
          t_laplace_samples += (time.time() - _time)

          if bootstrap_method == 'parametric_mixture':
            # sample from mixture distribution
            _time = time.time()
            # actually faster to just replicate the sample matrix and then
            # take subindices from that instead of sampling from null every
            # iteration of the loop below. this seems inefficient, especially
            # if given very small PEPs, but still better than sampling every iteration.
            # could probably optimize the size of the null sample matrix by
            # looking at predicted false positive rates, but for now we're
            # just going with worst case scenario and assuming for all false positives.
            null_samples = norm.rvs(size=(k * num_obs)).reshape(k, num_obs)
            # shift and scale sampled RTs by mean and std of null dists
            null_samples = (null_samples * null_dists[exp_ids[i],1]) + null_dists[exp_ids[i],0]
            t_null_samples += (time.time() - _time)

            _time = time.time()
            for j in range(0, num_obs): # for each observation in the matrix
              # take a chunk of the coin flip pool
              fp = (coin_flip_pool[coin_counter:(coin_counter + k)] < peps[i][j]).astype(bool)
              coin_counter += k
              # overwrite original samples with samples from null distribution
              samples[fp, j] = null_samples[fp, j]
            t_loop_indexing += (time.time() - _time)

          _time = time.time()
          # now take the median of each row and store it in mu_k
          mu_k[i] = np.median(samples, axis=1)
          # or take the weighted mean
          #weights = ((1 - peps[i]) - (1 - config['pep_threshold'])) / config['pep_threshold']
          #mu_k[i] = (np.sum(samples * weights, axis=1) / np.sum(weights))
          t_medians += (time.time() - _time)

        logger.debug('laplace sampling: {:.1f} ms'.format(t_laplace_samples*1000))
        logger.debug('coin flips: {:.1f} ms'.format(t_coin_flips*1000))
        logger.debug('null sampling: {:.1f} ms'.format(t_null_samples*1000))
        logger.debug('loop indexing: {:.1f} ms'.format(t_loop_indexing*1000))
        logger.debug('taking medians: {:.1f} ms'.format(t_medians*1000))

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
      
      _t_dist_building = time.time()
      # map of stan_peptide_id onto 1:num_peptides
      pep_inds = {ind: var for var, ind in enumerate(exp_peptides)}
      pep_inds = exp['stan_peptide_id'].map(pep_inds)

      # for each bootstrap iteration:
      for j in range(0, k):
        # evaluate the transformed RTs (predicted mus) on distributions
        # with the bootstrapped, estimated mus as the means.
        #rt_plus = rt_plus + laplace.pdf(exp['retention_time'], \
        #  loc=model['ref_to_rt'](exp, mu_k[:,j][pep_inds], params), \
        #  scale=exp['sigmaij'])

        rt_plus = rt_plus + laplace.pdf(exp['mu_pred'], \
          loc=mu_k[:,j][pep_inds], \
          scale=exp['sigma_pred'])

      # divide total likelihood by # of iterations to normalize to area of 1
      rt_plus = rt_plus / k

      logger.debug('distribution building: {:.1f} ms'.format((time.time() - _t_dist_building)*1000))

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
    df = df.drop(['pep_new', 'participated', 'exclude', 'mu', 'muij', 
      'rt_minus', 'rt_plus', 'sigmaij', 'residual',
      'input_id', 'exp_id', 'peptide_id', 'stan_peptide_id'], axis=1)

  # filter by PSM FDR?
  if 'psm_fdr_threshold' in config and type(config['psm_fdr_threshold']) == float:
    if config['psm_fdr_threshold'] <= 0:
      logger.warning('PSM FDR threshold equal to or below 0. Please provide a value between 0 and 1. Ignoring...')
    elif config['psm_fdr_threshold'] >= 1:
      logger.warning('PSM FDR threshold equal to or greater than 1. Please provide a value between 0 and 1. Ignoring...')
    else:
      to_remove = (df['dart_qval'] > config['psm_fdr_threshold'])
      logger.info('{}/{} ({:.2%}) PSMs removed at a threshold of {:.2%} FDR.'.format(np.sum(to_remove), df.shape[0], np.sum(to_remove) / df.shape[0], config['psm_fdr_threshold']))
      df = df[~to_remove].reset_index(drop=True)

  # filter by protein FDR?
  if 'protein_fdr_threshold' in config and type(config['protein_fdr_threshold']) == float:
    if 'razor_protein_fdr' in df.columns:
      if config['protein_fdr_threshold'] <= 0:
        logger.warning('Protein FDR threshold equal to or below 0. Please provide a value between 0 and 1. Ignoring...')
      elif config['protein_fdr_threshold'] >= 1:
        logger.warning('Protein FDR threshold equal to or greater than 1. Please provide a value between 0 and 1. Ignoring...')
      else:
        to_remove = ((df['razor_protein_fdr'] > config['protein_fdr_threshold']) | pd.isnull(df['razor_protein_fdr']))
        logger.info('{}/{} ({:.2%}) PSMs removed at a threshold of {:.2%} Protein FDR.'.format(np.sum(to_remove), df.shape[0], np.sum(to_remove) / df.shape[0], config['protein_fdr_threshold']))
        df = df[~to_remove].reset_index(drop=True)
    else:
      raise ConfigFileError('Protein FDR threshold specified, but no protein inference run with this analysis. Please set \"run_pi\" to true and fill out all the respective parameters.')
    
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
  if 'params_folder' in config and type(config['params_folder']) is str:
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
  df_adjusted = pd.concat([ \
    df_original.loc[~df_original['remove']].reset_index(drop=True), \
    df_new.drop(['id', 'input_id'], axis=1).reset_index(drop=True)], axis=1)

  # add rows of PSMs originally removed from analysis
  if np.sum(df_original['remove']) > 0:
    logger.info('Reattaching {} PSMs excluded from initial filters'.format(df_original['remove'].sum()))
    # store a copy of the columns and their order for later
    df_cols = df_adjusted.columns
    # concatenate data frames
    df_adjusted = pd.concat([ \
      df_adjusted, \
      df_original.loc[df_original['remove']]], 
      axis=0, ignore_index=True, sort=True)
    # pd.concat reindexes the order of the columns, 
    # so just order it back to what it used to be
    df_adjusted = df_adjusted.reindex(df_cols, axis=1)

  # sort by ID, and reset index
  df_adjusted = df_adjusted.sort_values(['id'])
  df_adjusted = df_adjusted.reset_index(drop=True)

  # add residual RT (alignment error) column
  df_adjusted['residual'] = np.abs(\
    df_adjusted[config['col_names']['retention_time']] - df_adjusted['muij'])

  # add dart_PEP column - which is pep_new, with the NaNs filled in
  # with the old PEPs.
  df_adjusted['dart_PEP'] = df_adjusted['pep_new']
  df_adjusted['dart_PEP'][pd.isnull(df_adjusted['pep_new'])] = \
    df_adjusted[config['col_names']['pep']][pd.isnull(df_adjusted['pep_new'])]
  # make sure that updated PEP does not exceed 1
  df_adjusted['dart_PEP'][df_adjusted['dart_PEP'] > 1] = 1

  # add q-value (FDR) column
  # rank-sorted, cumulative sum of PEPs is expected number of false positives
  # q-value is just that vector divided by # of observations, to get FDR
  logger.info('Calculating FDR (q-values)')
  
  # q-value, without fixing # of false positives to a discrete number
  #df_adjusted['q-value'] = \
  #  ( \
  #    np.cumsum(df_adjusted['dart_PEP'][np.argsort(df_adjusted['dart_PEP'])]) / \
  #    np.arange(1, df_adjusted.shape[0]+1) \
  #  )[np.argsort(np.argsort(df_adjusted['dart_PEP']))]
  
  # q-value, by fixing # of false positives to a discrete number
  # for now, set all null PEPs to 1. we'll remember the index and set them back to nan later
  null_peps = pd.isnull(df_adjusted['dart_PEP'])
  if null_peps.sum() > 0:
    df_adjusted['dart_PEP'][null_peps] = 1

  # get the index order of sorted PEPs
  pep_order = np.argsort(df_adjusted['dart_PEP'])
  # Take the ceiling of the cumulative sum of the sorted PEPs to get the pessimistic
  # estimate of the number of false positives when selecting at that level.
  # because using ceiling, PSMs with different PEPs but within the same relative interval
  # will get the same "num_fp" value.
  num_fp = np.ceil(np.cumsum(df_adjusted['dart_PEP'][pep_order])).astype(int)
  # count the number of occurrences of num_fp and sum them up to get the sample size for each
  # discrete false positive # threshold
  fp_counts = np.cumsum(num_fp.value_counts().sort_index()).values
  # divide # of false positivies by sample size to get q-value. sorting the index order brings
  # the order of values back to their original form
  df_adjusted['dart_qval'] = (num_fp / fp_counts[num_fp-1]).values[np.argsort(pep_order.values)]

  # set null PEPs and q-values back to nan
  if null_peps.sum() > 0:
    df_adjusted['dart_PEP'][null_peps] = np.nan
    df_adjusted['dart_qval'][null_peps] = np.nan

  # rename 'remove' column - which indicates whether or not the PSM participated in the
  # DART-ID alignment and update
  df_adjusted['participated'] = ~df_adjusted['remove']

  ## Run protein inference (fido)?
  if 'run_pi' in config and config['run_pi']:
    logger.info('Running protein inference with Fido...')

    # build fido options into a dict (parameter_map)
    parameter_map = {
      'gamma': config['pi_gamma'] if 'pi_gamma' in config else None,
      'alpha': config['pi_alpha'] if 'pi_alpha' in config else None,
      'beta':  config['pi_beta']  if 'pi_beta'  in config else None,

      'connected_protein_threshold':  config['pi_connected_protein_thresh'],
      'omit_clean_peptide_name':     ~config['pi_clean_peptide_name'],
      'all_psms':                     config['pi_use_all_psms'],
      'group_proteins':               config['pi_group_proteins'],
      'prune_low_scores':             config['pi_prune_low_scores'],
      'parameter_accuracy':           config['pi_parameter_accuracy'],

      'proteins_column':         config['col_names']['proteins'],
      'protein_delimiter':       config['pi_protein_delimiter'],
      'leading_protein_column':  config['col_names']['leading_protein'],
      'decoy_tag':               config['pi_decoy_tag'],

      'sequence_column':         config['col_names']['sequence'],
      #'error_prob_column':       config['col_names']['pep']
      'error_prob_column':       'dart_PEP',

      # pass in output folder so fido can save some intermediate and output files
      'output': config['output']
    }
    logger.debug('parameter_map for fido:')
    logger.debug(str(parameter_map))

    # run fido subroutine
    df_adjusted = run_internal(df_adjusted, parameter_map)

    logger.info('Fido finished')
    logger.info('FDR for PSM\'s razor protein, from protein inference, placed in \"razor_protein_fdr\" column')

  # print figures?
  if config['print_figures']:
    figures(df_adjusted, config, params)

  # overwrite PEP?
  # if true, then store old PEP in "Spectra PEP" column,
  # and put the dart PEP in "PEP" column.
  # then drop the pep_new and dart_PEP columns
  if config['overwrite_pep']:
    logger.info('Overwriting PEP column with new PEP. Saving old PEP in \"Spectra PEP\" column.')
    df_adjusted['Spectra PEP'] = df_adjusted[config['col_names']['pep']]
    df_adjusted[config['col_names']['pep']] = df_adjusted['dart_PEP']
    df_adjusted = df_adjusted.drop(['pep_new', 'dart_PEP'], axis=1)

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




