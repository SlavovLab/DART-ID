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
from scipy.stats import norm, lognorm

logger = logging.getLogger('root')

def update(dfa, params, config):
  dff = dfa[~(dfa['exclude'])]
  dff = dff.reset_index(drop=True)

  logger.info('{} / {} ({:.2%}) confident, alignable observations (PSMs) after filtering.'.format(dff.shape[0], dfa.shape[0], dff.shape[0] / dfa.shape[0]))

  # refactorize peptide id into stan_peptide_id, 
  # to preserve continuity when feeding data into STAN
  dff['stan_peptide_id'] = dff['sequence'].map({
    ind: val for val, ind in enumerate(dff['sequence'].unique())})

  num_experiments = dff['exp_id'].max() + 1
  num_observations = dff.shape[0]
  num_peptides = dff['peptide_id'].max() + 1
  exp_names = np.sort(dff['raw_file'].unique())
  mean_log_rt = np.mean(np.log(dff['retention_time']))
  sd_log_rt = np.std(np.log(dff['retention_time']))
  max_rt = dff['retention_time'].max()
  pep_id_list = dff['peptide_id'].unique()

  model = get_model_from_config(config)

  # output table
  df_new = pd.DataFrame()

  for i, e in enumerate(np.sort(dff['exp_id'].unique())):
    exp_name = exp_names[i]
    logger.info('Updating PEPs for experiment #{} - {}'.format(i+1, exp_name))
    
    exp = dfa[dfa['exp_id'] == e]
    exp = exp.reset_index(drop=True)
    
    # not all peptides in this experiment have data from the model
    # we can only update those that have that data. others will not be touched
    exp_matches = np.isin(exp['peptide_id'].values, pep_id_list)
    exp_f = exp[exp_matches]
    exp_f = exp_f.reset_index(drop=True)

    # ensure that PEP does not exceed 1
    # will result in incorrect negative densities when applying mixture model
    exp_f['pep'][exp_f['pep'] > 1] = 1

    # convert peptide_id to stan_peptide_id
    exp_f['stan_peptide_id'] = exp_f['peptide_id'].map({
      ind: val for val, ind in enumerate(pep_id_list)
    })
    exp_peptides = exp_f['stan_peptide_id'].unique()
    exp_f['mu'] = params['peptide']['mu'].values[exp_f['stan_peptide_id']]

    # get muij from mu, using the monotone transformation parameters from this experiment
    exp_f['muij'] = models[model]['muij_func'](exp_f, i, params)
    # get sigmaij from data and parameters
    exp_f['sigmaij'] = models[model]['sigmaij_func'](exp_f, i, params)

    #                                    P(RT|PSM-)*P(PSM-)
    # PEP.new = P(PSM-|RT) =  ---------------------------------------
    #                         P(RT|PSM-)*P(PSM-) + P(RT|PSM+)*P(PSM+)
    #                         
    # PSM+ <- PSM is Correct
    # PSM- <- PSM is Incorrect
    
    # P(RT|PSM-) = probability of peptides RT, given that PSM is incorrect
    #           estimate empirical density of RTs over the experiment
    
    rt_minus = models[model]['rt_minus_func'](exp_f)

    # P(PSM-) = probability that PSM is incorrect (PEP)
    # P(PSM+) = probability that PSM is correct (1-PEP)
    
    # P(RT|PSM+) = probability that given the correct ID, the RT falls in the
    #           normal distribution of RTs for that peptide, for that experiment
    #
    # this is defined in fit_RT3d.stan as a mixture between 2 normal distributions
    # where one distribution for the peptide RT is weighted by 1-PEP
    # and the other distribution for all RTs is weighted by PEP
    # -- summing to a total density of 1

    rt_plus = models[model]['rt_plus_func'](exp_f)

    # now we can update the PEP
    #                                    P(RT|PSM-)*P(PSM-)
    # PEP.new = P(PSM-|RT) =  ---------------------------------------
    #                         P(RT|PSM-)*P(PSM-) + P(RT|PSM+)*P(PSM+)
    # + | PSM = Correct
    # - | PSM = Incorrect
    pep_new = (rt_minus * exp_f['pep']) / \
      ((rt_minus * exp_f['pep']) + (rt_plus * (1 - exp_f['pep'])))

    # for PSMs for which we have alignment/update data
    exp_new = pd.DataFrame({
        'rt_minus':          rt_minus.tolist(),
        'rt_plus':           rt_plus.tolist(),
        'mu':                exp_f['mu'].values.tolist(),
        'muij':              exp_f['muij'].values.tolist(),
        'sigmaij':           exp_f['sigmaij'].values.tolist(),
        'pep_new':           pep_new.tolist(),

        'id':                exp_f['id'],
        'exp_id':            exp_f['exp_id'],
        'peptide_id':        exp_f['peptide_id'],
        'stan_peptide_id':   exp_f['stan_peptide_id'],
        'input_id':          exp_f['input_id'],
        'exclude':           exp_f['exclude']
    })
    # for PSMs without alignment/update data
    exp_new = exp_new.append(pd.DataFrame({
        'rt_minus':          np.nan,
        'rt_plus':           np.nan,
        'mu':                np.nan,
        'muij':              np.nan,
        'sigmaij':           np.nan,
        'pep_new':           np.nan,

        'id':                exp['id'][~(exp_matches)],
        'exp_id':            exp['exp_id'][~(exp_matches)],
        'peptide_id':        exp['peptide_id'][~(exp_matches)],
        'stan_peptide_id':   np.nan,
        'input_id':          exp['input_id'][~(exp_matches)],
        'exclude':           exp['exclude'][~(exp_matches)]
    }))
    # append to master DataFrame and continue
    df_new = df_new.append(exp_new)

  # reorder by ID and reset the index
  df_new = df_new.sort_values('id')
  df_new = df_new.reset_index(drop=True)

  return df_new

def write_output(df, out_path, config):
  # remove diagnostic columns, unless they are specified to be kept
  if not config['add_diagnostic_cols']:
    df = df.drop(['input_exclude', 'exclude', 'mu', 'muij', 
      'rt_minus', 'rt_plus', 'sigmaij', 
      'input_id', 'exp_id', 'peptide_id', 'stan_peptide_id'], axis=1)
    
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

  df_adjusted = pd.concat([df_original.loc[~df_original['input_exclude']].reset_index(drop=True), df_new.drop(['id', 'input_id'], axis=1).reset_index(drop=True)], axis=1)

  # add rows of removed experiments (done with the --remove-exps options)
  if config['exclude'] is not None or config['include'] is not None:
    logger.info('Reattaching {} PSMs excluded with the \
      experiment blacklist'.format(df_original['input_exclude'].sum()))
    # store a copy of the columns and their order for later
    df_cols = df_adjusted.columns
    # concatenate data frames
    df_adjusted = pd.concat([df_adjusted, df_original.loc[df_original['input_exclude']]], axis=0, ignore_index=True)
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




