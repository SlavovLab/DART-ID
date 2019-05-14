#!/usr/bin/env python3
# coding: utf-8

# converter.py - converters from search-engine inputs to compatible outputs
# for alignment and update

import argparse
import logging
import numpy as np
import os
import pandas as pd
import re

from functools import reduce
from dart_id.exceptions import ConfigFileError, FilteringError
from dart_id.helper import add_global_args, read_config_file, init_logger

logger = logging.getLogger('root')

# all filter funcs take in the df, config object, and the filter object
# as inputs, and output the filter vector (True/False), where True means
# to filter out that particular row

def filter_exclude_filename(df, config, _filter):
  # remove experiments from blacklist
  if 'expr' not in _filter or len(_filter['expr']) == 0:
    raise ConfigFileError('Experiment exclusion by raw file name provided, but exclusion expression is defined incorrectly. Please provide a regular expression, with the \"expr\" key, to match against raw file names, that determines whether or not that raw file is included in the alignment procedure.')
  
  # see if any raw file names match the user-provided expression
  exclude_exps = list(filter(lambda x: re.search(r'' + _filter['expr'] + '', x), 
    df['raw_file'].unique()))
  # remove excluded rows from the sparse dataframe,
  # but keep them in the original data frame, so that we can stitch together
  # the final output later
  exclude_exps = df['raw_file'].isin(exclude_exps).values
  df = df[~exclude_exps]
  logger.info('Filtering out {} observations matching \"{}\"'.format(np.sum(exclude_exps), _filter['expr']))

  return exclude_exps

def filter_include_filename(df, config, _filter):
  # exclusively include experiments from whitelist
  if 'expr' not in _filter or len(_filter['expr']) == 0:
    raise ConfigFileError('Experiment inclusion by raw file name provided, but inclusion expression is missing or defined incorrectly. Please provide a regular expression, with the \"expr\" key, to match against raw file names, that determines whether or not that raw file is included in the alignment procedure.')
  
  # get matches for this expression
  include_exps = list(filter(lambda x: re.search(r'' + _filter['expr'] + '', x), df['raw_file'].unique()))
  # only keep rows that are in these raw file matches
  include_exps = df['raw_file'].isin(include_exps).values
  df = df[include_exps]
  logger.info('Keeping {} observations out of {} matching inclusion expression \"{}\"'.format(np.sum(include_exps), 
      df.shape[0], _filter['expr']))

  # filter out the opposite of the included experiments
  return ~include_exps

def filter_uniprot_exclusion_list(df, config, _filter):
  """
  Filter proteins from exclusion list using UniProt IDs
  """

  exclusion_list = []

  # parse exclusion list
  # if exclusion_list param is a path, then load the IDs from that path
  if 'file' in _filter and _filter['file'] is not None:
    # load UniProt IDs from file line-by-line
    # first expand user or any vars
    _filter['file'] = os.path.expanduser(_filter['file'])
    _filter['file'] = os.path.expandvars(_filter['file'])

    # open the exclusion list file and read in the UniProt IDs, line by line
    try:
      with open(_filter['file'], 'r') as f:
        logger.info('Loading UniProt IDs from exclusion list file {} ...'.format(_filter['file']))
        exclusion_list = [line.rstrip('\n') for line in f]
        logger.info('Loaded {} proteins from exclusion list.'.format(len(exclusion_list)))
    except EnvironmentError:
      raise ConfigFileError('Exclusion list file {} not found. Please provide a path to a file with UniProt IDs separated by line'.format(_filter['file']))

  elif 'list' in _filter and len(_filter['list']) > 0:
    # load UniProt IDs from the configuration file
    exclusion_list = _filter['list']
    logger.info('Loading {} UniProt IDs from exclusion list as defined in config file'.format(len(exclusion_list)))
  else:
    raise ConfigFileError('No exclusion list file or list of UniProt IDs provided. Please provide a path to a file with UniProt IDs separated by line with the \"file\" key, or provide a python list of UniProt IDs with the \"list\" key. If not using a UniProt ID exclusion list, then comment out the \"uniprot_exclusion\" key from the filter list.')

  # filter exclusion list
  if len(exclusion_list) > 0:
    logger.info('UniProt IDs from exclusion list: {}'.format(exclusion_list))

    # we could only match the excluded IDs to the razor protein,
    # but we can be more strict and match the blacklisted IDs to the entire protein
    # string, containing all possible proteins
    pat = reduce((lambda x, y: x + '|' + y), exclusion_list)
    blacklist_filter = df['proteins'].str.contains(pat)
    blacklist_filter[pd.isnull(blacklist_filter)] = False

    logger.info('Filtering out {} PSMs from the exclusion list'.format(np.sum(blacklist_filter)))
    return blacklist_filter
  else:
    raise ConfigFileError('Exclusion list found and loaded, but no UniProt IDs found. Check the format of the file, or the list in the config file.')

def filter_contaminant(df, config, _filter):
  """
  Filter contaminants, as marked by the search engine
  Looking for a contaminant tag in the leading_protein column
  """

  # load the tag in from the config file
  CON_TAG = _filter['tag']

  if CON_TAG is None:
    raise ConfigFileError('Contaminant filter has no tag specified. Please provide an expression with the \"tag\" key that matches a contaminant tag from the search engine output. If you don\'t want to apply this filter, then comment out the \"contaminant\" key from the filter list.')
  elif type(CON_TAG) is not str:
    raise ConfigFileError('Contaminant tag {} is not of type string. Please provide an expression with the \"tag\" key that matches a contaminant tag from the search engine output.'.format(CON_TAG))
  elif len(CON_TAG) <= 0:
    raise ConfigFileError('Contaminant tag empty. Please provide an expression with the \"tag\" key that matches a contaminant tag from the search engine output.')

  # search for the tag in the 'proteins' column
  filter_con = df['proteins'].str.contains(CON_TAG)
  filter_con[pd.isnull(filter_con)] = False

  logger.info('Filtering out {} PSMs as contaminants with tag \"{}\"'.format(np.sum(filter_con), CON_TAG))
  return filter_con

def filter_decoy(df, config, _filter):
  """
  Filter decoys, as marked by the search engine
  Looking for a decoy tag in the leading_protein column
  """

  # load the tag in from the config file
  REV_TAG = _filter['tag']

  if REV_TAG is None:
    raise ConfigFileError('Decoy filter has no tag specified. Please provide an expression with the \"tag\" key that matches a decoy tag from the search engine output. If you don\'t want to apply this filter, then comment out the \"decoy\" key from the filter list.')
  elif type(REV_TAG) is not str:
    raise ConfigFileError('Decoy tag {} is not of type string. Please provide an expression with the \"tag\" key that matches a decoy tag from the search engine output.'.format(REV_TAG))
  elif len(REV_TAG) <= 0:
    raise ConfigFileError('Decoy tag empty. Please provide an expression with the \"tag\" key that matches a decoy tag from the search engine output.')

  # search for the tag in the 'leading protein' column
  filter_rev = df['leading_protein'].str.contains(REV_TAG)
  filter_rev[pd.isnull(filter_rev)] = False

  logger.info('Filtering out {} PSMs as decoys with tag \"{}\"'.format(np.sum(filter_rev), REV_TAG))
  return filter_rev

def filter_retention_length(df, config, _filter):
  """
  Filter by retention length, which is a measure of the peak width
  during chromatography.
  """

  # input checks
  if 'value' not in _filter or _filter['value'] is None:
    raise ConfigFileError('No value provided to the retention_length filter. If you don\'t want to apply this filter, then comment out the \"retention_length\" key from the filter list. Or, please provide a threshold with the \"value\" key, and specify whether or not this is an absolute time in minutes with \"dynamic\"=false or \"dynamic\"=true if the filter is a fraction of the experiment run-time. If \"dynamic\"=true, then \"value\" must be between 0 and 1.')
  if 'dynamic' not in _filter or type(_filter['dynamic']) is not bool:
    raise ConfigFileError('Incorrect value provided to the \"dynamic\" field of the retention_length filter. Please provide a bool, either true or false. If \"dynamic\"=false, then \"value\" is the threshold in minutes. If \"dynamic\"=true, then \"value\" must be between 0 and 1.')
  
  if _filter['dynamic']:
    # use the dynamic filter, where the value is a proportion
    # of the max RT (the run-time) of that raw file
    
    # only allow values between 0 and 1 for the dynamic filter
    if _filter['value'] > 1 or _filter['value'] <= 0:
      raise ConfigFileError('Dynamic retention_length filter {} is above 1 or below 0. Please provide a number between 0 and 1, which is the fraction of the max RT for each experiment. e.g., 0.01 means that 1%% of the max RT will be used as the retention_length threshold.'.format(_filter['value']))

    logger.info('Using dynamic retention length of {} * run-time (max RT) for each experiment'.format(_filter['value']))

    # get the max RT for each raw file, reindex to the same dimension as the
    # retention_length column, and then multiply by the filter value
    max_rts = df.groupby('raw_file')['retention_time'].max().values
    
    filter_rtl = max_rts[df['raw_file'].map({ 
      ind: val for val, ind in enumerate(np.sort(df['raw_file'].unique()))
    })] * _filter['value']

    filter_rtl = (df['retention_length'] > filter_rtl)

  else:
    # use a constant filter for the retention length
    logger.info('Using constant retention length (in RT) of {} for all raw files.'.format(_filter['value']))

    # only allow values between 0 and max(RT)
    if _filter['value'] <= 0 or _filter['value'] > np.max(df['retention_time']):
      raise ConfigFileError('retention_length filter {} is not defined or incorrectly defined. Please provide a decimal number between 0.0 and max(RT).'.format(_filter['value']))
    
    filter_rtl = (df['retention_length'] > _filter['value'])

  if _filter['dynamic']:
    logger.info('Filtering out {} PSMs with retention length greater than {:.4f} * max(exp_RT) of each raw file.'.format(np.sum(filter_rtl), _filter['value']))
  else:
    logger.info('Filtering out {} PSMs with retention length greater than {:.4f}'.format(np.sum(filter_rtl), _filter['value']))

  return filter_rtl

def filter_smears(df, config, _filter):
  """
  Filter out "smears". even confidently identified PSMs can have bad chromatography,
  and in that case it is unproductive to include them into the alignment.
  In theory, this should be made redundant by the retention length filter, but
  some PSMs still slip through the cracks of that, possibly because the search engine
  cannot adequately track the elution peak?
  """

  # input checking
  if 'value' not in _filter or _filter['value'] is None:
    raise ConfigFileError('No value provided to the \"smears\" filter. If you don\'t want to apply this filter, then comment out the \"smears\" key from the filter list. Or, please provide a threshold with the \"value\" key, and specify whether or not this is an absolute time in minutes with \"dynamic\"=false or \"dynamic\"=true if the filter is a fraction of the experiment run-time. If \"dynamic\"=true, then \"value\" must be between 0 and 1.')
  if 'dynamic' not in _filter or type(_filter['dynamic']) is not bool:
    raise ConfigFileError('Incorrect value provided to the \"dynamic\" field of the \"smears\" filter. Please provide a bool, either true or false. If \"dynamic\"=false, then \"value\" is the threshold in minutes. If \"dynamic\"=true, then \"value\" must be between 0 and 1.')

  logger.info('Determining RT spread of peptides within each experiment...')
  # for each experiment-peptide pair, get the range of retention times
  # this is the step that could take a long time
  # TODO: optimize this?
  smears = df.groupby(['raw_file', 'sequence'])['retention_time'].apply(np.ptp)
  
  if _filter['dynamic']:
    # use the dynamic filter, where the value is a proportion
    # of the max RT (the run-time) of that raw file
    
    if _filter['value'] > 1 or _filter['value'] <= 0:
      raise ConfigFileError('Dynamic smear filter {} is above 1 or below 0. Please provide a number between 0 and 1, which is the fraction of the max RT for each experiment. e.g., 0.01 means that 1%% of the max RT will be used as the retention_length threshold.'.format(_filter['value']))

    logger.info('Using dynamic smear length (in RT) of {:.4f} * run-time (max RT) for each experiment'.format(_filter['value']))

    max_rts = df.groupby('raw_file')['retention_time'].max().values

    smear_pair_inds = smears.index.to_frame()['raw_file'].values
    smear_pair_inds = pd.Series(smear_pair_inds).map({
      ind: val for val, ind in enumerate(np.sort(df['raw_file'].unique())) })

    # get the (raw_file, sequence) tuples for PSMs with a range above the threshold
    smears = smears[smears > (max_rts[smear_pair_inds] * _filter['value'])].index.values

  else:
    # use a constant filter for the retention length
    logger.info('Using constant smear length (in RT) of {:.4f} for all raw files.'.format(_filter['value']))

    if _filter['value'] <= 0:
      raise ConfigFileError('Smear filter {:.4f} is not defined or incorrectly defined. Please provide a decimal number between 0.0 and max(RT).'.format(_filter['value']))

    # get the (exp_id, peptide_id) tuples for PSMs with a range above the threshold
    smears = smears[smears > _filter['value']].index.values
  
  # map the tuples back to the original data frame, and set smears to be excluded
  smears = pd.Series(list(zip(df['raw_file'], df['sequence']))).isin(smears)

  if _filter['dynamic']:
    logger.info('Filtering out {} PSMs with an intra-experiment RT spread greater than {:.4f} * max(exp_RT) for each raw file.'.format(smears.sum(), _filter['value']))
  else:
    logger.info('Filtering out {} PSMs with an intra-experiment RT spread greater than {:.4f}'.format(smears.sum(), _filter['value']))

  return smears.values

# dictionary of all filter functions
filter_funcs = {
  'exclude_filename':   filter_exclude_filename,
  'include_filename':   filter_include_filename,
  'uniprot_exclusion':  filter_uniprot_exclusion_list,
  'contaminant':        filter_contaminant,
  'decoy':              filter_decoy,
  'retention_length':   filter_retention_length,
  'smears':             filter_smears
}

# columns required for each filter to run
# will skip the filter if this column does not exist in the input file
required_cols = {
  'exclude_filename':   [],
  'include_filename':   [],
  'uniprot_exclusion':  ['proteins'],
  'contaminant':        ['proteins'],
  'decoy':              ['leading_protein'],
  'retention_length':   ['retention_length'],
  'smears':             []
}

def convert(df, config):

  cols = []
  col_names = []
  
  # loop thru all columns listed in the config file
  for col in list(config['col_names'].keys()):
    if config['col_names'][col] is None: 
      logger.debug('Column \"{}\" is left empty in the config file. Skipping...'.format(col))
      continue

    # check if the column specified in the config file exists in the df or not
    if config['col_names'][col] not in df.columns:
      # this is probably grounds to kill the program
      raise ConfigFileError('Column {} of value {} not found in the input file. Please check that this column exists, or leave the field for {} empty in the config file.'.format(col, config['col_names'][col], col))

    # keep track of the column and the column name
    cols.append(config['col_names'][col])
    col_names.append(col)

  # take the subset of the input file, and also rename the columns
  dfa = df[cols]
  dfa.columns = col_names

  return dfa


def filter_psms(df, config):
  logger.info('Filtering PSMs...')

  # load the filtering functions specified by the input config
  # types of filter functions depends on what stage of filtering this is:
  # removing observations or merely excluding from alignment
  filters = config['filters']

  # each filter has a specified required column from the dataframe
  # make sure these columns exist before proceeding
  for i, f in enumerate(filters):
    # make sure filter has required field of 'name'
    if 'name' not in f or f['name'] is None or len(f['name']) == 0:
      raise ConfigFileError('Filter {} does not have the required field name or it is empty. Please provide a name (type) for this filter. Available filters are listed in the example config file, under example/config.yaml'.format(str(f)))

    # for each required column in the filter, check if it exists
    for j in required_cols[f['name']]:
      if j not in df.columns:
        raise ConfigFileError('Filter {} required a data column {}, but this was not found in the input dataframe.'.format(f['name'], j))

  # by default, filter out nothing. we'll use binary ORs (|) to
  # gradually add more and more observations to this filter out blacklist
  df['remove'] = np.repeat(False, df.shape[0])

  # run all the filters specified by the list in the input config file
  # all filter functions are passed df, and the run configuration
  # after each filter, append it onto the exclusion master list with a bitwise OR
  # if the filter function returns None, then just ignore it.
  for i, f in enumerate(filters):
    e = filter_funcs[f['name']](df, config, f)
    if e is not None:
      df['remove'] = (df['remove'] | e)

  return df

def process_files(config):

  # create our output data frames
  df_original = pd.DataFrame()
  df = pd.DataFrame()

  # iterate through each input file provided.
  for i, f in enumerate(config['input']):
    # first expand user or any vars
    f = os.path.expanduser(f)
    f = os.path.expandvars(f)

    logger.info('Reading in input file #{} | {} ...'.format(i+1, f))

    # load the input file with pandas
    # 
    # have a variable low memory option depending on the input type.
    # MaxQuant, for example, has a structure that forces pandas out of its
    # optimal low memory mode, and we have to specify it here.
    dfa = pd.read_csv(f, sep='\t', low_memory=config['low_memory'])

    # keep track of where observations came from. this is _not_ the raw file ID
    # but instead the ID from which input file it originated from, so that if
    # we need to split these observations up by input file in the future we can do so
    dfa['input_id'] = i

    # append a copy of dfa into df_original, because the conversion process will heavily
    # modify dfa. we need to keep a copy of the original dataframe in order to append
    # the new columns back onto it later.
    # re-index columns with '[dfa.columns.tolist()]' to preserve the general column order
    df_original = df_original.append(dfa, sort=True)[dfa.columns.tolist()]

    # if this input data already has DART-ID columns in it, then drop them,
    # since they cause problems later
    dart_cols = ['rt_minus', 'rt_plus', 'mu', 'muij', 'sigmaij', 'pep_new', 'exp_id', 'peptide_id', 'stan_peptide_id', 'exclude', 'residual', 'pep_updated', 'q-value']
    # print a warning if we see any
    if np.any(df_original.columns.isin(dart_cols)):
      logger.warning('Columns {} are recognized as DART-ID output columns. Removing these columns before proceeding. In the future, please input original input data files, not output files from DART-ID.'.format(np.array_str(df_original.columns[df_original.columns.isin(dart_cols)])))

      # drop existing dart cols
      for col in dart_cols:
        if col in df_original.columns:
          logger.debug('Removing column {}'.format(col))
          df_original = df_original.drop(col, axis=1)

    logger.info('Converting {} ({} PSMs)...'.format(f, dfa.shape[0]))

    # convert - takes subset of columns and renames them
    dfa = convert(dfa, config)

    # need to reset the input_id after the conversion process
    dfa['input_id'] = i
    # append to master dataframe
    df = df.append(dfa)

  # modify columns?
  # append the ion charge to the sequence
  # also make sure the charge column is specified and exists
  if config['add_charge_to_sequence'] and 'charge' in df.columns:
    logger.info('Appending charge to peptide sequence, to align different charge states separately.')
    df['sequence'] = df['sequence'] + '_' + df['charge'].apply(str)

  # create a unique ID for each PSM to help with stiching the final result together
  # after all of our operations
  df['id'] = range(0, df.shape[0])
  df_original['id'] = range(0, df.shape[0])

  # by default, exclude nothing from the original experiment
  df_original['remove'] = np.repeat(False, df_original.shape[0])

  # if the input already has an 'remove' column, then skip this step
  if 'remove' in config['col_names'] and config['col_names']['remove'] is not None:
    df['remove'] = df['remove'].astype(bool)
  else: # otherwise, run the filters
    df = filter_psms(df, config)

  # apply non-optional filters, PEP threshold and requirement that
  # sequence is observed in at least n experiments (num_experiments)
  
  # first check if PEP threshold is valid
  if 'pep_threshold' not in config or config['pep_threshold'] is None:
    raise ConfigFileError('PEP threshold not defined. Please provide a PEP threshold between 0 and 1 with the \"pep_threshold\" key. If including all PEPs, then set this value to 1.')
  # only allow values between 0 and 1
  if config['pep_threshold'] <= 0 or config['pep_threshold'] > 1:
    raise ConfigFileError('PEP threshold {} is incorrectly defined. Please provide a decimal number between 0.0 and 1.0.'.format(config['pep_threshold']))

  # remove any observations with null pep
  null_pep = pd.isnull(df['pep'])
  if np.sum(null_pep) > 0:
    df['remove'] = ((df['remove']) | (null_pep))
    logger.info('Removing {} PSMs with no PEP entry.'.format(np.sum(null_pep)))

  # check if num_experiments option is valid
  if 'num_experiments' not in config or config['num_experiments'] is None:
    raise ConfigFileError('No value provided to the filter for the number of raw files that a peptide must be observed in. Please provide an integer greater than or equal to 1 with the \"num_experiments\" key.')
  if type(config['num_experiments']) is not int or config['num_experiments'] < 1:
    raise ConfigFileError('Incorrect value of type {} provided to the filter for the number of raw files that a peptide must be observed in. Please provide an integer greater than or equal to 1 with the \"num_experiments\" key.'.format(str(type(config['num_experiments']))))
  num_exps = len(df['raw_file'].unique())
  if config['num_experiments'] > num_exps:
    raise ConfigFileError('Number of experiments filter threshold {} is greater than the number of experiments in the input list. Please provide an integer greater than or equal to 1 and less than the number of experiments with the \"num_experiments\" key.'.format(config['num_experiments']))
  if config['num_experiments'] == 1:
    logger.warning('Filter for number of raw files that a peptide must be observed in is set to 1. The alignment will proceed but this may result in non-informative canonical RTs and high residuals. It is recommended that this parameter is at least 3.')
  
  # count the number of experiments a peptide is observed in, but filter out
  # 1) PSMs removed from previous filters
  # 2) PSMs with PEP > pep_threshold
  exps_per_pep = df[-((df['remove']) | (df['pep'] >= config['pep_threshold']))].groupby('sequence')['raw_file'].unique().apply((lambda x: len(x)))
  # map values to DataFrame. peptides without any value will get NaN,
  # which will then be assigned to 0.
  exps_per_pep = df['sequence'].map(exps_per_pep)
  exps_per_pep[pd.isnull(exps_per_pep)] = 0

  # flag these sequences for removal as well
  logger.info('Removing {} PSMs from peptide sequences not observed confidently in more than {} experiments'.format(np.sum(exps_per_pep < config['num_experiments']), config['num_experiments']))
  df['remove'] = (df['remove'] | (exps_per_pep < config['num_experiments']))

  # check that every experiment has at least n PSMs available for alignment.
  # if not, then exclude them from alignment
  psms_per_exp = df.groupby('raw_file')['remove'].apply(lambda x: np.sum(x < config['pep_threshold']))
  exclude_exps = psms_per_exp.index.values[psms_per_exp < config['min_psms_per_experiment']]
  
  logger.warning('Experiments {} have < {} confident PSMs (PEP < {}) remaining after filtering. All PSMs belonging to these experiments will be excluded from the retention time alignment'.format(np.array_str(exclude_exps), config['min_psms_per_experiment'], config['pep_threshold']))

  # exclude experiments without enough PSMs
  df['remove'] = (df['remove'] | df['raw_file'].isin(exclude_exps))

  # recalculate exps_per_pep, since we removed some experiments and this
  # number will change based on the set of experiments we consider
  logger.info('Recalculating number of confident peptides across experiments...')

  exps_per_pep = df[-((df['remove']) | (df['pep'] >= config['pep_threshold']))].groupby('sequence')['raw_file'].unique().apply((lambda x: len(x)))
  exps_per_pep = df['sequence'].map(exps_per_pep)
  exps_per_pep[pd.isnull(exps_per_pep)] = 0
  logger.info('Additional {} PSMs from peptide sequences not observed confidently in more than {} experiments flagged for removal.'.format(np.sum(exps_per_pep < config['num_experiments']) - np.sum(df['remove']), config['num_experiments']))
  df['remove'] = (df['remove'] | (exps_per_pep < config['num_experiments']))

  ## --------------
  ## DONE FILTERING
  ## --------------

  # flag the observations in df_original that were removed
  df_original['remove'] = df['remove']
  # remove the flagged observations from the dataframe, and reset index
  df = df[df['remove']==False].reset_index(drop=True)

  # map peptide and experiment IDs
  # sort experiment IDs alphabetically - or else the order is by 
  # first occurrence of an observation of that raw file
  
  # if experiment or peptide IDs are already provided, then skip this step
  if 'exp_id' not in config['col_names'] or config['col_names']['exp_id'] is None:
    df['exp_id'] = df['raw_file'].map({
      ind: val for val, ind in enumerate(np.sort(df['raw_file'].unique()))})
  logger.info('{} experiments (raw files) loaded'.format(np.max(df['exp_id'])+1))

  if 'peptide_id' not in config['col_names'] or config['col_names']['peptide_id'] is None:
    df['peptide_id'] = df['sequence'].map({
      ind: val for val, ind in enumerate(df['sequence'].unique())})
  logger.info('{} peptide sequences loaded'.format(np.max(df['peptide_id'])+1))

  # EXCLUSION = PSM does not participate in alignment, but will participate in 
  #             confidence update since the PSM's associated peptide will get
  #             parameters from the alignment. 
  #             This is NOT the same as "remove", which means that the PSM's
  #             associated peptide does not have enough PSMs to participate
  #             in alignment and therefore receive parameters.
  
  # flag non-confident PSMs for exclusion from alignment process
  df['exclude'] = (df['pep'] >= config['pep_threshold'])
  logger.info('Excluding {} / {} ({:.2%}) PSMs from alignment process after filtering at PEP threshold of {}'.format(np.sum(df['pep'] >= config['pep_threshold']), df.shape[0], np.sum(df['pep'] >= config['pep_threshold']) / df.shape[0], config['pep_threshold']))

  # only take the four required columns (+ the IDs) with us
  # the rest were only needed for filtering and can be removed
  df = df[['sequence', 'raw_file', 'retention_time', 'pep', 
           'exp_id', 'peptide_id', 'input_id', 'id', 'exclude']]

  # sort by peptide_id, exp_id
  df = df.sort_values(['peptide_id', 'exp_id'])

  return df, df_original

def main():
  # load command-line args
  parser = argparse.ArgumentParser()  
  add_global_args(parser)
  args = parser.parse_args()

  # load config file
  # this function also creates the output folder
  config = read_config_file(args)

  # initialize logger
  init_logger(config['verbose'], os.path.join(config['output'], 'converter.log'), config['log_file'])

  # process all input files (converts and filters)
  df, df_original = process_files(config)
  
  #logger.info('{} / {} ({:.2%}) observations pass filters and will be used for alignment'.format(np.sum(~df['exclude']), 
  #    df_original.shape[0], np.sum(~df['exclude']) / df_original.shape[0]))

  # if no output is specified, throw an error
  if config['save_combined_output'] == False and config['save_separate_output'] == False:
    raise ConfigFileError('No output format specified. Either set \"save_combined_output\" to true, or set \"save_separate_output\" to true.')

  # write to file
  if config['save_combined_output']:
    # if combining input files, then write to one combined file
    out_path = os.path.join(config['output'], config['combined_output_name'])
    logger.info('Combining input file(s) and writing adjusted data file to {} ...'.format(out_path))
    df.to_csv(out_path, sep='\t', index=False)
  
  if config['save_separate_output']:
    # if keeping input files separate, then use 'input_id' to retain the
    # order in which the input files were passed in
    logger.info('Saving output to separate files...')
    for i, f in enumerate(config['input']):
      out_path = os.path.join(config['output'], 
        os.path.splitext(os.path.basename(f))[0] + config['output_suffix'] + \
          '_' + str(i) + '.txt')
      logger.info('Saving input file {} to {}'.format(i, out_path))
      df_a = df.loc[df['input_id'] == i]
      df_a.to_csv(out_path, sep='\t', index=False)

if __name__ == '__main__':
  main()
