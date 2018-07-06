#!/usr/bin/env python3
# coding: utf-8

import argparse
import importlib
import logging
import numpy as np
import os
import pandas as pd
import re

from dart_id.helper import *

logger = logging.getLogger('root')

def main():
  # load command-line args
  parser = argparse.ArgumentParser()  
  add_global_args(parser, add_config_file=False)

  parser.add_argument('--include', type=str, default=None, help='Inclusion expression (regular expression) that is evaluated against all raw file names to determine whether or not that raw file is included.')
  parser.add_argument('--raw-file-col', type=str, default='Raw file', help='Name of the raw file column.')
  args = parser.parse_args()

  # initialize logger
  init_logger(args.verbose, '', log_to_file=False)
  logger = logging.getLogger('root')

  df = pd.DataFrame()

  # iterate through each input file provided.
  for i, f in enumerate(args.input):

    logger.info('Reading in input file #{} | {} ...'.format(i+1, f.name))

    # load the input file with pandas
    # 
    # have a variable low memory option depending on the input type.
    # MaxQuant, for example, has a structure that forces pandas out of its
    # optimal low memory mode, and we have to specify it here.
    dfa = pd.read_csv(f.name, sep='\t', low_memory=False)

    # keep track of where observations came from. this is _not_ the raw file ID
    # but instead the ID from which input file it originated from, so that if
    # we need to split these observations up by input file in the future we can do so
    dfa['input_id'] = i

    # append a copy of dfa into df, because the conversion process will heavily
    # modify dfa. we need to keep a copy of the original dataframe in order to append
    # the new columns back onto it later.
    # re-index columns with '[dfa.columns.tolist()]' to preserve the general column order
    df = df.append(dfa)[dfa.columns.tolist()]

  logger.info('Raw files:')
  logger.info(df[args.raw_file_col].unique())

  if args.include is not None:
    if len(args.include) == 0:
      raise Exception('Experiment inclusion expression by raw file name provided, but expression is defined incorrectly.')
    else:
      # see if any raw file names match the user-provided expression
      include_exps = list(filter(lambda x: re.search(r'' + args.include + '', x), df[args.raw_file_col].unique()))

      logger.info('{} raw files match the inclusion expression \"{}\"'.format(len(include_exps), args.include))
      logger.info(include_exps)
      logger.info('Retaining {} observations out of {} total'.format(np.sum(df[args.raw_file_col].isin(include_exps)), df.shape[0]))

      include_exps = df[args.raw_file_col].isin(include_exps).values
      df = df.loc[include_exps].reset_index(drop=True)

  # if combining input files, then write to one combined file
  logger.info('Combining input file(s) and writing adjusted data file to {} ...'.format(args.output))
  df.to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
  main()