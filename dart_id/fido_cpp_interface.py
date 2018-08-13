#!/usr/bin/env python3
# coding: utf-8

## ===============================
## NOT CURRENTLY BEING USED
## ===============================

import argparse
import importlib
import logging
import numpy as np
import os
import pandas as pd
import re

from dart_id.helper import *

logger = logging.getLogger('root')

def write_outputs(df, out_folder):

  #with open(out_file, 'w') as f:

  graph_text = ''
  targets = []
  decoys = []

  prot_list = []
  prot_cols = ['Leading proteins', 'Proteins']

  for i in range(0, df.shape[0]):

    prot_list = []
    for col in prot_cols:
      prots = df[col].loc[i]
      if pd.isnull(prots): continue
      #prot_list = prot_list + [str.join('_', prot.split('|')[0:2]) for prot in prots.split(';')]
      prot_list = prot_list + [prot for prot in prots.split(';')]

    prot_list = np.unique(prot_list)
    if len(prot_list) < 1: continue

    graph_text += ('e ' + df['Sequence'].loc[i] + '\n')

    for prot in prot_list:
      graph_text += ('r ' + prot + '\n')

      if 'REV__' in prot:
        # add to decoy list
        decoys.append(prot)
      else:
        # add to targets list
        targets.append(prot)
        

    graph_text += ('p ' + str(1-df['PEP'].loc[i]) + '\n')
    
  targets = np.unique(targets)
  decoys = np.unique(decoys)

  # write to strings
  target_decoy_text = '{ '
  for i, t in enumerate(targets):  
    if i == len(targets)-1:
      target_decoy_text += (t + ' }')
    else: 
      target_decoy_text += (t + ' , ')

  target_decoy_text += ('\n{ ')
  for i, d in enumerate(decoys):
    if i == len(targets)-1:
      target_decoy_text += (d + ' }')
    else: 
      target_decoy_text += (d + ' , ')


  # write psm graph
  with open(os.path.join(out_folder, 'psm_graph.txt'), 'w') as f:
    f.write(graph_text)
    f.close()

  # write target decoy file
  with open(os.path.join(out_folder, 'target_decoy.txt'), 'w') as f:
    f.write(target_decoy_text)
    f.close()



def main():
  # load command-line args
  parser = argparse.ArgumentParser()  
  add_global_args(parser, add_config_file=False)
  #parser.add_argument()
  args = parser.parse_args()

  # initialize logger
  init_logger(args.verbose, '', log_to_file=False)
  logger = logging.getLogger('root')

  df = pd.DataFrame()

  # iterate through each input file provided.
  for i, f in enumerate(args.input):
    logger.info('Reading in input file #{} | {} ...'.format(i+1, f.name))

    # load the input file
    dfa = pd.read_csv(f.name, sep='\t', low_memory=False)
    # append a copy of dfa into df, because the conversion process will heavily
    # modify dfa. we need to keep a copy of the original dataframe in order to append
    # the new columns back onto it later.
    # re-index columns with '[dfa.columns.tolist()]' to preserve the general column order
    df = df.append(dfa)[dfa.columns.tolist()]

  df = df.sort_values(by=['Sequence']).reset_index(drop=True)

  # write psm graph file
  write_outputs(df, args.output)
  

if __name__ == '__main__':
  main()