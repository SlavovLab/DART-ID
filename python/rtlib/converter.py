#!/usr/bin/env python3
# coding: utf-8

# Converters from search-engine inputs to compatible outputs
# for alignment and update

import argparse
import json
import logging
import numpy as np
import os
import pandas as pd
import pkg_resources
import re
import sys

from functools import reduce
from rtlib.helper import intersect, add_version_arg, add_config_file_arg

logger = logging.getLogger()

# all filter functions take in the df, dfa, config object, and command line args
# as inputs, and output the exclude vector (True/False), where True means
# to exclude that particular row

def filter_uniprot_exclusion_list(df, dfa, config, args):
  """
  Filter proteins from exclusion list
  """
  prots_col = "proteins"
  # if proteins column is not present, then use the leading protein instead
  if "proteins" not in dfa.columns:
    prots_col = "leading_protein"

  exclusion_list = args.exclusion_list

  # parse exclusion list
  # if exclusion_list param is a path, then load the IDs from that path
  if exclusion_list is not None:
    # account for being provided an IO object
    if type(exclusion_list) is not str:
      exclusion_list = exclusion_list.name
    # load UniProt IDs from file line-by-line
    logger.info("Loading UniProt IDs from exclusion list {} ...".format(exclusion_list))
    exclusion_list = [line.rstrip('\n') for line in open(exclusion_list)]
    logger.info("{} proteins loaded from exclusion list.".format(len(exclusion_list)))
  else:
    exclusion_list = []

  # filter exclusion list
  if len(exclusion_list) > 0:    
    # we could only match the excluded IDs to the razor protein,
    # but we can be more strict and match the blacklisted IDs to the entire protein
    # string, containing all possible proteins
    pat = reduce((lambda x, y: x + "|" + y), exclusion_list)
    blacklist_filter = dfa["proteins"].str.contains(pat)
    blacklist_filter[pd.isnull(blacklist_filter)] = False

    logger.info("Filtering out {} PSMs from the exclusion list".format(np.sum(blacklist_filter)))
    return blacklist_filter
  else:
    logger.warning("Exclusion list empty or has incorrect format. Please list UniProt IDs separated by line breaks. Skipping exclusion list filter...")
    return None

def filter_contaminant(df, dfa, config, args):
  """
  Filter contaminants, as marked by the search engine
  Looking for a contaminant tag in the leading_protein column
  """
  CON_TAG = config["filters"]["contaminant"]["tag"]
  filter_con = dfa["proteins"].str.contains(CON_TAG)
  filter_con[pd.isnull(filter_con)] = False

  logger.info("Filtering out {} PSMs as contaminants with tag \"{}\"".format(np.sum(filter_con), CON_TAG))
  return filter_con

def filter_decoy(df, dfa, config, args):
  """
  Filter decoys, as marked by the search engine
  Looking for a decoy tag in the leading_protein column
  """
  REV_TAG = config["filters"]["decoy"]["tag"]
  filter_rev = dfa["leading_protein"].str.contains(REV_TAG)
  filter_rev[pd.isnull(filter_rev)] = False

  logger.info("Filtering out {} PSMs as decoys with tag \"{}\"".format(np.sum(filter_rev), REV_TAG))
  return filter_rev

def filter_retention_length(df, dfa, config, args):
  """
  Filter by retention length, which is a measure of the peak width
  during chromatography.
  """
  filter_retention_length = args.filter_retention_length

  if filter_retention_length is None or filter_retention_length == 0:
    logger.info("Retention length filter not defined. Skipping retention length filter...")
    return None

  if type(filter_retention_length) is not float or filter_retention_length < -1 or filter_retention_length > np.max(dfa["retention_time"]):
    logger.warning("Retention length filter {} is not defined or incorrectly defined. Please provide a decimal number between 0.0 and max(RT). If set to 0, then this filter will be skipped. If set to -1, then this filter will be set to its default, max(RT) / 60. Skipping retention length filter...".format(filter_retention_length))
    return None

  # default is -1, which means that we will set it to max(RT) / 60
  if filter_retention_length == -1.0:
    filter_retention_length = np.max(dfa["retention_time"]) / 60

  filter_rtl = (df[config["col_names"]["retention_length"]] > filter_retention_length)

  logger.info("Filtering out {} PSMs with retention length greater than {:.2f}".format(np.sum(filter_rtl), filter_retention_length))
  return filter_rtl

def filter_pep(df, dfa, config, args):
  """
  Filter by PEP (Posterior Error Probability), 
  measured from spectra and a search engine
  """
  filter_pep = args.filter_pep

  if filter_pep is None or filter_pep == -1:
    logger.info("PEP filter not defined. Skipping PEP filter...")
    return None

  if type(filter_pep) is not float or filter_pep <= 0 or filter_pep > 1:
    logger.warning("PEP filter {} is not defined or incorrectly defined. Please provide a decimal number between 0.0 and 1.0. Skipping PEP filter...".format(filter_pep))
    return None

  filter_pep = (dfa["pep"] > args.filter_pep)

  logger.info("Filtering out {} PSMs with PEP greater than {:.2f}".format(np.sum(filter_pep), args.filter_pep))
  return filter_pep


filter_funcs = {
  "uniprot_exclusion":  filter_uniprot_exclusion_list,
  "contaminant":        filter_contaminant,
  "decoy":              filter_decoy,
  "retention_length":   filter_retention_length,
  "pep":                filter_pep
}


def convert(df, config, args):

  # first load the required sequence/modified sequence, raw file, retention time,
  # and pep column into a new dataframe, dfa
  col_names = config["col_names"]

  # use canonical sequence, or modified sequence?
  seq_column = "sequence"
  # make sure the modified sequence column exists as well
  if not args.use_unmodified_sequence:
    if type(col_names["modified_sequence"]) is str:
      seq_column = "modified_sequence"
    else:
      raise ValueError("Modified Sequence selected but either input file type does not support modified sequences, or the input file type is missing the modified sequence column.")
  else:
    logger.info("Using unmodified peptide sequence instead of modified peptide sequence")

  # get the four required columns from the input config, and load into dfa
  cols = [col_names[seq_column], col_names["raw_file"], col_names["retention_time"], col_names["pep"]]
  dfa = df[cols]
  # rename columns
  dfa = dfa.rename(columns={dfa.columns[0]:"sequence", dfa.columns[1]:"raw_file", dfa.columns[2]:"retention_time", dfa.columns[3]:"pep"})

  ## Begin Filtering ----------
  
  # by default, exclude nothing. we'll use binary ORs (|) to
  # gradually add more and more observations to this exclude blacklist
  exclude = np.repeat(False, df.shape[0])

  # load the filtering functions specified by the input config
  filters = list(config["filters"].keys())

  # each filter has a specified required column from the dataframe
  # make sure these columns exist before proceeding
  for f in filters:
    f_obj = config["filters"][f]
    # the required_cols options is allowed to be empty or non-existent
    if "required_cols" not in f_obj or (type(f_obj["required_cols"]) is list and len(f_obj["required_cols"]) == 0):
      continue
    for i in f_obj["required_cols"]:
      if col_names[i] not in df.columns:
        raise ValueError("Filter {} required a column {}, but this was not found in the input dataframe.".format(f, i))

  # for some filters, we're going to need the leading protein or protein column
  # if so, then pre-emptively load the protein into df. we'll remove it later
  # after the filtering processes
  filters_with_proteins = ["uniprot_exclusion", "contaminant", "decoy"]
  if len(intersect(filters, filters_with_proteins)) > 0:
    # grab the leading razor protein
    if type(col_names["leading_protein"]) is str:
      dfa["leading_protein"] = df[col_names["leading_protein"]]
    # grab all proteins
    if type(col_names["proteins"]) is str:
      dfa["proteins"] = df[col_names["proteins"]]
  
  # run all the filters specified by the list in the input config file
  # all filter functions are passed df, dfa, the input config, and command-line args
  # after each filter, append it onto the exclusion master list with a bitwise OR
  # if the filter function returns None, then just ignore it.
  for f in filters:
    e = filter_funcs[f](df, dfa, config, args)
    if e is not None:
      exclude = (exclude | e)

  dfa["exclude"] = exclude

  return dfa


def process_files(args):

  # load input file types
  input_types = pkg_resources.resource_stream("rtlib", "/".join(("config", "input_types.json")))
  input_types = json.load(input_types)["input_types"]
  
  # get the input config of the specified file type
  config = input_types[args.type]
  logger.info("Parsing input files as the output from {}".format(config["name"]))

  # create our output data frames
  df_original = pd.DataFrame()
  df = pd.DataFrame()

  # iterate through each input file provided.
  for i, f in enumerate(args.input):
    logger.info("Reading in input file #{} | {} ...".format(i, f.name))

    # load the input file with pandas
    # 
    # have a variable low memory option depending on the input type.
    # MaxQuant, for example, has a structure that forces pandas out of its
    # optimal low memory mode, and we have to specify it here.
    dfa = pd.read_csv(f, sep="\t", low_memory=config["low_memory"])

    # keep track of where observations came from. this is _not_ the raw file ID
    # but instead the ID from which input file it originated from, so that if
    # we need to split these observations up by input file in the future we can do so
    dfa["input_id"] = i

    # append a copy of dfa into df_original, because the conversion process will heavily
    # modify dfa. we need to keep a copy of the original dataframe in order to append
    # the new columns back onto it later.
    df_original = df_original.append(dfa)

    # check whether the columns specified in the input types configuration
    # are present in the loaded dataframe
    col_names = config["col_names"]
    for j in list(col_names.keys()):
      col = col_names[j]

      # prep our not found error
      not_found_error = ValueError("Column \"{}\" of column type \"{}\" not present in input file of type \"{}\". Please check that the input file is of the specified type. If column names are still not matching up, contact the package author.".format(col, j, config["name"]))

      # each column can either be a:
      # - string
      # - list of strings
      # - empty list
      if type(col) is str:
        # if the column name is a string, then there is only one possible
        # column name in the dataframe for this column.
        if col not in dfa.columns:
          raise not_found_error

      elif type(col) is list and len(col) >  0:
        # if the column name is a list of strings, then check each one
        # if one is found, then overwrite the config dict's col_name list
        # with the found string
        # not going to test if two or more match here. the first one will 
        # be taken and that's it
        found = False
        for col_ in col:
          if col_ in dfa.columns:
            found = True
            config["col_names"][j] = col_
        if not found:
          raise not_found_error

      elif type(col) is list and len(col) == 0:
        # if the column name is an empty list, then
        # it doesn't exist in this input type. skip...
        pass
      
      else:
        # there's something wrong with the input types config file
        raise ValueError("Unknown column configuration in input_types.json. Please contact the package author.")

    logger.info("Converting {} ({} PSMs)...".format(f.name, dfa.shape[0]))

    dfa = convert(dfa, config, args)

    dfa["input_id"] = i
    df = df.append(dfa)

  
  # create a unique ID for each PSM to help with stiching the final result together
  # after all of our operations
  df["id"] = range(0, df.shape[0])
  df_original["id"] = range(0, df.shape[0])

  # remove experiments from blacklist
  # by default, exclude nothing from the original experiment
  df_original["input_exclude"] = np.repeat(False, df_original.shape[0])
  if args.remove_exps is not None and len(args.remove_exps) > 0:
    exclude_exps = list(filter(lambda x: re.search(r""+args.remove_exps+"", x), df["raw_file"].unique()))
    logger.info("Filtering out {} observations matching \"{}\"".format(np.sum(df["raw_file"].isin(exclude_exps)), args.remove_exps))

    # remove excluded rows from the sparse dataframe,
    # but keep them in the original data frame, so that we can stitch together
    # the final output later
    exclude_exps = df["raw_file"].isin(exclude_exps).values
    df = df[~exclude_exps]

    # keep track of which experiments were excluded in 
    df_original["input_exclude"] = exclude_exps
  
  df = df.reset_index(drop=True)

  # map peptide and experiment IDs
  # sort experiment IDs alphabetically - or else the order is by 
  # first occurrence of an observation of that raw file
  df["exp_id"] = df["raw_file"].map({ind: val for val, ind in enumerate(np.sort(df["raw_file"].unique()))})
  df["peptide_id"] = df["sequence"].map({ind: val for val, ind in enumerate(df["sequence"].unique())})

  # filter for occurence in number of experiments
  if args.filter_num_exps is not None and args.filter_num_exps >= 2:
    # only want to do this operation on PSMs that aren't already marked to be filtered out
    # first, take subset on remaining PSMs, then count number 
    # of unique experiments for each of them
    exps_per_pep = df[-(df["exclude"])].groupby("peptide_id")["exp_id"].unique().apply((lambda x: len(x)))
    # map values to DataFrame. peptides without any value will get NaN,
    # which will then be assigned to 0.
    exps_per_pep = df["peptide_id"].map(exps_per_pep)
    exps_per_pep[np.isnan(exps_per_pep)] = 0

    logger.info("Filtering out {} PSMs that have less than {} occurrences in different experiments.".format((exps_per_pep < args.filter_num_exps).sum(), args.filter_num_exps))
    df["exclude"] = (df["exclude"] | (exps_per_pep < args.filter_num_exps))

  if args.filter_smear_threshold != 0:
    # filter out "smears". even confidently identified PSMs can have bad chromatography,
    # and in that case it is unproductive to include them into the alignment.
    # 
    # TODO: there might also be merit to excluding these observations from the PEP update
    # process as well, given that the spectral PEP is below a very 
    # conservative threshold (1% or maybe even lower)
    logger.info("Determining RT spread of peptides within each experiment...")
    # for each experiment-peptide pair, get the range of retention times
    smears = df.groupby(["exp_id", "peptide_id"])["retention_time"].apply(lambda x: np.ptp(x))
    # get the (exp_id, peptide_id) tuples for PSMs with a range above the threshold
    max_rts = df.groupby("exp_id")["retention_time"].max().values
    smears = smears[smears > max_rts[smears.index.to_frame()["exp_id"].values] * args.filter_smear_threshold].index.values

    # map the tuples back to the original data frame, and set smears to be excluded
    smears = pd.Series(list(zip(df["exp_id"], df["peptide_id"]))).isin(smears)
    logger.info("Filtering out {} PSMs with an intra-experiment RT spread of above max(exp_RT) / {} minutes".format(smears.sum(), args.filter_smear_threshold))
    df["exclude"] = (df["exclude"] | (smears))

  # sort by peptide_id, exp_id
  df = df.sort_values(["peptide_id", "exp_id"])

  return df, df_original
  

def add_converter_args(parser):
  """
  Tack on arguments needed to parse search-engine input file

  Parameters
  ----------
  parser: argparser object

  Returns
  -------
  argparser object

  """
  parser.add_argument("-t", "--type", type=str, default="RTLib", 
    choices=["MQ", "PD", "RTLib"],
    help="""
    Search engine type. e.g.:
    MQ    (MaxQuant), 
    PD    (ProteomeDiscoverer),
    RTLib (RTLib format - from the "convert" command)
    Default: RTLib
    """)
  parser.add_argument("-v", "--verbose", action="store_true", default=False)
  parser.add_argument("--filter-retention-length", type=float, default=-1.0, help="Filter PSMs by retention time length (Peak width, in minutes). Default: max(RT) / 60")
  parser.add_argument("--filter-pep", type=float, default=0.5, help="Filter PSMs by PEP (Posterior Error Probability). Default: 0.5")
  parser.add_argument("--filter-num-exps", type=int, default=3, help="Filter by occurence of peptide in number of experiments. i.e., if a peptide is not observed in at least this number of experiments over the entire set, then it will be filtered out before alignment. Default: 3")
  parser.add_argument("--filter-smear-threshold", type=float, default=0.03, help="Filter out peptides that have a intra-experiment RT range of this number * max(experiment_RT). Set to 0 to skip this step.")
  parser.add_argument("-e", "--exclusion-list", type=argparse.FileType("r"), help="Path to exclusion list file - UniProt IDs separated by lines. Excluding contaminants that are missed by search engines can be crucial for alignment. If STAN is exceeding a reasonable iteration limit, or if the alignment residuals are unexpectedly large, consider filtering out proteins with peptides with inconsistent RTs.")
  parser.add_argument("--remove-exps", type=str, default=None, help="Regular expression of experiments (raw files) to remove from the data. Regular expression is case-sensitive. Default: None")
  parser.add_argument("-m", "--use-unmodified-sequence", action="store_true", default=False, help="Use unmodified sequence instead of modified peptide sequence. Default: False")


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("input", type=argparse.FileType("r"), nargs="+", help="Input file(s) from search engine")
  parser.add_argument("-o", "--output", help="Path to converted file. Default: prints to stdout")

  add_converter_args(parser)
  add_version_arg(parser)
  add_config_file_arg(parser)

  args = parser.parse_args()

  # set up logger
  for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler) 
   
  logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
  logger = logging.getLogger()

  if args.verbose: logger.setLevel(logging.DEBUG)
  else: logger.setLevel(logging.WARNING)

  consoleHandler = logging.StreamHandler(sys.stdout)
  consoleHandler.setFormatter(logFormatter)
  logger.addHandler(consoleHandler)

  logger.info(" ".join(sys.argv[1:]))

  df, df_original = process_files(args)
  
  logger.info("{} / {} ({:.2%}) observations pass criteria and will be used for alignment".format(df.shape[0] - df["exclude"].sum(), df.shape[0], (df.shape[0] - df["exclude"].sum()) / df.shape[0]))
  
  logger.info("Saving converted data to {} ...".format(args.output))
  df.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
  main()



