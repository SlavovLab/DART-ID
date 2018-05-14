#!/usr/bin/env python3
# coding: utf-8

# Converters from search-engine inputs to compatible outputs
# for alignment and update

import argparse
import logging
import numpy as np
import os
import pandas as pd
import re
import sys

from functools import reduce

logger = logging.getLogger()

def convert_mq(df, use_unmodified_sequence=False, filter_pep=0.5, filter_retention_length=-1, include_contaminants=False, include_decoys=False, exclusion_list=[]):
  """
  Parameters
  ----------
  df: pandas.DataFrame
    Input from pd.readcsv(path, sep='\t'), where path is the path of
    the input MaxQuant evidence file

  use_unmodified_sequence: bool (Default: False)
    Use the modified peptide sequence instead of the canonical sequence.
    Assumes that TMT tags are excluded from modifications (modifications
    are Oxidation, Methylation, etc.)

  filter_pep: float (Default: 0.5)
    Mark PSMs with PEP > filter_pep for filtering.
    Set to None to ignore PEP

  filter_retention_length: int (Default: -1 --> max(RT) / 60)
    Mark PSMs with retention length > filter_retention_length for filtering. 
    Set to None to ignore retention length

  include_contaminants: bool (Default: Fakse)
    Include contaminants (CON__ proteins).

  include_decoys: bool (Default: False)
    Include decoys (REV__ proteins). Don't know why this would
    ever be set to true, but here's the option.

  exclusion_list: list|string (Default: [])
    Blacklist of UniProt protein IDs to mark for filtering.
    If exclusion_list is a string, will treat as a path and load UniProt IDs
    from file. Protein IDs should be separated by lines

  Returns
  -------
  pandas.DataFrame
    Truncated, processed DataFrame

  """

  # first set all column names to lowercase so we avoid upper/lowercase
  # shenanigans between MaxQuant output versions
  # replace spaces with underscores
  df.columns = [re.sub(r"\s", "_", str.lower(c)) for c in df.columns]

  # use canonical sequence, or modified sequence?
  seq_column = "modified_sequence"
  if use_unmodified_sequence: seq_column = "sequence"

  # grab only the columns we need for the alignment process
  dfa = df[[seq_column, "raw_file", "retention_time", "pep"]]

  # alias modified sequence as sequence moving forwards
  if use_unmodified_sequence:
    logger.info("Using unmodified peptide sequence instead of modified peptide sequence")
  else:
    dfa = dfa.rename(columns={ "modified_sequence": "sequence" })

  ## Begin Filtering ----------
  
  # by default, exclude nothing. we'll use binary ORs (|) to
  # gradually add more and more observations to this exclude blacklist
  dfa["exclude"] = False
  
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


  # see if the user even needs to filter the proteins out from the data
  # if not, don't even bother pulling out the protein data
  if len(exclusion_list) > 0 or filter_contaminants or filter_decoys:

    # grab the leading razor protein for filtering
    prots = []
    if "leading_razor_protein" in df.columns:
      prots = df["leading_razor_protein"]
    elif "proteins" in df.columns:
      # if the leading razor protein column does not exist, then take
      # the first protein from the Proteins column
      logger.info("No razor proteins column found. Extracting razor proteins from 'Proteins' column")
      prots = [x.split(';')[0] if type(x) is str else "" for x in ev.Proteins]
    else:
      raise IOError("No protein column found. Please provide file with either \"Leading razor protein\", or \"Proteins\" as a column.")
      return None

    # contaminant and decoy tags
    # not sure if these will ever change between MaxQuant versions,
    # but will leave them as vars if they do
    CON_TAG = "CON*"
    REV_TAG = "REV*"

    # convert prots to series
    prots = pd.Series(prots)

    # filter contaminants and decoys from the MaxQuant evidence file
    # these should always be on
    if not include_contaminants:
      logger.info("Filtering contaminants - with tag: {}".format(CON_TAG))
      filter_con = prots.str.contains(CON_TAG)
      filter_con[pd.isnull(filter_con)] = False
      dfa["exclude"] = (dfa["exclude"] | filter_con)
    if not include_decoys:
      logger.info("Filtering decoys - with tag: {}".format(REV_TAG))
      filter_rev = prots.str.contains(REV_TAG)
      filter_rev[pd.isnull(filter_rev)] = False
      dfa["exclude"] = (dfa["exclude"] | filter_rev)

    # filter exclusion list
    if len(exclusion_list) > 0:
      logger.info("Filtering PSMs with proteins in exclusion list...")
      
      # we could only match the excluded IDs to the razor protein,
      # but we can be more strict and match the blacklisted IDs to the entire protein
      # string, containing all possible proteins
      pat = reduce((lambda x, y: x + "|" + y), exclusion_list)
      blacklist_filter = df["proteins"].str.contains(pat)
      blacklist_filter[pd.isnull(blacklist_filter)] = False
      dfa["exclude"] = (dfa["exclude"] | blacklist_filter)


  # filter for retention length, if specified
  # default is -1, which means that we will set it to max(RT) / 60
  if filter_retention_length is -1:
    filter_retention_length = np.max(dfa["retention_time"]) / 60
  if filter_retention_length is not None or filter_retention_length > 0:
    logger.info("Filtering PSMs with Retention Length greater than {:.2f} ...".format(filter_retention_length))
    dfa["exclude"] = (dfa["exclude"] | (df["retention_length"] > filter_retention_length))

  # filter for PEP, if specified
  if filter_pep is not None and filter_pep > 0:
    logger.info("Filtering PSMs with PEP greater than {} ...".format(filter_pep))
    dfa["exclude"] = (dfa["exclude"] | (df["pep"] > filter_pep))

  return dfa


def convert_pd(df):
  """
  Parameters
  ----------
  df: pandas.DataFrame

  """
  return df

def process_files(args):
  df_original = pd.DataFrame()
  df = pd.DataFrame()
  for i, f in enumerate(args.input):
    logger.info("Reading in input file {} ...".format(f.name))
    dfa = pd.read_csv(f, sep="\t", low_memory=False)

    # keep track of where observations came from
    dfa["input_id"] = i
    df_original = df_original.append(dfa)

    logger.info("Converting " + f.name + " ...")
    if args.type == None:
      logger.info("Using pre-converted input file for alignment...")
      # check if input file has the necessary columns
      cols = ["sequence", "raw_file", "retention_time", "pep", "exclude", "exp_id", "peptide_id"]
      if np.sum(np.isin(cols, dfa.columns.values)) != len(cols):
        # we're missing some columns here. break out
        raise IOError("Input file missing required columns. Did you forget to specify the input file type?")
        return df
    elif args.type == "MQ":
      dfa = convert_mq(dfa, use_unmodified_sequence=args.use_unmodified_sequence, filter_retention_length=args.filter_retention_length, filter_pep=args.filter_pep, include_contaminants=args.include_contaminants, include_decoys=args.include_decoys, exclusion_list=args.exclusion_list)
    elif args.type == "PD":
      raise NotImplementedError("ProteomeDiscoverer not supported yet.")
    else:
      raise ValueError("""
        Input type not recognized. Acceptable types are:\n
        'MQ' - MaxQuant\n
        'PD' - ProteomeDiscoverer\n
        """)

    # keep track of where observations came from
    dfa["input_id"] = i
    df = df.append(dfa)

  # create a unique ID for each PSM to help with stiching the final result together
  # after all of our operations
  df["id"] = range(0, df.shape[0])
  df_original["id"] = range(0, df.shape[0])

  # remove experiments from blacklist
  # by default, exclude nothing from the original experiment
  df_original["exclude"] = False
  if args.remove_exps is not None and len(args.remove_exps) > 0:
    exclude_exps = list(filter(lambda x: re.search(r""+args.remove_exps+"", x), df["raw_file"].unique()))
    logger.info("Removing {} observations matching {} ...".format(np.sum(df["raw_file"].isin(exclude_exps)), args.remove_exps))

    # remove excluded rows from the sparse dataframe,
    # but keep them in the original data frame, so that we can stitch together
    # the final output later
    exclude_exps = df["raw_file"].isin(exclude_exps).values
    df = df[~exclude_exps]

    # keep track of which experiments were excluded in 
    df_original["exclude"] = exclude_exps
  
  df = df.reset_index(drop=True)

  # map peptide and experiment IDs
  # sort experiment IDs alphabetically - or else the order is by 
  # first occurrence of an observation of that raw file
  df["exp_id"] = df["raw_file"].map({ind: val for val, ind in enumerate(np.sort(df["raw_file"].unique()))})
  df["peptide_id"] = df["sequence"].map({ind: val for val, ind in enumerate(df["sequence"].unique())})

  # filter for occurence in number of experiments
  if args.filter_num_exps is not None and args.filter_num_exps >= 2:
    logger.info("Removing peptides with less than {} occurrences in different experiments...".format(args.filter_num_exps))
    # only want to do this operation on PSMs that aren't already marked to be filtered out
    # first, take subset on remaining PSMs, then count number 
    # of unique experiments for each of them
    exps_per_pep = df[-(df["exclude"])].groupby("peptide_id")["exp_id"].unique().apply((lambda x: len(x)))
    # map values to DataFrame. peptides without any value will get NaN,
    # which will then be assigned to 0.
    exps_per_pep = df["peptide_id"].map(exps_per_pep)
    exps_per_pep[np.isnan(exps_per_pep)] = 0
    df["exclude"] = (df["exclude"] | (exps_per_pep < args.filter_num_exps))

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
  parser.add_argument("-t", "--type", type=str, default=None, 
    choices=["MQ", "PD"],
    help="""
    Search engine type. e.g.:
    MQ (MaxQuant), 
    PD (ProteomeDiscoverer).
    Default: None (Input assumed to be already converted)
    """)
  parser.add_argument("-v", "--verbose", action="store_true", default=False)
  parser.add_argument("--include-contaminants", action="store_true", default=False, help="Filter contaminant proteins, as marked so by the search engine. Default: False")
  parser.add_argument("--include-decoys", action="store_true", default=False, help="Filter decoy matches. Default: False")
  parser.add_argument("--filter-retention-length", type=int, default=-1, help="Filter PSMs by retention time length (Peak width, in minutes). Default: max(RT) / 60")
  parser.add_argument("--filter-pep", type=float, default=0.5, help="Filter PSMs by PEP (Posterior Error Probability). Default: 0.5")
  parser.add_argument("--filter-num-exps", type=int, default=3, help="Filter by occurence of peptide in number of experiments. i.e., if a peptide is not observed in at least this number of experiments over the entire set, then it will be filtered out before alignment. Default: 3")
  parser.add_argument("-e", "--exclusion-list", type=argparse.FileType("r"), help="Path to exclusion list file - UniProt IDs separated by lines. Excluding contaminants that are missed by search engines can be crucial for alignment. If STAN is exceeding a reasonable iteration limit, or if the alignment residuals are unexpectedly large, consider filtering out proteins with peptides with inconsistent RTs.")
  parser.add_argument("--remove-exps", type=str, default=None, help="Regular expression of experiments (raw files) to remove from the data. Regular expression is case-sensitive. Default: None")
  parser.add_argument("-m", "--use-unmodified-sequence", action="store_true", default=False, help="Use unmodified sequence instead of modified peptide sequence. Default: False")


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("input", type=argparse.FileType("r"), nargs="+", help="Input file(s) from search engine")
  parser.add_argument("-o", "--output", help="Path to converted file. Default: prints to stdout")

  add_converter_args(parser)

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

  #if args.type == None:
  #  raise TypeError("No input file type specified. Is the input already in the converted format?")

  df, df_original = process_files(args)

  logger.info("{} / {} ({:.2%}) observations pass criteria and will be used for alignment".format(df.shape[0] - df["exclude"].sum(), df.shape[0], (df.shape[0] - df["exclude"].sum()) / df.shape[0]))
  
  logger.info("Saving converted data to {} ...".format(args.output))
  df.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
  main()



