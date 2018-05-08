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

verbose = False

def convert_mq(df, use_modified_sequence=True, filter_pep=0.5, filter_retention_length=5, filter_contaminants=True, filter_decoys=True, exclusion_list=[]):
  """
  Parameters
  ----------
  df: pandas.DataFrame
    Input from pd.readcsv(path, sep='\t'), where path is the path of
    the input MaxQuant evidence file

  use_modified_sequence: bool (Default: True)
    Use the modified peptide sequence instead of the canonical sequence.
    Assumes that TMT tags are excluded from modifications (modifications
    are Oxidation, Methylation, etc.)

  filter_pep: float (Default: 0.5)
    Mark PSMs with PEP > filter_pep for filtering.
    Set to None to ignore PEP

  filter_retention_length: int (Default: 5)
    Mark PSMs with retention length > filter_retention_length for filtering. 
    Set to None to ignore retention length

  filter_contaminants: bool (Default: True)
    Mark contaminants (CON__ proteins) for filtering.

  filter_decoys: bool (Default: True)
    Mark decoys (REV__ proteins) for filtering. Don't know why this would
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
  seq_column = "sequence"
  if use_modified_sequence: seq_column = "modified_sequence"

  # grab only the columns we need for the alignment process
  dfa = df[[seq_column, "raw_file", "retention_time", "pep"]]

  # alias modified sequence as sequence moving forwards
  if use_modified_sequence:
    logging.info("Using modified peptide sequence instead of peptide sequence")
    dfa = dfa.rename(columns={"modified_sequence":"sequence"})

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
    logging.info("Loading UniProt IDs from exclusion list " + exclusion_list + " ...")
    exclusion_list = [line.rstrip('\n') for line in open(exclusion_list)]
    logging.info(str(len(exclusion_list)) + " proteins loaded from exclusion list.")
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
      logging.info("No razor proteins column found. Extracting razor proteins from 'Proteins' column")
      prots = [x.split(';')[0] if type(x) is str else "" for x in ev.Proteins]
    else:
      logging.error("No protein column found. Please provide file with either \"Leading razor protein\", or \"Proteins\" as a column.")
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
    if filter_contaminants:
      logging.info("Filtering contaminants - with tag: " + CON_TAG)
      filter_con = prots.str.contains(CON_TAG)
      filter_con[pd.isnull(filter_con)] = False
      dfa["exclude"] = (dfa["exclude"] | filter_con)
    if filter_decoys:
      logging.info("Filtering decoys - with tag: " + REV_TAG)
      filter_rev = prots.str.contains(REV_TAG)
      filter_rev[pd.isnull(filter_rev)] = False
      dfa["exclude"] = (dfa["exclude"] | filter_rev)

    # filter exclusion list
    if len(exclusion_list) > 0:
      logging.info("Filtering PSMs with proteins in exclusion list...")
      
      # we could only match the excluded IDs to the razor protein,
      # but we can be more strict and match the blacklisted IDs to the entire protein
      # string, containing all possible proteins
      pat = reduce((lambda x, y: x + "|" + y), exclusion_list)
      blacklist_filter = df["proteins"].str.contains(pat)
      blacklist_filter[pd.isnull(blacklist_filter)] = False
      dfa["exclude"] = (dfa["exclude"] | blacklist_filter)


  # filter for retention length, if specified
  if filter_retention_length is not None or filter_retention_length > 0:
    logging.info("Filtering PSMs with Retention Length greater than " + str(filter_retention_length) + " ...")
    dfa["exclude"] = (dfa["exclude"] | (df["retention_length"] > filter_retention_length))

  # filter for PEP, if specified
  if filter_pep is not None and filter_pep > 0:
    logging.info("Filtering PSMs with PEP greater than " + str(filter_pep) + " ...")
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
    logging.info("Reading in input file " + f.name + " ...")
    dfa = pd.read_csv(f, sep="\t", low_memory=False)

    df_original = df_original.append(dfa)

    logging.info("Converting " + f.name + " ...")
    if args.type == None:
      logging.info("Using pre-converted input file for alignment...")
      # check if input file has the necessary columns
      cols = ["sequence", "raw_file", "retention_time", "pep", "exclude", "exp_id", "peptide_id"]
      if np.sum(np.isin(cols, dfa.columns.values)) != len(cols):
        # we're missing some columns here. break out
        logging.error("Input file missing required columns. Did you forget to specify the input file type?")
        return df
    elif args.type == "MQ":
      dfa = convert_mq(dfa, use_modified_sequence=args.use_modified_sequence, filter_retention_length=args.filter_retention_length, filter_pep=args.filter_pep, filter_contaminants=args.filter_contaminants, filter_decoys=args.filter_decoys, exclusion_list=args.exclusion_list)
    elif args.type == "PD":
      logging.error("ProteomeDiscoverer not supported yet.")
    else:
      logging.error("""
        Input type not recognized. Acceptable types are:\n
        'MQ' - MaxQuant\n
        'PD' - ProteomeDiscoverer\n
        """)

    # keep track of where observations came from
    dfa["input_id"] = i

    df = df.append(dfa)

  # remove experiments from blacklist
  if args.remove_exps is not None and len(args.remove_exps) > 0:
    exclude_exps = list(filter(lambda x: re.search(r""+args.remove_exps+"", x), df["raw_file"].unique()))
    logging.info("Removing " + str(np.sum(df["raw_file"].isin(exclude_exps))) + " observations matching " + args.remove_exps + " ...")

    keep_exps = ~(df["raw_file"].isin(exclude_exps).values)
    df = df[keep_exps]  
    df_original = df_original[keep_exps] 

    #df["exclude"] = (df["exclude"] | df["raw_file"].isin(exclude_exps))
  
  df = df.reset_index(drop=True)

  # map peptide and experiment IDs
  # sort experiment IDs alphabetically - or else the order is by 
  # first occurrence of an observation of that raw file
  df["exp_id"] = df["raw_file"].map({ind: val for val, ind in enumerate(np.sort(df["raw_file"].unique()))})
  df["peptide_id"] = df["sequence"].map({ind: val for val, ind in enumerate(df["sequence"].unique())})

  # filter for occurence in number of experiments
  if args.filter_num_exps is not None and args.filter_num_exps >= 2:
    logging.info("Removing peptides with less than " + str(args.filter_num_exps) + " occurrences in different experiments...")
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
  parser.add_argument("-fc", "--filter-contaminants", action="store_true", default=True, help="Filter contaminant proteins, as marked so by the search engine. Default: True")
  parser.add_argument("-fd", "--filter-decoys", action="store_true", default=True, help="Filter decoy matches. Default: True")
  parser.add_argument("-fr", "--filter-retention-length", type=int, default=5, help="Filter PSMs by retention time length (Peak width, in minutes). Default: 5 mins")
  parser.add_argument("-fp", "--filter-pep", type=float, default=0.5, help="Filter PSMs by PEP (Posterior Error Probability). Default: 0.5")
  parser.add_argument("-fe", "--filter-num-exps", type=int, default=2, help="Filter by occurence of peptide in number of experiments. i.e., if a peptide is not observed in at least this number of experiments over the entire set, then it will be filtered out before alignment.")
  parser.add_argument("-e", "--exclusion-list", type=argparse.FileType("r"), help="Path to exclusion list file - UniProt IDs separated by lines. Excluding contaminants that are missed by search engines can be crucial for alignment. If STAN is exceeding a reasonable iteration limit, or if the alignment residuals are unexpectedly large, consider filtering out proteins with peptides with inconsistent RTs.")
  parser.add_argument("-re", "--remove-exps", type=str, default=None, help="Regular expression of experiments (raw files) to remove from the data. Regular expression is case-sensitive. Default: None")
  parser.add_argument("-m", "--use-modified-sequence", action="store_true", default=True, help="Use modified sequence instead of unmodified peptide sequence. Default: True")

if __name__ == "__main__":
  print("nice")
  parser = argparse.ArgumentParser()
  parser.add_argument("input", type=argparse.FileType("r"), nargs="+",
    help="Input file(s) from search engine")
  parser.add_argument("-o", "--output",
    help="Path to converted file. Default: prints to stdout")

  add_converter_args(parser)

  args = parser.parse_args()

  
  if args.verbose:
    verbose = True

  if args.type == None:
    logging.error("No input file type specified. Is the input already in the converted format?")

  # set up logger
  for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler) 
   
  logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
  rootLogger = logging.getLogger()
  rootLogger.setLevel(logging.DEBUG)

  consoleHandler = logging.StreamHandler(sys.stdout)
  consoleHandler.setFormatter(logFormatter)
  rootLogger.addHandler(consoleHandler)

  logging.info(" ".join(sys.argv[1:]))

  df, df_original = process_files(args)

  logging.info(str(df.shape[0] - df["exclude"].sum()) + " / " + str(df.shape[0]) + " (" + "{:.2f}".format((df.shape[0] - df["exclude"].sum()) / df.shape[0] * 100) + "%) " + "observations pass criteria and will be used for alignment")
  
  logging.info("Saving converted data to " + args.output + " ...")
  df.to_csv(args.output, sep="\t", index=False)
