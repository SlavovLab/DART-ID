#!/usr/bin/env python3
# coding: utf-8

import argparse
import logging
import numpy as np
import os
import pandas as pd
import re

from rtlib.helper import *

logger = logging.getLogger("root")

def main():
  # load command-line args
  parser = argparse.ArgumentParser()  
  add_global_args(parser, add_config_file=False)

  parser.add_argument("--include", type=str, default=None, help="Inclusion expression (regular expression) that is evaluated against all raw file names to determine whether or not that raw file is included.")
  parser.add_argument("--raw-file-col", type=str, default="Raw file", help="Name of the raw file column.")
  args = parser.parse_args()

  # initialize logger
  init_logger(args.verbose, "", log_to_file=False)
  logger = logging.getLogger("root")

  df = pd.DataFrame()

  # iterate through each input file provided.
  for i, f in enumerate(args.input):

    logger.info("Reading in input file #{} | {} ...".format(i+1, f.name))

    # load the input file with pandas
    # 
    # have a variable low memory option depending on the input type.
    # MaxQuant, for example, has a structure that forces pandas out of its
    # optimal low memory mode, and we have to specify it here.
    dfa = pd.read_csv(f.name, sep="\t", low_memory=False)

    # keep track of where observations came from. this is _not_ the raw file ID
    # but instead the ID from which input file it originated from, so that if
    # we need to split these observations up by input file in the future we can do so
    dfa["input_id"] = i

    # append a copy of dfa into df, because the conversion process will heavily
    # modify dfa. we need to keep a copy of the original dataframe in order to append
    # the new columns back onto it later.
    # re-index columns with "[dfa.columns.tolist()]" to preserve the general column order
    df = df.append(dfa)[dfa.columns.tolist()]

  logger.info("Raw files:")
  logger.info(df[args.raw_file_col].unique())

  if args.include is not None:
    if len(args.include) == 0:
      raise Exception("Experiment inclusion expression by raw file name provided, but expression is defined incorrectly.")
    else:
      # see if any raw file names match the user-provided expression
      include_exps = list(filter(lambda x: re.search(r"" + args.include + "", x), df[args.raw_file_col].unique()))

      logger.info("{} raw files match the inclusion expression \"{}\"".format(len(include_exps), args.include))
      logger.info(include_exps)
      logger.info("Retaining {} observations out of {} total".format(np.sum(df[args.raw_file_col].isin(include_exps)), df.shape[0]))

      include_exps = df[args.raw_file_col].isin(include_exps).values
      df = df.loc[include_exps].reset_index(drop=True)

  # if combining input files, then write to one combined file
  logger.info("Combining input file(s) and writing adjusted data file to {} ...".format(args.output))
  df.to_csv(args.output, sep="\t", index=False)



def mq2pin():
  # load command-line args
  parser = argparse.ArgumentParser()  
  add_global_args(parser)
  args = parser.parse_args()

  # load config file
  # this function also creates the output folder
  config = read_config_file(args, create_output_folder=False)

  # initialize logger
  init_logger(config["verbose"], "", log_to_file=False)
  logger = logging.getLogger("root")

  df = pd.DataFrame()

  # iterate through each input file provided.
  for i, f in enumerate(config["input"]):
    # first expand user or any vars
    f = os.path.expanduser(f)
    f = os.path.expandvars(f)

    logger.info("Reading in input file #{} | {} ...".format(i+1, f))

    dfa = pd.read_csv(f, sep="\t", low_memory=False)

    # take only the columns we need
    cols = config["mq2pin_cols"]

    for key in list(cols.keys()):
      if cols[key] is None: cols.pop(key, None)

    dfa = dfa[list(cols.values())]
    dfa.columns = list(cols.keys())

    df = df.append(dfa)

  # unique ID that we'll use for SpecId
  df["id"] = range(0, df.shape[0])

  # remove rows that weren't MS/MSed
  df = df[~pd.isnull(df["mass_calibration"])].reset_index(drop=True)

  pin = pd.DataFrame()

  pin["SpecId"] = df["id"]

  # target or decoy
  pin["Label"] = df["protein"].str.contains("REV__").values.astype(int)
  pin["Label"].loc[pin["Label"]==1] = -1
  pin["Label"].loc[pin["Label"]==0] = 1

  # adjust scan numbers so that the scan numbers from different
  # experiments don't overlap
  max_scan_nums = df.groupby("raw_file")["scan_num"].apply((lambda x: np.max(x)))
  max_scan_nums = max_scan_nums[np.argsort(max_scan_nums.index.values)]
  max_scan_nums = np.cumsum(max_scan_nums) - max_scan_nums[0]

  pin["ScanNr"] = df["scan_num"] + df["raw_file"].map(max_scan_nums)

  # ExpMass - measured mass of precursor ion (measured m/z * charge)
  pin["ExpMass"] = df["mz"] * df["charge"]
  # CalcMass - theoretical mass of precursor ion
  pin["CalcMass"] = df["mass"]

  # --doc features, RT and mass calibration
  pin["RT"] = df["retention_time"]
  pin["MassCalibration"] = df["mass_calibration"]

  # precursor properties
  pin["Mass"] = df["mass"]

  # Score - most important feature
  pin["Score"] = df["score"]
  pin["DeltaScore"] = df["delta_score"]

  # sequence characteristics
  # peptide length
  pin["PepLen"] = df["length"]
  # charge states
  pin["Charge1"] = (df["charge"] == 1).values.astype(int)
  pin["Charge2"] = (df["charge"] == 2).values.astype(int)
  pin["Charge3"] = (df["charge"] == 3).values.astype(int)

  # enzymatic performance - trypsin
  # don't have data on n-terminus cleavage...
  pin["enzC"] = df["sequence"].str.slice(-1).isin(["R", "K"]).values.astype(int)
  # missed cleavages = number of enzymatic sites 
  pin["enzInt"] = df["missed_cleavages"]

  # lastly, the peptide and protein
  pin["Peptide"] = df["sequence"]
  # need to surround peptide with flanking amino acids
  # we don't have this data, so just append alanines on both sides
  pin["Peptide"] = "A." + pin["Peptide"] + ".A"

  pin["Protein"] = df["protein"]
  # extract UniProt IDs from protein string
  pin["Protein"] = pin["Protein"].str.split("|").apply((lambda x: x[1] if len(x) is 3 else x[0]))
  # mark reverse ones differently still
  pin["Protein"].loc[pin["Label"]==-1] = ("REV_" + pin["Protein"].loc[pin["Label"]==-1])

  
  # write headers and weights
  with open(config["output"], "w") as f:

    # headers
    for i, col in enumerate(pin.columns):
      f.write(col)
      if i != (len(pin.columns)-1):
        f.write("\t")
    f.write("\n")

    # weights
    weights = [
    # empty for specid, label, scannr, expmass, calcmass, 
    # and doc features (RT and mass calibration) - 7 blanks total
    "DefaultDirection", "-", "-", "-", "-", "-", "-",
    "0", # mass
    "1", # score
    "1.5", # delta score
    "-0.573", # peptide length
    "0.0335", "0.149", "-0.156", # charge states
    "0", "0" # enzymatic features
    # skip peptide and protein
    ]
    for i, w in enumerate(weights):
      f.write(w)
      if i != (len(weights)-1):
        f.write("\t")
    f.write("\n")

    f.close()

  logger.info("Writing percolator-in (pin) file to {} ...".format(config["output"]))
  pin.to_csv(config["output"], sep="\t", header=False, index=False, mode="a")
  

if __name__ == "__main__":
  main()