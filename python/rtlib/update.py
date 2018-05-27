#!/usr/bin/env python3
# coding: utf-8

import argparse
import logging
import numpy as np
import os
import pandas as pd
import sys
import time

from rtlib.align import add_alignment_args, align
from rtlib.converter import add_converter_args, process_files
#from rtlib.helper import load_params_from_file, add_version_arg, add_config_file_arg
from rtlib.helper import *
from scipy.stats import norm, lognorm, laplace

logger = logging.getLogger()

def update(dfa, params):
  dff = dfa[-(dfa["exclude"])]
  dff = dff.reset_index(drop=True)

  logger.info("{} / {} ({:.2%}) confident, alignable observations (PSMs) after filtering.".format(dff.shape[0], dfa.shape[0], dff.shape[0] / dfa.shape[0]))

  # refactorize peptide id into stan_peptide_id, 
  # to preserve continuity when feeding data into STAN
  dff["stan_peptide_id"] = dff["sequence"].map({ind: val for val, ind in enumerate(dff["sequence"].unique())})

  num_experiments = dff["exp_id"].max() + 1
  num_observations = dff.shape[0]
  num_peptides = dff["peptide_id"].max() + 1
  exp_names = dff["raw_file"].unique()
  mean_log_rt = np.mean(np.log(dff["retention_time"]))
  sd_log_rt = np.std(np.log(dff["retention_time"]))
  max_rt = dff["retention_time"].max()
  pep_id_list = dff["peptide_id"].unique()

  # output table
  df_new = pd.DataFrame()

  for i, e in enumerate(np.sort(dff["exp_id"].unique())):
    exp_name = exp_names[i]
    logger.info("Experiment #{} - {}".format(i, exp_name))
    
    exp = dfa[dfa["exp_id"]==e]
    exp = exp.reset_index(drop=True)
    
    # not all peptides in this experiment have data from the model
    # we can only update those that have that data. others will not be touched
    exp_matches = np.isin(exp["peptide_id"].values, pep_id_list)
    exp_f = exp[exp_matches]
    exp_f = exp_f.reset_index(drop=True)

    # convert peptide_id to stan_peptide_id
    exp_f["stan_peptide_id"] = exp_f["peptide_id"].map({ind: val for val, ind in enumerate(pep_id_list)})
    exp_peptides = exp_f["stan_peptide_id"].unique()
    exp_f["mu"] = params["peptide"]["mu"].values[exp_f["stan_peptide_id"]]

    # get mu from muij, using the linear regression parameters from this experiment
    exp_f["muij"] = 0
    # if the mu is before the split point, only account for the first segment
    exp_f["muij"][exp_f["mu"] < params["exp"]["split_point"][i]] = params["exp"]["beta_0"][i] + (params["exp"]["beta_1"][i] * exp_f["mu"])
    # if the mu is after the split point, account for both segments
    exp_f["muij"][exp_f["mu"] >= params["exp"]["split_point"][i]] = params["exp"]["beta_0"][i] + (params["exp"]["beta_1"][i] * params["exp"]["split_point"][i]) + (params["exp"]["beta_2"][i] * (exp_f["mu"] - params["exp"]["split_point"][i]))

    # get sigmaij from the sigma_intercept and sigma_slope parameters for this experiment
    exp_f["sigmaij"] = params["exp"]["sigma_intercept"][i] + params["exp"]["sigma_slope"][i] / 100 * exp_f["mu"]


    # PEP.new = P(-|RT) = P(RT|-)*P(-) / (P(RT|-)*P(-) + P(RT|+)*P(+)
    # + <- PSM=Correct
    # - <- PSM=Incorrect
    
    # P(RT|-) = probability of peptides RT, given that PSM is incorrect
    #           calculated from the uniform density from 0 to max(RT)
    #exp.rt.minus <- 1 / max(exp.f$`Retention time`) #experiment-specific
    
    # Fit3d, normal density over all retention times
    rt_mean = np.mean(exp_f["retention_time"])
    rt_std = np.std(exp_f["retention_time"])
    exp_rt_minus = norm.pdf(exp_f["retention_time"], loc=rt_mean, scale=rt_std)

    # P(-) = probability that PSM is incorrect (PEP)
    # P(+) = probability that PSM is correct (1-PEP)
    
    # P(RT|+) = probability that given the correct ID, the RT falls in the
    #           normal distribution of RTs for that peptide, for that experiment
    #
    # this is defined in fit_RT3.stan as a mixture between 2 normal distributions
    # where one distribution for the peptide RT is weighted by 1-PEP
    # and the other distribution for all RTs is weighted by PEP
    # -- summing to a total density of 1
    
    # ensure that pep does not exceed 1
    # will result in incorrect negative densities when applying mixture model
    exp_f["pep"][exp_f["pep"] > 1] = 1

    # Fit3d - mixture between two normal densities
    comp1 = exp_f["pep"] * norm.pdf(exp_f["retention_time"], loc=rt_mean, scale=rt_std)
    comp2 = (1 - exp_f["pep"]) * norm.pdf(exp_f["retention_time"], loc=exp_f["muij"], scale=exp_f["sigmaij"])
    exp_rt_plus = comp1 + comp2

    # now we can update the PEP
    # PEP.new = P(-|RT) = P(RT|-)*P(-) / (P(RT|-)*P(-) + P(RT|+)*P(+)
    # + | PSM = Correct
    # - | PSM = Incorrect
    pep_new = (exp_rt_minus * exp_f["pep"]) / ((exp_rt_minus * exp_f["pep"]) + (exp_rt_plus * (1 - exp_f["pep"])))

    # for PSMs for which we have alignment/update data
    exp_new = pd.DataFrame({
        "rt_minus":          exp_rt_minus.tolist(),
        "rt_plus":           exp_rt_plus.tolist(),
        "mu":                exp_f["mu"].values.tolist(),
        "muij":              exp_f["muij"].values.tolist(),
        "sigmaij":           exp_f["sigmaij"].values.tolist(),
        "pep_new":           pep_new.tolist(),

        "id":                exp_f["id"],
        "exp_id":            exp_f["exp_id"],
        "peptide_id":        exp_f["peptide_id"],
        "stan_peptide_id":   exp_f["stan_peptide_id"],
        "input_id":          exp_f["input_id"]
    })
    # for PSMs without alignment/update data
    exp_new = exp_new.append(pd.DataFrame({
        "rt_minus":          np.nan,
        "rt_plus":           np.nan,
        "mu":                np.nan,
        "muij":              np.nan,
        "sigmaij":           np.nan,
        "pep_new":           np.nan,

        "id":                exp["id"][~(exp_matches)],
        "exp_id":            exp["exp_id"][~(exp_matches)],
        "peptide_id":        exp["peptide_id"][~(exp_matches)],
        "stan_peptide_id":   np.nan,
        "input_id":          exp["input_id"][~(exp_matches)]
    }))
    # append to master DataFrame and continue
    df_new = df_new.append(exp_new)

  # reorder by ID and reset the index
  df_new = df_new.sort_values("id")
  df_new = df_new.reset_index(drop=True)

  return df_new

def write_output(df, out_path, args):
  # remove diagnostic columns, unless they are specified to be kept
  if not args.add_diagnostic_cols:
    df = df.drop(["input_exclude", "exclude", "mu", "muij", "rt_minus", "rt_plus", "sigmaij", "input_id", "exp_id", "peptide_id", "stan_peptide_id"], axis=1)

  df.to_csv(out_path, sep="\t", index=False)

def add_update_args(parser):
  parser.add_argument("-p", "--params-folder", type=str, default=None, help="Path to folder with params files from alignment")
  parser.add_argument("-c", "--combine-output", action="store_true", default=False, help="If multiple input files are provided, then combine them into one output file. This may throw errors if input files have different column formats. Default: False")
  parser.add_argument("--combined-name", type=str, default="ev_combined.txt", help="If combining inputs into one output file with --combine-output, then this is the name of that file. Default: ev_combined.txt")
  parser.add_argument("--output-suffix", type=str, default="_c", help="Suffix added to output file(s). Default: '_c'.")
  parser.add_argument("--add-diagnostic-cols", action="store_true", default=False, help="Add diagnostic columns to the output data frame. These include the canonical retention time, adjusted retention time, standard deviation of the RT for that peptide sequence, and other variables used in the bayesian update. Default: False")
  #parser.add_argument("--save-converted-input")

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("input", type=argparse.FileType("r"), nargs="+", help="Input file(s) from search engine")
  parser.add_argument("-o", "--output", type=str, default="./rt_update", help="Path to output data. Default: './rt_update'")
  
  add_converter_args(parser)
  add_alignment_args(parser)
  add_update_args(parser)
  add_global_args(parser)

  args = parser.parse_args()
  #args = get_args(parser)
  #print(vars(args))

  # create output folder
  if args.output is None or len(args.output) < 1:
    args.output = "./alignment"
  if not os.path.exists(args.output):
    os.makedirs(args.output)

  # set up logger
  logger = init_logger(args.verbose, os.path.join(args.output, "update.log"))

  df, df_original = process_files(args)

  logger.info("Finished converting files and filtering PSMs")

  # load params, either from a defined folder or from running the alignment
  params = {}
  if args.params_folder is None:
    logger.info("Running alignment...")
    params = align(df, filter_pep=args.filter_pep, mu_min=args.mu_min, rt_distortion=args.rt_distortion, prior_iters=args.prior_iters, stan_iters=args.stan_iters, stan_file=args.stan_file, save_params=args.save_params, print_figures=args.print_figures, output_path=args.output, verbose=args.verbose)
  else:
    params = load_params_from_file(args.params_folder)

  # now we have the params, run the update
  logger.info("Updating PEPs with alignment data...")
  df_new = update(df, params)

  # save the sparse combined input file?
  #df_new.to_csv(os.path.join(args.output, "df_converted.txt"), sep="\t", index=False)

  # add new columns to original DF, and remove the duplicate ID column
  logger.info("Concatenating results to original data...")

  df_adjusted = pd.concat([df_original.loc[~df_original["input_exclude"]].reset_index(drop=True), df_new.drop(["id", "input_id"], axis=1).reset_index(drop=True)], axis=1)

  # add rows of removed experiments (done with the --remove-exps options)
  if args.remove_exps is not None:
    # store a copy of the columns and their order for later
    df_cols = df_adjusted.columns
    # concatenate data frames
    df_adjusted = pd.concat([df_adjusted, df_original.loc[df_original["input_exclude"]]], axis=0, ignore_index=True)
    # pd.concat reindexes the order of the columns, 
    # so just order it back to what it used to be
    df_adjusted = df_adjusted.reindex(df_cols, axis=1)

  # sort by ID, and reset index
  df_adjusted = df_adjusted.sort_values(["id"])
  df_adjusted = df_adjusted.reset_index(drop=True)

  # add pep_updated column - which is pep_new, with the NaNs filled in
  # with the old PEPs.
  # this PEP columns will change between search engine output types, so make sure to
  # update this as well when the time comes...
  pep_col = "PEP"
  df_adjusted["pep_updated"] = df_adjusted["pep_new"]
  df_adjusted["pep_updated"][pd.isnull(df_adjusted["pep_new"])] = df_adjusted[pep_col][pd.isnull(df_adjusted["pep_new"])]

  # write to file
  if args.combine_output:
    out_path = os.path.join(args.output, args.combined_name)
    logger.info("Combining input file(s) and writing adjusted data file to {} ...".format(out_path))
    write_output(df_adjusted, out_path, args)
  else:
    logger.info("Saving output to separate files...")
    for i, f in enumerate(args.input):
      out_path = os.path.join(args.output, os.path.splitext(os.path.basename(f.name))[0] + args.output_suffix + str(i) + ".txt")
      logger.info("Saving input file {} to {}".format(i, out_path))
      df_a = df_adjusted.loc[df_adjusted["input_id"] == i]
      # save to file
      # other data formats might have a different separator, or have an index column
      write_output(df_a, out_path, args)

  logger.info("Done!")

if __name__ == "__main__":
  main()




