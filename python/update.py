#!/usr/bin/env python3
# coding: utf-8

import argparse
import logging
import numpy as np
import os
import pandas as pd
import sys
import time

from align import add_alignment_args, align
from converter import add_converter_args, convert_pd, convert_mq, process_files
from scipy.stats import norm, lognorm, laplace

verbose = False

def update(dfa, params):
  # create a unique ID for each PSM to help with stiching the final result together
  # after all of our operations
  dfa["id"] = range(0, dfa.shape[0])

  dff = dfa[-(dfa["exclude"])]
  dff = dff.reset_index(drop=True)

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
    logging.info("Experiment #" + str(i) + " - " + exp_name)
    
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
    
    def mu_to_muij(mu):
        if mu < params["exp"]["split_point"][i]:
            return params["exp"]["beta_0"][i] + (params["exp"]["beta_1"][i] * mu)
        else:
            return params["exp"]["beta_0"][i] + (params["exp"]["beta_1"][i] * params["exp"]["split_point"][i]) + (params["exp"]["beta_2"][i] * (mu - params["exp"]["split_point"][i]))
    
    exp_f["muij"] = exp_f["mu"].apply(mu_to_muij)
    exp_f["sigmaij"] = params["exp"]["sigma_intercept"][i] + params["exp"]["sigma_slope"][i] / 100 * exp_f["mu"]

    # PEP.new = P(-|RT) = P(RT|-)*P(-) / (P(RT|-)*P(-) + P(RT|+)*P(+)
    # + <- PSM=Correct
    # - <- PSM=Incorrect
    
    # P(RT|-) = probability of peptides RT, given that PSM is incorrect
    #           calculated from the uniform density from 0 to max(RT)
    #exp.rt.minus <- 1 / max(exp.f$`Retention time`) #experiment-specific
    
    # Fit3d, normal density over all retention times
    exp_rt_minus = norm.pdf(exp_f["retention_time"], loc=np.mean(exp_f["retention_time"]), scale=np.std(exp_f["retention_time"]))

    # P(-) = probability that PSM is incorrect (PEP)
    # P(+) = probability that PSM is correct (1-PEP)
    
    # P(RT|+) = probability that given the correct ID, the RT falls in the
    #           normal distribution of RTs for that peptide, for that experiment
    #
    # this is defined in fit_RT3.stan as a mixture between 2 normal distributions
    # where one distribution for the peptide RT is weighted by 1-PEP
    # and the other distribution for all RTs is weighted by PEP
    # -- summing to a total density of 1
    
    # define a function to apply to each unique peptide
    # x is a DataFrame where all stan_peptide_ids are shared
    def rt_plus(x):
      # ensure that pep does not exceed 1
      # will result in incorrect negative densities when applying mixture model
      pep = x["pep"]
      pep[pep > 1] = 1
      # Fit3d - mixture between two normal densities
      comp1 = (pep) * norm.pdf(x["retention_time"], loc=np.mean(exp_f["retention_time"]), scale=np.std(exp_f["retention_time"]))
      comp2 = (1 - pep) * norm.pdf(x["retention_time"], loc=x["muij"], scale=x["sigmaij"])
      y = comp1 + comp2
      return y.values.tolist()

    # apply the above function to each unique peptide
    exp_rt_plus = exp_f.groupby("stan_peptide_id")[["retention_time", "muij", "sigmaij", "pep"]].apply(rt_plus).values.tolist()
    # collapse the list of lists into just one list
    exp_rt_plus = np.array([st for row in exp_rt_plus for st in row])
    
    # sometimes rt.plus will go so low that it will round to 0
    # just round this back up to the smallest number R will handle
    exp_rt_plus[exp_rt_plus == 0] = np.finfo(float).eps

    exp_PEP = exp_f["pep"].values
    # sometimes MQ will output PEP > 1, which makes no sense, and will
    # result in negative values for our adjusted PEP
    # set all PEP > 1 to PEP = 1
    exp_PEP[exp_PEP > 1] = 1

    # now we can update the PEP
    # PEP.new = P(-|RT) = P(RT|-)*P(-) / (P(RT|-)*P(-) + P(RT|+)*P(+)
    # + | PSM = Correct
    # - | PSM = Incorrect
    pep_new = (exp_rt_minus * exp_PEP) / ((exp_rt_minus * exp_PEP) + (exp_rt_plus * (1 - exp_PEP)))
    
    # for PSMs for which we have alignment/update data
    exp_new = pd.DataFrame({
        "rt_minus": exp_rt_minus.tolist(),
        "rt_plus": exp_rt_plus.tolist(),
        "mu": exp_f["mu"].values.tolist(),
        "muij": exp_f["muij"].values.tolist(),
        "sigmaij": exp_f["sigmaij"].values.tolist(),
        "pep_new": pep_new.tolist(),
        "id": exp_f["id"]
    })
    # for PSMs without alignment/update data
    exp_new = exp_new.append(pd.DataFrame({
        "rt_minus": np.nan,
        "rt_plus": np.nan,
        "mu": np.nan,
        "muij": np.nan,
        "sigmaij": np.nan,
        "pep_new": np.nan,
        "id": exp["id"][~(exp_matches)]
    }))
    # append to master DataFrame and continue
    df_new = df_new.append(exp_new)

  # reorder by ID and reset the index
  df_new = df_new.sort_values("id")
  df_new = df_new.reset_index(drop=True)

  return df_new


def add_update_args(parser):
  parser.add_argument("-p", "--params-folder", type=str, default=None,
    help="Path to folder with params files from alignment")
  parser.add_argument("-c", "--combine-inputs", action="store_true", default=False, help="If multiple input files are provided, then combine them into one output file. This may throw errors if input files have different column formats. Default: False")

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  #parser.add_argument("input", help="Input file from search engine", type=str)
  parser.add_argument("input", type=argparse.FileType("r"), nargs="+",
    help="Input file(s) from search engine")
  parser.add_argument("-o", "--output", type=str, default="./rt_update",
    help="Path to output data. Default: './rt_update'")
  
  
  add_converter_args(parser)
  add_alignment_args(parser)
  add_update_args(parser)

  args = parser.parse_args()

  if args.verbose:
    verbose = True

  # create output folder
  if args.output is None or len(args.output) < 1:
    args.output = "./alignment"
  if not os.path.exists(args.output):
    os.makedirs(args.output)

  # set up logger
  for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler) 
   
  logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
  rootLogger = logging.getLogger()

  fileHandler = logging.FileHandler(os.path.join(args.output, "update.log"), mode="w")
  fileHandler.setFormatter(logFormatter)
  rootLogger.addHandler(fileHandler)

  consoleHandler = logging.StreamHandler()
  consoleHandler.setFormatter(logFormatter)
  rootLogger.addHandler(consoleHandler)

  logging.info(" ".join(sys.argv[1:]))

  df, df_original = process_files(args)

  logging.info("Finished converting files")

  # load params, either from a defined folder or from running the alignment
  params = {}
  if args.params_folder is None:
    logging.info("Running alignment...")
    params = align(df, filter_pep=args.filter_pep, mu_min=args.mu_min, rt_distortion=args.rt_distortion, prior_iters=args.prior_iters, stan_iters=args.stan_iters, stan_file=args.stan_file, print_figures=args.print_figures, output_path=args.output)
  else:
    # load parameters if they are specified in the command line
    logging.info("Using provided alignment parameters. Loading params...")
    param_files = ["exp_params.txt", "pair_params.txt", "peptide_params.txt"]
    for pf in param_files:
      pfp = os.path.join(args.params_folder, pf)
      if os.path.exists(pfp):
        with open(pfp, "rb") as f:
          try:
            params[pf.split("_")[0]] = pd.read_csv(pfp, sep="\t")
          except:
            logging.error("Error loading param file")

  # now we have the params, run the update
  logging.info("Updating PEPs with alignment data...")
  df_new = update(df, params)

  df_new.to_csv(os.path.join(args.output, "ev_new.txt"), sep="\t", index=False)

  # add new columns to original DF, and remove the duplicate ID column
  logging.info("Concatenating results to original data...")
  df_adjusted = pd.concat([df_original.reset_index(drop=True), df_new.drop(["id"], axis=1).reset_index(drop=True)], axis=1)

  # write to file
  logging.info("Writing adjusted data file...")
  df_adjusted.to_csv(os.path.join(args.output, "ev_adjusted.txt"), sep="\t", index=False)

  logging.info("Done!")




