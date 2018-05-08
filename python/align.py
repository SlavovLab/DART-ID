#!/usr/bin/env python3
# coding: utf-8

import argparse
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
import pystan
import re
import sys
import time

from converter import add_converter_args, convert_pd, convert_mq, process_files
from hashlib import md5
from scipy.stats import laplace
from sklearn import linear_model

pd.options.mode.chained_assignment = None
verbose = False

def align(dfa, filter_pep=0.5, mu_min=1, rt_distortion=10, prior_iters=10, stan_iters=2e4, stan_attempts=3, stan_file=None, print_figures=False, output_path=None):

  # take subset of confident observations to use for alignment
  dff = dfa[-(dfa["exclude"])]
  dff = dff.reset_index(drop=True)

  logging.info(str(dff.shape[0]) + " / " + str(dfa.shape[0]) + " (" + "{:.2f}".format(dff.shape[0]/dfa.shape[0]*100) + "%)" + " confident, alignable observations (PSMs) after filtering.")

  # refactorize peptide id into stan_peptide_id, 
  # to preserve continuity when feeding data into STAN
  dff["stan_peptide_id"] = dff["sequence"].map({ind: val for val, ind in enumerate(dff["sequence"].unique())})

  exp_names = np.sort(dff["raw_file"].unique())
  num_experiments = len(exp_names)
  num_observations = dff.shape[0]
  num_peptides = dff["stan_peptide_id"].max() + 1
  
  logging.info('Building peptide-experiment pairs...')

  # build a unique peptide-experiment ID
  pep_exp_all = dff["stan_peptide_id"].map(str) + " - " + dff["exp_id"].map(str)
  pep_exp_pairs = pep_exp_all.unique()
  num_pep_exp_pairs = len(pep_exp_pairs)

  # muij_map - maps pair ID back to dataframe index
  muij_map = pep_exp_all.map({ind: val for val, ind in enumerate(pep_exp_pairs)})
  # maps to experiment and peptide ID
  splt = pd.Series(pep_exp_pairs).str.split(" - ")
  muij_to_pep = splt.str.get(0).map(int)
  muij_to_exp = splt.str.get(1).map(int)

  logging.info(str(num_pep_exp_pairs) + " peptide-experiment pairs.")

  # build dictionary of data to feed to STAN
  # revert all data to primitive types to avoid problems later
  # STAN code is all 1-indexed, so add 1 to all indexed forms of data
  stan_data = {
    "num_experiments": num_experiments,
    "num_peptides": num_peptides,
    "num_pep_exp_pairs": num_pep_exp_pairs,
    "num_total_observations": num_observations,
    "muij_map": muij_map+1,
    "muij_to_pep": (muij_to_pep+1).tolist(),
    "muij_to_exp": (muij_to_exp+1).tolist(),
    "experiment_id": (dff["exp_id"]+1).tolist(),
    "peptide_id": (dff["stan_peptide_id"]+1).tolist(),
    "retention_times": dff["retention_time"].tolist(),
    "mean_log_rt": np.mean(np.log(dff["retention_time"])),
    "sd_log_rt": np.std(np.log(dff["retention_time"])),
    "mean_rt": np.mean(dff["retention_time"]),
    "sd_rt": np.std(dff["retention_time"]),
    "pep": dff["pep"].tolist(),
    "max_retention_time": dff["retention_time"].max()
  }

  logging.info("Initializing fit priors for " + str(num_peptides) + " peptides...")

  # get the average retention time for a peptide, weighting by PEP
  def get_mu(x):
      weights = ((1 - x["pep"].values) - (1 - filter_pep)) / filter_pep
      return np.sum(x["retention_time"].values * weights) / np.sum(weights)
  
  # apply the get_mu function on all peptides and add some distortion
  mu_init = dff.groupby("stan_peptide_id")[["pep", "retention_time"]].apply(get_mu).values + np.random.normal(0, rt_distortion, num_peptides)

  # negative or very low retention times not allowed. floor at 5 minutes
  mu_init[mu_init <= mu_min] = mu_min
  # canonical retention time shouldn't be bigger than largest real RT
  mu_max = dff["retention_time"].max()
  mu_init[mu_init > mu_max] = mu_max

  # take retention times and distort
  rt_distorted = dff["retention_time"] + np.random.normal(0, rt_distortion, len(dff["retention_time"]))
  # make sure distorted retention times stay within bounds of real ones
  rt_distorted[rt_distorted > dff["retention_time"].max()] = dff["retention_time"].max()
  rt_distorted[rt_distorted < dff["retention_time"].min()] = dff["retention_time"].min()

  # initialize priors for the segmented linear regression
  # first element of vector is beta_0, or the intercept
  # second element is beta_1 and beta_2, the slopes of the two segments
  beta_init = np.array((np.repeat(10, num_experiments), np.repeat(1, num_experiments)), dtype=float)

  logging.info("Optimizing priors with linear approximation for " + str(prior_iters) + " iterations.")

  mu_pred = np.zeros(num_peptides)
  # temporary data frame to quickly map over in the loop
  dft = pd.DataFrame(dict(stan_peptide_id=dff["stan_peptide_id"], exp_id=dff["exp_id"], pep=dff["pep"], retention_time=mu_init[dff["stan_peptide_id"]]))

  for i in range(0, prior_iters):
    # for each experiment, fit a simple linear regression
    # between the distorted RTs and the initial canonical retention times
    for j in np.sort(dff["exp_id"].unique()):
        idx     = (dff["exp_id"] == j)
        rt_cur  = rt_distorted[idx]
        mu_cur  = mu_init[dff["stan_peptide_id"][idx]]
        pep_cur = dff["pep"][idx]

        # for this experiment, run a linear regression (1 degree polyfit)
        # of the mus and the distorted RTs. store the linear regression params
        m, c = np.polyfit(mu_cur, rt_cur, 1, w=(1 - pep_cur))
        beta_init[(0,1), j] = [c, m]

    # calculate new set of canonical RTs based on linear regression params
    mu_pred = (rt_distorted - beta_init[0][dff["exp_id"]]) / beta_init[1][dff["exp_id"]] 
    # make sure new canonical RTs are within same range as distorted RTs
    mu_pred[mu_pred <= 0] = rt_distorted.min()
    mu_pred[mu_pred >= rt_distorted.max()] = rt_distorted.max()
    dft["retention_time"] = np.copy(mu_pred)

    mu_prev = np.copy(mu_init)

    # new set of priors for canonical RTs based on weighted combination of
    # this set of predicted canonical RTs
    mu_init = dft.groupby("stan_peptide_id")[["pep", "retention_time"]].apply(get_mu).values

    logging.info("Iter " + str(i + 1) + " | Avg. canonical RT shift: " + "{:.5f}".format(pow(np.sum(mu_prev - mu_init), 2) / len(mu_init)))

  # grab linear regression params
  # set beta_2 (slope of second segment) to the same as the slope of the 1st segment
  beta_0 = beta_init[0]
  beta_1 = beta_init[1]
  beta_2 = np.copy(beta_1)

  # apply lower bound of (-1.0 * min(beta_1) * min(muInit)) to beta_0
  # where (-1.0 * min(beta_1) * min(muInit)) is the lowest possible intercept
  # given the lowest possible mu and lowest possible beta_1
  beta_0[beta_0 <= (-1 * beta_1.min() * mu_init.min())] = (-1 * beta_1.min() * mu_init.min()) + 1e-3

  # apply upper bound to prior canonical RTs
  mu_init[mu_init >= dff["retention_time"].max()] = 0.95 * dff["retention_time"].max()

  # create prior list for STAN
  init_list = {
    "mu": mu_init,
    "beta_0": beta_0,
    "beta_1": beta_1,
    "beta_2": beta_2,
    "sigma_slope": np.repeat(0.1, num_experiments),
    "sigma_intercept": np.repeat(0.1, num_experiments),
    "split_point": np.repeat(np.median(mu_init), num_experiments)
  }

  # run STAN, store optimization parameters
  op = stan_optimize(stan_file=stan_file, stan_data=stan_data, init_list=init_list, stan_iters=stan_iters, stan_attempts=stan_attempts, verbose=verbose)

  # write parameters to file
  # exp_params - contains regression parameters for each experiment
  # peptide_params - contains canonical RT (mu) for each peptide sequence
  # pair_params - contains muij and sigmaij for each peptide-experiment pair
  #               not exactly necessary as these can be extrapolated from 
  #               exp_params and peptide_params
  # 
  logging.info("Writing STAN parameters to file...")
  exp_params = pd.DataFrame({ key: op[key] for key in ["beta_0", "beta_1", "beta_2", "split_point", "sigma_slope", "sigma_intercept"]})
  peptide_params = pd.DataFrame({ key: op[key] for key in ["mu"]})
  pair_params = pd.DataFrame({ key: op[key] for key in ["muij", "sigma_ij"]})

  # add exp_id to exp_params
  exp_params["exp_id"] = np.sort(dff["exp_id"].unique())

  # put the maps in the pair parameters file as well
  pair_params = pair_params.assign(muij_to_pep=muij_to_pep)
  pair_params = pair_params.assign(muij_to_exp=muij_to_exp)

  # write parameters to file, so operations can be done on alignment data without
  # the entire alignment to run again
  exp_params.to_csv(os.path.join(output_path, "exp_params.txt"), sep="\t", index=False)
  pair_params.to_csv(os.path.join(output_path, "pair_params.txt"), sep="\t", index=True, index_label="pair_id")
  peptide_params.to_csv(os.path.join(output_path, "peptide_params.txt"), sep="\t", index=True, index_label="peptide_id")

  # write parameters to dict for further operations
  params = {}
  params["exp"] = exp_params
  params["pair"] = pair_params
  params["peptide"] = peptide_params

  # generate alignment plots?
  if print_figures:
    logging.info("Generating Alignment Figures...")
    generate_figures(params, exp_names, muij_map, dff["retention_time"], dff["pep"], output_path)

  return params

def stan_optimize(stan_file, stan_data, init_list, stan_iters, stan_attempts, verbose):
  sm = StanModel_cache(model_file=stan_file)

  # sometimes STAN will error out due to bad RNG or bad priors
  # set a limit on how many times we will try this stan configuration 
  # before running it all the way back again
  num_tries = stan_attempts
  counter = 1
  op = None
  while op is None and counter <= num_tries:
    try:
      logging.info("Starting STAN Model | Attempt #" + str(counter) + " ...")
      start = time.time()
      op = sm.optimizing(data=stan_data, init=init_list, iter=stan_iters, verbose=verbose)
      logging.info("STAN Model Finished. Run time: " + "{:.3f}".format(time.time() - start) + " seconds.")
    except RuntimeError as e:
      logging.error(str(e))
      counter = counter + 1

  if op is None:
    raise Exception("Maximum number of tries exceeded for STAN. Please re-run process or choose different parameters.")

  return op

def generate_figures(params, exp_names, muij_map, retention_times, pep, output_path):
  # split PEP into 10 bins, for coloring points later
  pep_col_code = pd.cut(pep, 10)
  # make figures folder if it doesn't exist yet
  figures_path = os.path.join(output_path, "alignment_figs")
  if not os.path.exists(figures_path):
    logging.info("Path for alignment figures " + figures_path + "does not exist. Creating...")
    os.makedirs(figures_path)

  # generate figures for each experiment
  for exp in params["exp"]["exp_id"].values:
    logging.info("Generating Summary for Experiment " + str(exp) + " | " + exp_names[exp])

    exp_params = params["exp"].iloc[exp]
    exp_indices = params["pair"]["muij_to_exp"] == exp

    # predicted RTs and SDs
    muijs = params["pair"]["muij"][exp_indices].values
    sigmas = params["pair"]["sigma_ij"][exp_indices].values
    # same as above, but with duplicates for same peptide-experiment pairs
    predicted = params["pair"]["muij"][muij_map][exp_indices[muij_map]]
    predicted_sd = params["pair"]["sigma_ij"][muij_map][exp_indices[muij_map]]

    # canonical RTs
    mus = params["peptide"]["mu"][params["pair"]["muij_to_pep"][muij_map][exp_indices[muij_map]]]

    # observed values
    observed = retention_times.values[exp_indices[muij_map]]
    obs_peps = pep.values[exp_indices[muij_map]]
    obs_code = pep_col_code.values[exp_indices[muij_map]]
    residual = observed - predicted

    # plot the 2-segment linear fit of mus to observed RTs
    plt.subplot(221)
    plt.scatter(mus, observed, s=1, color="black")
    plt.plot([0, exp_params["split_point"]],
             [exp_params["beta_0"], (exp_params["split_point"] * exp_params["beta_1"]) + exp_params["beta_0"]],
            color="red")
    plt.plot([exp_params["split_point"], 300], 
             [(exp_params["split_point"] * exp_params["beta_1"]) + exp_params["beta_0"], (exp_params["split_point"] * exp_params["beta_1"]) + ((300-exp_params["split_point"]) * exp_params["beta_2"]) + exp_params["beta_0"]],
            color="green")
    plt.plot(np.repeat(exp_params["split_point"], 2), [-100, 300], color="blue", linestyle="dashed")
    plt.axis([0, mus.max() + 10, exp_params["beta_0"]-10, observed.max() + 10])
    plt.title(exp_names[exp])
    plt.xlabel("Canonical RT (min)")
    plt.ylabel("Observed RT (min)")

    # plot residual vs observed RTs
    plt.subplot(222)
    plt.scatter(predicted, observed, color="black", s=1)
    plt.plot([0, 300], [0, 300])
    plt.plot(np.repeat(exp_params["split_point"], 2), [-100, 300], color="blue", linestyle="dashed")
    plt.axis([0, predicted.max()+10, 0, observed.max()+10])
    plt.xlabel("Predicted RT (min)")
    plt.ylabel("Observed RT (min)")

    # plot residuals, quantiles of residuals, and color points by PEP
    plt.subplot(223)
    plt.scatter(predicted, residual, s=1, c=pep_col_code.cat.codes.values[exp_indices[muij_map]])
    plt.plot([0, 300], [0, 0], color="blue")
    plt.plot(np.repeat(exp_params["split_point"], 2), [-100, 300], color="blue", linestyle="dashed")
    plt.plot(predicted.values[np.argsort(predicted)], laplace.ppf(0.025, loc=0, scale=predicted_sd)[np.argsort(predicted)], color="red")
    plt.plot(predicted.values[np.argsort(predicted)], laplace.ppf(0.975, loc=0, scale=predicted_sd)[np.argsort(predicted)], color="red")
    plt.axis([predicted.min()-5, predicted.max()+5, residual.min()-5, residual.max()+5])
    cbar = plt.colorbar()
    cbar.ax.set_yticklabels(pep_col_code.cat.categories.values)
    plt.xlabel("Predicted RT (min)")
    plt.ylabel("Residual RT (min)")

    # add some space between subplots
    plt.subplots_adjust(hspace=0.3, wspace=0.3)

    # finalize and save figure
    fig = plt.gcf()
    fig.set_size_inches(7, 7)
    fname = os.path.join(figures_path, str(exp) + "_" + exp_names[exp] + ".png")
    logging.info("Saving figure to " + fname + " ...")
    fig.savefig(fname, dpi=160)
    
    plt.close()
    fig.clf()

# taken from: https://pystan.readthedocs.io/en/latest/avoiding_recompilation.html#automatically-reusing-models
def StanModel_cache(model_file=None, model_name=None):
  if model_file is None:
    model_file = "./fit_RT3d.stan"
  elif type(model_file) is not str:
    model_file = model_file.name

  if model_name is None:
    model_name = os.path.basename(model_file)
  # cleanse model name for use as a var name in the C++ code
  model_name = re.sub(r"_|-|\.|\.stan|\.STAN|\.Stan", "", model_name)

  model_code = ""
  with open(model_file, "r") as f:
    model_code = f.read()

  # Use just as you would `stan`
  code_hash = md5(model_code.encode("ascii")).hexdigest()
  cache_fn = "cached-model-{}-{}.pkl".format(model_name, code_hash)

  try:
      # load cached model from file
      sm = pickle.load(open(cache_fn, "rb"))
  except:
      # compile model
      sm = pystan.StanModel(model_name=model_name, model_code=model_code)
      with open(cache_fn, "wb") as f:
          # save model to file
          pickle.dump(sm, f)
  else:
      logging.info("Using cached StanModel: " + model_name + "_" + code_hash)

  return sm

def add_alignment_args(parser):
  parser.add_argument("-s", "--stan-file", type=argparse.FileType("r"), default=None, help="Path to STAN fit file. Leave empty to use the provided STAN configuration.")
  parser.add_argument("-mm", "--mu-min", type=float, default=1.0, help="Minimum value of canonical retention time for a peptide, when first calculating priors. Any canonical RT below this value will floor to this value.")
  parser.add_argument("-rd", "--rt-distortion", type=float, default=10.0, help="Distortion of retention times before alignment, in minutes. Increase this number if the STAN alignment is too close to optima before alignment. Default: 10. 10 dminutes is for a normal 180 min. proteomics run, so reduce this respective to your average run-time.")
  parser.add_argument("-ip", "--prior-iters", type=int, default=10, help="Number of iterations for generating priors before STAN alignment. Increase this number to get closer to the optima before alignment. Default: 10")
  parser.add_argument("-is", "--stan-iters", type=int, default=1e5, help="Number of max iterations passed onto STAN. Optimization may terminate normally before reaching this number of iterations. For larger amounts of data, this may need to be increased in order to terminate the optimization normally at the gradient threshold. Default: 20,000")
  parser.add_argument("-as", "--stan-attempts", type=int, default=3, help="Number of attempts to run STAN with. Consider changing parameters to improve generation of priors before trying to run STAN more times.")
  parser.add_argument("-f", "--print-figures", action="store_true", default=False,
    help="Print alignment figures for each experiment that is aligned. Each alignment figure is stored as an image, with a summary of the segmented fit, as well as the residuals of the fit.")


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  #parser.add_argument("input", help="Input file from search engine", type=str)
  parser.add_argument("input", type=argparse.FileType("r"), nargs="+",
    help="Input file(s) from search engine")
  parser.add_argument("-o", "--output", type=str, default="./alignment",
    help="Path to output data. Default: './alignment'")

  add_converter_args(parser)
  add_alignment_args(parser)

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

  fileHandler = logging.FileHandler(os.path.join(args.output, "alignment.log"), mode="w")
  fileHandler.setFormatter(logFormatter)
  rootLogger.addHandler(fileHandler)

  consoleHandler = logging.StreamHandler()
  consoleHandler.setFormatter(logFormatter)
  rootLogger.addHandler(consoleHandler)

  logging.info(" ".join(sys.argv[1:]))

  logging.info("Beginning alignment procedure.")

  df, df_original = process_files(args)

  logging.info("Finished converting files")

  align(df, filter_pep=args.filter_pep, mu_min=args.mu_min, rt_distortion=args.rt_distortion, prior_iters=args.prior_iters, stan_iters=args.stan_iters, stan_attempts=args.stan_attempts, stan_file=args.stan_file, print_figures=args.print_figures, output_path=args.output)

  #print(df)
  #sm = StanModel_cache()


