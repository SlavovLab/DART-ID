# coding: utf-8

import logging
import os
import pandas as pd

logger = logging.getLogger()

# set functions taken from https://www.saltycrane.com/blog/2008/01/how-to-find-intersection-and-union-of/

#def unique(a):
#    """ return the list with duplicate elements removed """
#    return list(set(a))

def intersect(a, b):
    """ return the intersection of two lists """
    return list(set(a) & set(b))

def union(a, b):
    """ return the union of two lists """
    return list(set(a) | set(b))

def create_fig_folder(output_path, fig_folder):
  figures_path = os.path.join(output_path, fig_folder)
  if not os.path.exists(figures_path):
    logger.info("Path for figures folder {} does not exist. Creating...".format(figures_path))
    os.makedirs(figures_path)
  return figures_path

def load_params_from_file(params_folder):
  # load parameters if they are specified in the command line
  params = {}
  logger.info("Using provided alignment parameters. Loading params from {}...".format(params_folder))
  param_files = ["exp_params.txt", "pair_params.txt", "peptide_params.txt"]
  for pf in param_files:
    pfp = os.path.join(params_folder, pf)
    if os.path.exists(pfp):
      with open(pfp, "rb") as f:
        try:
          params[pf.split("_")[0]] = pd.read_csv(pfp, sep="\t")
        except:
          logger.error("Error loading param file")

  return params