# coding: utf-8

import argparse
import json
import logging
import os
import pandas as pd
import pkg_resources
import sys

from rtlib.version import __version__

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

def init_logger(verbose, log_file_path):
  # set up logger
  for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler) 
   
  logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
  logger = logging.getLogger()

  if verbose: logger.setLevel(logging.DEBUG)
  else: logger.setLevel(logging.WARNING)

  fileHandler = logging.FileHandler(log_file_path, mode="w")
  fileHandler.setFormatter(logFormatter)
  logger.addHandler(fileHandler)

  consoleHandler = logging.StreamHandler()
  consoleHandler.setFormatter(logFormatter)
  logger.addHandler(consoleHandler)

  logger.info(" ".join(sys.argv[1:]))

  return logger

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

def add_version_arg(parser):
  parser.add_argument("--version", action="version", version="%(prog)s {version}".format(version=__version__), help="Display the program's version")

def add_config_file_arg(parser):
  parser.add_argument("--config-file", required=True, type=argparse.FileType('r', encoding='UTF-8'), help="Path to config file. See example/config_example.json")

def add_global_args(parser):
  add_version_arg(parser)
  add_config_file_arg(parser)

def read_default_config_file():
  # load input file types
  default_config = pkg_resources.resource_stream("rtlib", "/".join(("config", "default.json")))
  default_config = json.load(default_config)
  return default_config

def read_config_file(args):
  pass

def get_args(parser):
  args = parser.parse_args()

  # check if config file exists
  if args.config_file is not None:
    # read default config file for comparison
    pass

  return args