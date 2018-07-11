# coding: utf-8

import argparse
import logging
import os
import pandas as pd
import pkg_resources
import sys
import yaml

from dart_id.exceptions import *
from dart_id.version import __version__
from shutil import copyfile

logger = logging.getLogger('root')

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

def init_logger(verbose, log_file_path, log_to_file=True):
  # set up logger
  for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler) 
   
  logFormatter = logging.Formatter('%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s')
  logger = logging.getLogger('root')

  if verbose: logger.setLevel(logging.DEBUG)
  else: logger.setLevel(logging.WARNING)

  if log_to_file:
    fileHandler = logging.FileHandler(log_file_path, mode='w')
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)

  consoleHandler = logging.StreamHandler()
  consoleHandler.setFormatter(logFormatter)
  logger.addHandler(consoleHandler)

  logger.info(' '.join(sys.argv[1:]))

  return logger

def create_fig_folder(output_path, fig_folder):
  figures_path = os.path.join(output_path, fig_folder)
  if not os.path.exists(figures_path):
    logger.info('Path for figures folder {} does not exist. Creating...'.format(figures_path))
    os.makedirs(figures_path)
  return figures_path

def load_params_from_file(params_folder):
  # first expand user or any vars
  params_folder = os.path.expanduser(params_folder)
  params_folder = os.path.expandvars(params_folder)

  # load parameters if they are specified in the command line
  params = {}
  logger.info('Using provided alignment parameters. Loading params from {}...'.format(params_folder))
  param_files = ['exp_params.txt', 'pair_params.txt', 'peptide_params.txt']
  for pf in param_files:
    pfp = os.path.join(params_folder, pf)
    if os.path.exists(pfp):
      with open(pfp, 'rb') as f:
        try:
          params[pf.split('_')[0]] = pd.read_csv(pfp, sep='\t')
        except:
          logger.error('Error loading param file')
    else:
      raise ConfigFileError('Params file {} does not exist'.format(pfp))

    logger.info('Loaded \"{}\" params file.'.format(pf.split('_')[0]))

  return params

def add_global_args(parser, add_config_file=True):
  parser.add_argument('-i', '--input', type=argparse.FileType('r'), nargs='+', default=None, help='Input file(s) from search engine')
  #parser.add_argument('-o', '--output', help='Path to converted file. Default: prints to stdout')
  parser.add_argument('-o', '--output', type=str, default=None, help='Path to output data. Default: None')
  parser.add_argument('-v', '--verbose', action='store_true', default=False)
  parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=__version__), help='Display the program\'s version')

  if add_config_file:
    parser.add_argument('--config-file', required=True, type=argparse.FileType('r', encoding='UTF-8'), help='Path to config file. See example/config_example.yaml')

def read_default_config_file():
  # load input file types
  default_config = pkg_resources.resource_stream('rtlib', '/'.join(('config', 'default.yaml')))
  default_config = yaml.load(default_config)
  return default_config

def read_config_file(args, create_output_folder=True):
  with open(args.config_file.name, 'r') as f:
    config = yaml.load(f)
    
  # override config file's input, output, and verbose options
  # if they were specified on the command-line
  if args.input is not None:
    if config['input'] is not None:
      logger.info('Input files specified in both the config file and the command line. Using command-line input files instead.')
    config['input'] = [f.name for f in args.input]

  if args.output is not None:
    if config['output'] is not None:
      logger.info('Output folder specified in both the config file and the command line. Using command-line output folder instead.')
    config['output'] = args.output

  if args.verbose:
    if not config['verbose']:
      logger.info('Verbose specified in the command-line but not the config file. Setting verbose to true anyways.')
      config['verbose'] = args.verbose

  # make sure that we have inputs and outputs before continuing
  if config['input'] is None:
    raise ConfigFileError('No input files specified, in either the config file or the command line. Please provide input files.')

  if config['output'] is None:
    raise ConfigFileError('No output folder specified, in either the config file or the command line. Please provide output folder.')

  # expand user or any vars
  config['output'] = os.path.expanduser(config['output'])
  config['output'] = os.path.expandvars(config['output'])

  # create output folder
  if not os.path.exists(config['output']) and create_output_folder:
    logger.info('Output folder does not yet exist. Creating...')
    os.makedirs(config['output'])

  # copy config file to output folder
  if create_output_folder:
    logger.info('Copying config file to output folder')
    copyfile(args.config_file.name, os.path.join(config['output'], os.path.basename(args.config_file.name)))

  return config

# TODO: do this dumb input check every time the config file is loaded
def check_config_inputs(config):
  pass
    

