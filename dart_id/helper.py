# coding: utf-8

import argparse
import logging
import numpy as np
import os
import pandas as pd
import pkg_resources
import sys
from yaml import load as yaml_load, dump as yaml_dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


from dart_id.exceptions import ConfigFileError
from dart_id.version import __version__
from jsonschema import validate, ValidationError, Draft7Validator
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
   
  #logFormatter = logging.Formatter('%(asctime)s [%(threadName)-5.5s] [%(levelname)-8.8s]  %(message)s')
  logFormatter = logging.Formatter('%(asctime)s [%(levelname)-.8s]  %(message)s', '%Y-%m-%d %H:%M:%S')
  logger = logging.getLogger('root')

  if   verbose == 3: logger.setLevel(logging.DEBUG)
  elif verbose == 2: logger.setLevel(logging.INFO)
  elif verbose == 1: logger.setLevel(logging.WARNING)
  elif verbose == 0: logger.setLevel(logging.ERROR)
  else: logger.setLevel(logging.WARNING) # default

  if log_to_file:
    fileHandler = logging.FileHandler(log_file_path, mode='w')
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)

  consoleHandler = logging.StreamHandler()
  consoleHandler.setFormatter(logFormatter)
  logger.addHandler(consoleHandler)

  logger.info(' '.join(sys.argv[0:]))

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
      try:
        params[pf.split('_')[0]] = pd.read_csv(pfp, sep='\t')
      except:
        logger.error('Error loading param file')
    else:
      raise ConfigFileError('Params file {} does not exist'.format(pfp))

    logger.info('Loaded \"{}\" params file.'.format(pf.split('_')[0]))

  return params

def add_global_args(parser, add_config_file=True):
  parser.add_argument('-i', '--input', type=argparse.FileType('r'), nargs='+', default=None, help='Input file(s) from search engine output (e.g., MaxQuant evidence.txt). Not required if input files are specified in the config file')
  #parser.add_argument('-o', '--output', help='Path to converted file. Default: prints to stdout')
  parser.add_argument('-o', '--output', type=str, default=None, help='Path to output folder')
  parser.add_argument('-v', '--verbose', type=int, default=1, help='Level of logging. 0 = ERROR, 1 = WARNING (default), 2 = INFO, 3 = DEBUG')
  parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=__version__), help='Display the program\'s version')

  if add_config_file:
    parser.add_argument('-c', '--config-file', required=True, type=argparse.FileType('r', encoding='UTF-8'), help='Path to config file (required). See example/config_example.yaml')

def read_default_config_file():
  # load input file types
  default_config = pkg_resources.resource_stream('dart_id', '/'.join(('config', 'defaults.yaml')))
  default_config = yaml_load(default_config, Loader=Loader)
  return default_config

def read_config_file(args, create_output_folder=True):
  # load defaults
  config = read_default_config_file()

  # override defaults with user config file
  with open(args.config_file.name, 'r') as f:
    config.update(yaml_load(f, Loader=Loader))
    
  # override config file's input, output, and verbose options
  # if they were specified on the command-line
  if args.input is not None:
    if config['input'] is not None:
      logger.warning('Input files specified in both the config file and the command line. Using command-line input files instead.')
    config['input'] = [f.name for f in args.input]

  if args.output is not None:
    if 'output' in config and config['output'] is not None:
      logger.warning('Output folder specified in both the config file and the command line. Using command-line output folder instead.')
    config['output'] = args.output

  if args.verbose:
    if 'verbose' in config:
      logger.info('Overwriting verbosity level in configuration file with the one provided on the command-line.')
      config['verbose'] = args.verbose

  # make sure that we have inputs and outputs before continuing
  # the jsonschema validator will catch this as well but here we can print
  # a more descriptive error message
  if 'input' not in config or config['input'] is None:
    raise ConfigFileError('No input files specified, in either the config file or the command line. Please provide input files.')

  if 'output' not in config or config['output'] is None:
    raise ConfigFileError('No output folder specified, in either the config file or the command line. Please provide output folder.')

  ### --------------------
  ### VALIDATE CONFIG FILE
  ### --------------------

  schema = pkg_resources.resource_stream('dart_id', '/'.join(('config', 'schema.yaml')))
  schema = yaml_load(schema, Loader=Loader)

  v = Draft7Validator(schema)
  errors = sorted(v.iter_errors(config), key=str)
  
  for error in errors:
    logger.error("""Configuration file error:
  In field: {}
  Error: {}
  Field description: {}
""".format(
      ' --> '.join(['\'' + str(x) + '\'' for x in error.path]),
      error.message,
      error.schema['description'] if 'description' in error.schema else 'No field description provided'
    ))

    #for suberror in sorted(error.context, key=lambda e: e.schema_path):
    #  print('suberror')
    #  print(list(suberror.schema_path), suberror.message, sep=", ")

  if len(errors) > 0:
    error_msg = '{} error(s) from configuration file. Please read the validation error messages carefully and fix the configuration file'.format(len(errors))
    raise ConfigFileError(error_msg)

  ### ====================================================
  ### ADVANCED CONFIGURATION FILE VALIDATION
  ### --------------------------------------
  ### apply rules too complex for the jsonschema validator
  ### ====================================================

  # ...

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

# pulled from: https://stackoverflow.com/a/29677616
def weighted_quantile(values, quantiles, sample_weight=None, values_sorted=False, old_style=False):
  """ Very close to numpy.percentile, but supports weights.
  NOTE: quantiles should be in [0, 1]!
  :param values: numpy.array with data
  :param quantiles: array-like with many quantiles needed
  :param sample_weight: array-like of the same length as `array`
  :param values_sorted: bool, if True, then will avoid sorting of initial array
  :param old_style: if True, will correct output to be consistent with numpy.percentile.
  :return: numpy.array with computed quantiles.
  """
  
  values = np.array(values)
  quantiles = np.array(quantiles)

  if sample_weight is None:
    sample_weight = np.ones(len(values))
  sample_weight = np.array(sample_weight)
  assert np.all(quantiles >= 0) and np.all(quantiles <= 1), 'quantiles should be in [0, 1]'

  if not values_sorted:
    sorter = np.argsort(values)
    values = values[sorter]
    sample_weight = sample_weight[sorter]

  weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
  if old_style:
    # To be convenient with np.percentile
    weighted_quantiles -= weighted_quantiles[0]
    weighted_quantiles /= weighted_quantiles[-1]
  else:
    weighted_quantiles /= np.sum(sample_weight)
  return np.interp(quantiles, weighted_quantiles, values)

def weighted_median(values, weights):
  order = np.argsort(values)
  values = values[order]
  weights = weights[order]

  #weighted_quantiles = (np.cumsum(weights) - (0.5 * weights)) / np.sum(weights)
  
  return np.interp([0.5], 
    (np.cumsum(weights) - (0.5 * weights)) / np.sum(weights),
    values)
    
def convert_numpy_scalar(x):
  if type(x) == float or type(x) == int:
    return x
  else: return x.tolist()

