#!/usr/bin/env python3
# coding: utf-8

import argparse
import json
import logging
import numpy as np
import os
import pandas as pd
import pkg_resources
import sys
import time

from dart_id.align import align
from dart_id.converter import process_files
from dart_id.figure_gen import exp_alignment, residuals, newpeps
from dart_id.helper import add_global_args, read_config_file, init_logger, load_params_from_file
from shutil import copyfile
from string import Template

logger = logging.getLogger('root')

def figures(df, config=None, params=None):

  logger.info('Saving figures to {}'.format(config['output']))

  fig_data = {}

  try:
    fig_data['alignment'] = exp_alignment.gen(df, config, params, config['output'])
  except: pass

  try:
    fig_data['residuals'] = residuals.gen(df, config, params, config['output'])
  except: pass
  
  try:
    fig_data['newpeps'] = newpeps.gen(df, config, params, config['output'])
  except: pass

  generate_html(fig_data, config['output'])

  logger.info('Figure generation complete')

def generate_html(fig_data, output_path):
  logger.info('Generating HTML...')

  # build dict for variable injection into the HTML file
  fig_data = { 'data': json.dumps(fig_data) }

  # load HTML template and inject data variable
  logger.info('Loading template HTML and injecting variables...')
  template = pkg_resources.resource_string('dart_id', '/'.join((
    'figure_resources', 'template.html')))
  template = template.decode('utf-8')
  template = Template(template).safe_substitute(fig_data)

  # write HTML template to output directory
  html_out = os.path.join(output_path, 'figures.html')
  logger.info('Writing template HTML file to {}'.format(html_out))
  with open(html_out, 'w') as template_out:
    template_out.write(template)

  # move resource files (CSS and JS)
  logger.info('Moving HTML resource files')
  resource_files = [
    'bootstrap.min.css', 'bootstrap.min.js', 
    'jquery-3.3.1.slim.min.js', 'styles.css'] #,
    #'jquery.dataTables.min.css', 'jquery.dataTables.min.js'
  #]
  for i in resource_files:
    copyfile(pkg_resources.resource_filename(
      'dart_id', '/'.join(('figure_resources', i))), 
    os.path.join(output_path, 'figures', i))


def main():
  # load command-line args
  parser = argparse.ArgumentParser()  
  add_global_args(parser)
  args = parser.parse_args()

  # load config file
  # this function also creates the output folder
  config = read_config_file(args)

  # initialize logger
  init_logger(config['verbose'], os.path.join(config['output'], 'figures.log'))

  # load first input file and replace home and user vars
  f = config['input'][0]
  f = os.path.expanduser(f)
  f = os.path.expandvars(f)

  # read in input files
  df = pd.read_csv(f, sep='\t', low_memory=False)
  params = load_params_from_file(config['params_folder'])

  # TODO: check that the columns we need are present

  # generate figures
  figures(df, config=config, params=params)

if __name__ == '__main__':
  main()