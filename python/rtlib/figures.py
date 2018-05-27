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

from rtlib.align import add_alignment_args, align
from rtlib.converter import add_converter_args, process_files
from rtlib.figure_gen import *
from rtlib.helper import load_params_from_file
from scipy.stats import norm, lognorm, laplace
from shutil import copyfile
from string import Template

logger = logging.getLogger()

def generate_html():
  logger.info("Generating HTML...")

def figures(df, config=None, params=None, output_path=None):
  #logger.info("hello world")
  #print(df)
  logger.info("Saving figures to {}".format(output_path))

  fig_data = {}

  #fig_data["alignment"] = exp_alignment.gen(df, config, params, output_path)
  #fig_data["residuals"] = residuals.gen(df, config, params, output_path)
  #fig_data["newpeps"] = newpeps.gen(df, config, params, output_path)

  generate_html(fig_data, output_path)

def generate_html(fig_data, output_path):
  #fig_data = {
  #  "data": json.dumps(fig_data)
  #}
  #print(fig_data["data"])
  fig_data = {
    "data": '{"alignment": ["/Users/albert/git/RTLib/Alignments/FP17_20180520_2/figures/alignment_0_180413S_X_FP17A.png", "/Users/albert/git/RTLib/Alignments/FP17_20180520_2/figures/alignment_1_180413S_X_FP17B.png", "/Users/albert/git/RTLib/Alignments/FP17_20180520_2/figures/alignment_2_180413S_X_FP17C.png", "/Users/albert/git/RTLib/Alignments/FP17_20180520_2/figures/alignment_3_180413S_X_FP17F.png", "/Users/albert/git/RTLib/Alignments/FP17_20180520_2/figures/alignment_4_180413S_X_FP17J.png", "/Users/albert/git/RTLib/Alignments/FP17_20180520_2/figures/alignment_5_180413S_X_FP17D.png", "/Users/albert/git/RTLib/Alignments/FP17_20180520_2/figures/alignment_6_180413S_X_FP17G.png", "/Users/albert/git/RTLib/Alignments/FP17_20180520_2/figures/alignment_7_180413S_X_FP17H.png", "/Users/albert/git/RTLib/Alignments/FP17_20180520_2/figures/alignment_8_180413S_X_FP17I.png", "/Users/albert/git/RTLib/Alignments/FP17_20180520_2/figures/alignment_9_180413S_X_FP17K.png", "/Users/albert/git/RTLib/Alignments/FP17_20180520_2/figures/alignment_10_180413S_X_FP17E.png"], "residuals": ["/Users/albert/git/RTLib/Alignments/FP17_20180520_2/figures/residuals_violin.png"], "newpeps": ["/Users/albert/git/RTLib/Alignments/FP17_20180520_2/figures/pep_new_scatterplot.png", "/Users/albert/git/RTLib/Alignments/FP17_20180520_2/figures/fold_change_ids.png", "/Users/albert/git/RTLib/Alignments/FP17_20180520_2/figures/ms1_intensity_validation.png"]}'
  }

  # load template
  template = pkg_resources.resource_string("rtlib", "/".join(("figure_resources", "template.html")))
  template = template.decode('utf-8')
  template = Template(template).safe_substitute(fig_data)

  with open(os.path.join(output_path, "figures.html"), "w") as template_out:
    template_out.write(template)

  resource_files = ["bootstrap.min.css", "bootstrap.min.js", "jquery-3.3.1.slim.min.js", "styles.css"]
  for i in resource_files:
    copyfile(pkg_resources.resource_filename("rtlib", "/".join(("figure_resources", i))), os.path.join(output_path, "figures", i))


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("input", type=str, help="Updated evidence file with diagnostic columns")
  parser.add_argument("-p", "--params-folder", type=str, default=None, help="Path to folder with params files from alignment")
  parser.add_argument("-o", "--output", type=str, default="./rt_update", help="Path to output data. Default: './rt_update'")

  add_converter_args(parser)

  args = parser.parse_args()


  # load input file types
  input_types = pkg_resources.resource_stream("rtlib", "/".join(("config", "input_types.json")))
  input_types = json.load(input_types)["input_types"]
  
  # get the input config of the specified file type
  config = input_types[args.type]
  logger.info("Parsing input files as the output from {}".format(config["name"]))

  df = pd.read_csv(args.input, sep="\t", low_memory=False)
  params = load_params_from_file(args.params_folder)

  # check whether the columns specified in the input types configuration
  # are present in the loaded dataframe
  col_names = config["col_names"]
  for j in list(col_names.keys()):
    col = col_names[j]

    # prep our not found error
    not_found_error = ValueError("Column \"{}\" of column type \"{}\" not present in input file of type \"{}\". Please check that the input file is of the specified type. If column names are still not matching up, contact the package author.".format(col, j, config["name"]))

    # each column can either be a:
    # - string
    # - list of strings
    # - empty list
    if type(col) is str:
      # if the column name is a string, then there is only one possible
      # column name in the dataframe for this column.
      if col not in df.columns:
        raise not_found_error

    elif type(col) is list and len(col) >  0:
      # if the column name is a list of strings, then check each one
      # if one is found, then overwrite the config dict's col_name list
      # with the found string
      # not going to test if two or more match here. the first one will 
      # be taken and that's it
      found = False
      for col_ in col:
        if col_ in df.columns:
          found = True
          config["col_names"][j] = col_
      if not found:
        raise not_found_error

    elif type(col) is list and len(col) == 0:
      # if the column name is an empty list, then
      # it doesn't exist in this input type. skip...
      pass
    
    else:
      # there's something wrong with the input types config file
      raise ValueError("Unknown column configuration in input_types.json. Please contact the package author.")

  figures(df, config=config, params=params, output_path=args.output)

if __name__ == "__main__":
  main()