#!/bin/sh
# coding: utf-8

if ! [ -x "$(command -v conda)" ]; then
  echo 'Error: conda is not installed.' >&2
  exit 1
fi

if conda env list | grep -q dart; then
  echo 'Found dart environment:'
  conda env list | grep dart
  echo 'Switching to dart environment...'
  env_name=conda env list | grep -o '^dart[^\s\n]*\s'
  source activate $env_name
else
  echo 'dart environment not found. creating...'

  # create conda environment for development
  conda env create -f environment.yml -n dart & source activate dart
fi
