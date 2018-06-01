#!/usr/bin/env python3
# coding: utf-8

import io
import os

from os import path
from setuptools import setup, find_packages, Command

name = "rtlib"

setup(
  name=name,
  version="1.0.0",
  description="RT Alignment and Updating for MS/MS Data",
  url="https://github.com/blahoink/RTLib",
  author="Albert Chen",
  author_email="chen.alb@husky.neu.edu",
  license="MIT",
  classifiers=[
    "Development Status :: 3 - Alpha",
    "Environment :: Console",
    "Intended Audience :: Developers",
    "Intended Audience :: Science/Research",
    
    "License :: OSI Approved :: MIT License",

    "Programming Language :: Python :: 3 :: Only",
    "Programming Language :: Python :: 3.4",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",

    "Topic :: Scientific/Engineering :: Bio-Informatics"
  ],
  keywords="bayesian retention time psm peptide spectra update",
  project_urls={
    "Documentation": "https://github.com/blahoink/RTLib",
    "Source": "https://github.com/blahoink/RTLib",
    "Tracker": "https://github.com/blahoink/RTLib/issues",
    "Lab Page": "https://web.northeastern.edu/slavovlab/"
  },
  packages=find_packages(
    package_dir={'':'python'},
    exclude=["docs", "example", "__pycache__", ".ipynb_checkpoints"]),
  #libraries=[],
  #setup_requires=[],
  install_requires=[
    "pyyaml>=3.12",
    "pystan>=2.17.1.0",
    "numpy>=1.14.3",
    "scipy>=1.0.0",
    "pandas>=0.22.0",
    "matplotlib>=2.1.2"
  ],
  #test_requires=[
  #  "pytest"
  #],
  extras_require={

  },
  package_data={
    "rtlib": [
      "fits/*",             # STAN fits
      "figure_gen/*",       # figure generation scripts
      "figure_resources/*"  # figure generation resources
    ]
  },
  # data outside the package
  # data_files=[('my_data', ['data/data_file'])],
  
  entry_points={
    "console_scripts": [
      ("rtlib_convert=rtlib.converter:main"),
      ("rtlib_align=rtlib.align:main"),
      ("rtlib_update=rtlib.update:main"),
      ("rtlib_fiures=rtlib.figures:main")
    ]
  }

)