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
  description="",
  long_description="""hey""",
  long_description_content_type="text/markdown",
  url="github.com/asdfasdf",
  author="",
  author_email="",
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
  keywords="bayesian retention time psm peptide spectra",
  project_urls={
    "Documentation": "",
    "Source": "",
    "Tracker": "",
    "Lab Page": ""
  },
  packages=find_packages(exclude=["docs", "example"]),
  #libraries=[],
  #setup_requires=[],
  install_requires=[
    "pystan>=2.17.1.0",
    "numpy>=1.14.3",
    "scipy>=1.0.0",
    "pandas>=0.22.0",
    "matplotlib>=2.1.2"
  ],
  tests_requires=[
    "pytest"
  ],
  extras_require={

  },
  package_data={
    "rtlib": [
      "fits/*" # STAN fits
    ]
  },
  # data outside the package
  # data_files=[('my_data', ['data/data_file'])],
  
  entry_points={
    "console_scripts": [
      ("rtlib_convert", "rtlib.converter:main"),
      ("rtlib_align",   "rtlib.align:main"),
      ("rtlib_update",  "rtlib.update:main")
    ]
  }

)