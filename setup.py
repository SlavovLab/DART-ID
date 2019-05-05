#!/usr/bin/env python3
# coding: utf-8

import io
import os

from dart_id.version import __version__
from os import path
from setuptools import setup, find_packages, Command

name = 'dart_id'
test_name = 'dart_id-fresca'

version=__version__
test_version='2.0.4'

# read the contents of your README file
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
  name=name,
  version=version,
  description='RT Alignment and Peptide ID Confidence Updating for LC-MS/MS Data',
  long_description=long_description,
  long_description_content_type='text/markdown',
  url='https://github.com/SlavovLab/DART-ID',
  author='Albert Chen',
  author_email='chen.alb@husky.neu.edu',
  license='MIT',
  classifiers=[
    'Development Status :: 4 - Beta',
    'Environment :: Console',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    
    'License :: OSI Approved :: MIT License',

    'Programming Language :: Python :: 3 :: Only',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',

    'Topic :: Scientific/Engineering :: Bio-Informatics'
  ],
  keywords='bayesian retention time psm peptide spectra update',
  project_urls={
    'Documentation': 'https://github.com/SlavovLab/DART-ID',
    'Source': 'https://github.com/SlavovLab/DART-ID',
    'Tracker': 'https://github.com/SlavovLab/DART-ID/issues',
    'Lab Page': 'https://web.northeastern.edu/slavovlab/'
  },
  packages=['dart_id'],
  #libraries=[],
  #setup_requires=[],
  install_requires=[
    'pyyaml>=3.12',
    'numpy>=1.14.3',
    'scipy>=1.0.0',
    'pandas>=0.22.0',
    'matplotlib>=2.1.2',
    'networkx>=2.1'
  ],
  #test_requires=[
  #  'pytest'
  #],
  extras_require={

  },
  include_package_data=True,
  # specified in MANIFEST.in instead
  #package_data={
  #  'dart_id': [
  #    'models/*',            # STAN models
  #    'figure_gen/*',        # figure generation scripts
  #    'figure_resources/*',  # figure generation resources
  #    'fido/*'               # fido scripts
  #  ]
  #},
  # data outside the package
  #data_files=data_files,
  
  entry_points={
    'console_scripts': [
      ('dart_id=dart_id.update:main'),
      ('dart_id_convert=dart_id.converter:main'),
      ('dart_id_align=dart_id.align:main'),
      ('dart_id_update=dart_id.update:main'),
      ('dart_id_figures=dart_id.figures:main')#,
#      ('dart_id_collate=dart_id.collate:main')
    ]
  }
)
