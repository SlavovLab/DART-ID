#!/usr/bin/env python3
# coding: utf-8

import logging
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dart_id.helper import *
from scipy.stats import norm, lognorm, laplace, gaussian_kde

logger = logging.getLogger('root')

def gen(df, config, params, output_path):
  figures_path = create_fig_folder(output_path, 'figures')
  fig_names = []
  col_names = config['col_names']

  # PEP vs. PEP.new scatterplot
  inds = (~pd.isnull(df['pep_new'])) & (df['pep_new'] > 1e-5) & (df[col_names['pep']] > 1e-5)

  f, ax = plt.subplots()
  hst = ax.hist2d(np.log10(df[col_names['pep']][inds]), np.log10(df['pep_new'][inds]), bins=(50, 50), cmap=plt.cm.Reds)
  ax.plot([-5, 0], [-5, 0], '-r')

  ax.grid()

  interval = (-5, 0)
  ax.set_xlim(interval)
  ax.set_ylim(interval)
  ax.set_xlabel('Spectral PEP', fontsize=16)
  ax.set_ylabel('DART-ID PEP', fontsize=16)
  interval = np.arange(-5, 1, 1)
  ax.set_xticks(interval)
  ax.set_yticks(interval)
  ax.set_xticklabels(['$10^{{{}}}$'.format(i) for i in interval], fontsize=12)
  ax.set_yticklabels(['$10^{{{}}}$'.format(i) for i in interval], fontsize=12)

  f.set_size_inches(7, 6)
  # plt.tight_layout()

  cbar = plt.colorbar(hst[3], ax=ax)
  cbar.set_label('Frequency', fontsize=16, labelpad=20, ha='center', va='top')
  cbar.ax.xaxis.set_label_position('top')

  fname = os.path.join(figures_path, 'pep_new_scatterplot.png')
  logger.info('Saving figure to {} ...'.format(fname))
  f.savefig(fname, dpi=160)
  fig_names.append(fname)

  plt.close()
  f.clf()

  # Fold-change increase
  num_points=100
  x = np.logspace(-5, 0, num=num_points)
  y = np.zeros(num_points)
  y2 = np.zeros(num_points)
  y3 = np.zeros(num_points)
  inds = ~pd.isnull(df['pep_new'])

  for i, j in enumerate(x):
      y[i] = np.sum(df['dart_PEP'] < j) / np.sum(df[col_names['pep']] < j)
      y2[i] = np.sum(df[col_names['pep']] < j) / df.shape[0]
      y3[i] = np.sum(df['dart_PEP'] < j) / df.shape[0]

  f, (ax1, ax2) = plt.subplots(2, 1)

  ax1.semilogx(x, (y*100)-100, '-b')
  # ax1.plot([np.min(x), np.max(x)], [0, 0], '-r', linestyle='dashed', linewidth=2)
  ax1.plot([1e-2, 1e-2], [-1000, 1000], '-k', linestyle='dashed')
  ax1.grid()
  ax1.set_xlim([3e-4, 3e-1])
  ax1.set_ylim([-25, np.max(y)*100-50])
  ax1.set_xlabel('PEP Threshold', fontsize=16)
  ax1.set_ylabel('Increase (%)', fontsize=16)
  ax1.set_title('Increase in confident PSMs', fontsize=16)

  ax2.semilogx(x, y2, '-b', linewidth=1, label='Spectra PEP')
  ax2.semilogx(x, y3, '-g', linewidth=1, label='DART-ID PEP')
  #ax2.fill_between(x, 0, y2)
  ax2.plot([1e-2, 1e-2], [-1000, 1000], '-k', linestyle='dashed')
  ax2.grid()
  ax2.set_xlim([3e-4, 3e-1])
  ax2.set_ylim([0, 1.05])
  ax2.set_xlabel('PEP Threshold', fontsize=16)
  ax2.set_ylabel('Fraction', fontsize=16)
  ax2.set_title('Fraction of confident PSMs', fontsize=16)
  ax2.legend(fontsize=16)
  plt.subplots_adjust(hspace=0.6, wspace=0.3)

  f.set_size_inches(5, 7)
  # plt.tight_layout()

  fname = os.path.join(figures_path, 'fold_change_ids.png')
  logger.info('Saving figure to {} ...'.format(fname))
  f.savefig(fname, dpi=160)
  fig_names.append(fname)

  plt.close()
  f.clf()

  # MS1 Intensity CV Validation
  '''
  prots = df[col_names['leading_protein']]

  # get the contaminant and decoy filter tags, if they exist
  filter_con = next((f for f in config['filters'] if f['name'] == 'contaminant'))
  if filter_con is not None: filter_con = filter_con['tag']
  else:                      filter_con = '!'
  filter_rev = next((f for f in config['filters'] if f['name'] == 'decoy'))
  if filter_rev is not None: filter_rev = filter_rev['tag']
  else:                      filter_rev = '!'

  prots = prots.loc[(~prots.str.contains(filter_con)) & 
                    (~prots.str.contains(filter_rev))]
  prots.reset_index(drop=True)

  prot_list = prots.value_counts()
  prot_list = prot_list.loc[prot_list > 50]

  pep_thresh = 0.05

  cvs = np.zeros((len(prot_list), 3))
  n_cvs = np.zeros((len(prot_list), 3))

  for i, j in enumerate(prot_list):
      prot_name = prot_list.index[i]
      intensities = df[col_names['intensity']][(df[col_names['leading_protein']] == prot_name) & 
                                    (~pd.isnull(df[col_names['intensity']])) & 
                                    (df[col_names['pep']] < pep_thresh)]
      cvs[i][0] = np.std(intensities) / np.mean(intensities)
      n_cvs[i][0] = len(intensities)
      
      intensities = df[col_names['intensity']][(df[col_names['leading_protein']] == prot_name) & 
                                    (~pd.isnull(df[col_names['intensity']])) & 
                                    (df['pep_new'] < pep_thresh) &
                                    (df[col_names['pep']] > pep_thresh)]
      cvs[i][1] = np.std(intensities) / np.mean(intensities)
      n_cvs[i][1] = len(intensities)
      
      dfa = df.loc[(~pd.isnull(df[col_names['intensity']])) & 
                   (df[col_names['pep']] < pep_thresh) &
                   (~df[col_names['leading_protein']].str.contains(filter_con)) & 
                   (~df[col_names['leading_protein']].str.contains(filter_rev))
                  ].sample(n=j)
      cvs[i][2] = np.std(dfa[col_names['intensity']]) / np.mean(dfa[col_names['intensity']])
      n_cvs[i][2] = dfa.shape[0]

  # remove pairs which have NAs in them
  na_rows = np.apply_along_axis((lambda x: np.isnan(x).any()), 1, cvs)
  na_rows = na_rows | (np.apply_along_axis((lambda x: (x == 1).any()), 1, n_cvs))
  cvs = cvs[~na_rows]
  n_cvs = n_cvs[~na_rows]

  x = np.linspace(0, 7.5, 1000)
  cvs_density = gaussian_kde(cvs[:,0])
  new_cvs_density = gaussian_kde(cvs[:,1])
  null_cvs_density = gaussian_kde(cvs[:,2])

  f, ax = plt.subplots()

  ax.plot(x, cvs_density(x), '-b', label='Spectra PEP < {}'.format(pep_thresh))
  ax.plot(x, new_cvs_density(x), '-g', label='Spectra+RT PEP < {}'.format(pep_thresh))
  ax.plot(x, null_cvs_density(x), '-r', label='Null')
  ax.set_xlabel('Intra-Protein MS1 Intensity CV ($\sigma$/$\mu$)', fontsize=16)
  ax.set_ylabel('Density', fontsize=16)
  ax.set_title('Intra-Protein MS1 Intensity CV\n Validation of Newly Boosted Observations', fontsize=16)
  ax.tick_params(labelsize=16)
  ax.legend(fontsize=16)

  f.set_size_inches(7, 7)
  plt.tight_layout()

  fname = os.path.join(figures_path, 'ms1_intensity_validation.png')
  logger.info('Saving figure to {} ...'.format(fname))
  f.savefig(fname, dpi=160)
  fig_names.append(fname)

  plt.close()
  f.clf()
  '''

  return fig_names

