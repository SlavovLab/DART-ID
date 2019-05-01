#!/usr/bin/env python3
# coding: utf-8

import logging
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dart_id.helper import *
from scipy.stats import norm, lognorm, laplace

logger = logging.getLogger('root')

def gen(df, config, params, output_path):
  figures_path = create_fig_folder(output_path, 'figures')
  fig_names = []

  num_experiments = len(df['exp_id'].unique())
  df['residual'] = df[config['col_names']['retention_time']] - df['muij']

  plots_per_row = 30
  if num_experiments < plots_per_row:
    plots_per_row = num_experiments

  num_rows = int(np.ceil(num_experiments / plots_per_row))

  resi = []

  for i in range(0, num_rows):
      ax = plt.subplot2grid((num_rows, 1), (i, 0))
      
      if (i + 1) * plots_per_row > num_experiments:
          resi = [df['residual'][(df['exp_id'] == i) & (~pd.isnull(df['residual']))] for i in range((i * plots_per_row), num_experiments)]
          ax.boxplot(resi, showfliers=False)
          ax.set_xticklabels(np.arange((i * plots_per_row), num_experiments, 1))
      else:
          resi = [df['residual'][(df['exp_id'] == i) & (~pd.isnull(df['residual']))] for i in range((i * plots_per_row), ((i + 1) * plots_per_row))]
          ax.boxplot(resi, showfliers=False)
          ax.set_xticklabels(np.arange((i * plots_per_row), ((i + 1) * plots_per_row), 1))
          
      #ax.violinplot(resi, showmedians=True, showextrema=True)
      #ax.boxplot(resi, showfliers=False)
      ax.set_xticks(np.arange(1, plots_per_row + 1, 1))

      ax.set_xlabel('Experiment Number')
      ax.set_ylabel('Residual RT (min)')
      
  #plt.subplots_adjust(hspace=0.6, wspace=0.3)
  #plt.tight_layout()

  # finalize and save figure
  f = plt.gcf()
  #f.text(0.5, 0, 'Experiment Number', fontsize=16, ha='center', va='center')
  #f.text(0.06, 0.5, 'Residual RT (min)', fontsize=16, ha='center', va='center', rotation='vertical')
  f.set_size_inches(12, num_rows * 2)

  fname = os.path.join(figures_path, 'residuals_violin.png')
  logger.info('Saving figure to {} ...'.format(fname))
  f.savefig(fname, dpi=160)
  fig_names.append(fname)

  plt.close()
  f.clf()

  return fig_names
