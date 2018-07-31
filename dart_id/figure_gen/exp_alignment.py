#!/usr/bin/env python3
# coding: utf-8

import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rtlib.helper import *
from scipy.stats import norm, lognorm

logger = logging.getLogger('root')

def gen(df, config, params, output_path):
  df = df[-(df['remove'])]

  figures_path = create_fig_folder(output_path, 'figures')
  fig_names = []

  exp_names = np.sort(df[config['col_names']['raw_file']].unique())
  num_experiments = len(exp_names)
  # ceil PEP to 1
  df[config['col_names']['pep']][df[config['col_names']['pep']] > 1] = 1
  # split PEP into 10 bins, for coloring points later
  pep_col_code = pd.cut(df[config['col_names']['pep']], 10)

  # generate figures for each experiment
  for exp in range(0, num_experiments):
    logger.info('Generating Summary for Experiment {} | {}'.format(exp+1, exp_names[exp]))

    exp_params = params['exp'].iloc[exp]
    #exp_indices = params['pair']['muij_to_exp'] == exp
    exp_inds = (df['exp_id'] == exp) & (~pd.isnull(df['pep_new']))

    predicted = df['muij'][exp_inds].values
    predicted_sd = df['sigmaij'][exp_inds].values
    mus = df['mu'][exp_inds].values

    # observed values
    observed = df[config['col_names']['retention_time']][exp_inds].values
    obs_peps = df[config['col_names']['pep']][exp_inds].values
    obs_code = pep_col_code[exp_inds].values
    residual = observed - predicted

    # plot the 2-segment linear fit of mus to observed RTs
    plt.subplot(221)
    plt.scatter(mus, observed, s=1, color='black')
    plt.plot([0, exp_params['split_point']],
             [exp_params['beta_0'], (exp_params['split_point'] * exp_params['beta_1']) + exp_params['beta_0']],
            color='red')
    plt.plot([exp_params['split_point'], 300], 
             [(exp_params['split_point'] * exp_params['beta_1']) + exp_params['beta_0'], (exp_params['split_point'] * exp_params['beta_1']) + ((300-exp_params['split_point']) * exp_params['beta_2']) + exp_params['beta_0']],
            color='green')
    plt.plot(np.repeat(exp_params['split_point'], 2), [-100, 300], color='blue', linestyle='dashed')
    plt.axis([0, mus.max() + 10, exp_params['beta_0']-10, observed.max() + 10])
    plt.title('Experiment {} - {}'.format(exp, exp_names[exp]))
    plt.xlabel('Canonical RT (min)')
    plt.ylabel('Observed RT (min)')

    # plot residual vs observed RTs
    plt.subplot(222)
    plt.scatter(predicted, observed, color='black', s=1)
    plt.plot([0, 300], [0, 300])
    plt.plot(np.repeat(exp_params['split_point'], 2), [-100, 300], color='blue', linestyle='dashed')
    plt.axis([0, predicted.max()+10, 0, observed.max()+10])
    plt.xlabel('Predicted RT (min)')
    plt.ylabel('Observed RT (min)')

    # plot residuals, quantiles of residuals, and color points by PEP
    plt.subplot(223)
    plt.scatter(predicted, residual, s=1, c=pep_col_code.cat.codes.values[exp_inds])
    plt.plot([0, 300], [0, 0], color='blue')
    plt.plot(np.repeat(exp_params['split_point'], 2), [-100, 300], color='blue', linestyle='dashed')
    plt.plot(predicted[np.argsort(predicted)], norm.ppf(0.025, loc=0, scale=predicted_sd)[np.argsort(predicted)], color='red')
    plt.plot(predicted[np.argsort(predicted)], norm.ppf(0.975, loc=0, scale=predicted_sd)[np.argsort(predicted)], color='red')
    plt.axis([predicted.min()-5, predicted.max()+5, residual.min()-5, residual.max()+5])
    cbar = plt.colorbar()
    cbar.ax.set_yticklabels(pep_col_code.cat.categories.values)
    plt.xlabel('Predicted RT (min)')
    plt.ylabel('Residual RT (min)')

    # add some space between subplots
    plt.subplots_adjust(hspace=0.3, wspace=0.3)

    # finalize and save figure
    fig = plt.gcf()
    fig.set_size_inches(7, 7)
    _fname = 'alignment_{}_{}.png'.format(str(exp), exp_names[exp])
    fname = os.path.join(figures_path, _fname)
    logger.info('Saving figure to {} ...'.format(fname))
    fig.savefig(fname, dpi=160)
    fig_names.append(os.path.join('figures', _fname))
    
    plt.close()
    fig.clf()

  return fig_names

