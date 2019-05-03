# coding: utf-8

# models.py - define models and their parameters

import logging
import numpy as np
import os
import pandas as pd

from dart_id.exceptions import ConfigFileError
from scipy.stats import norm, lognorm, laplace

logger = logging.getLogger('root')

def generate_inits_linear(dff, config):

  exp_names = np.sort(dff['raw_file'].unique())
  num_experiments = len(exp_names)
  num_peptides = dff['stan_peptide_id'].max() + 1

  # get the average retention time for a peptide, weighting by PEP
  def get_mu(x):
    weights = ((1 - x['pep'].values) - (1 - config['pep_threshold'])) / config['pep_threshold']
    return np.sum(x['retention_time'].values * weights) / np.sum(weights)

  # apply the get_mu function on all peptides and add some distortion
  mu_init = dff.groupby('stan_peptide_id')[['pep', 'retention_time']].apply(get_mu).values + np.random.normal(0, config['rt_distortion'], num_peptides)

  # negative or very low retention times not allowed. floor at 5 minutes
  mu_init[mu_init <= config['mu_min']] = config['mu_min']
  # canonical retention time shouldn't be bigger than largest real RT
  mu_max = dff['retention_time'].max()
  mu_init[mu_init > mu_max] = mu_max

  # take retention times and distort
  if config['rt_distortion'] > 0:
    logger.info('Distorting RTs by {} minutes for initial value generation'.format(config['rt_distortion']))

  rt_distorted = dff['retention_time'] + np.random.normal(0, config['rt_distortion'], len(dff['retention_time']))
  # make sure distorted retention times stay within bounds of real ones
  rt_distorted[rt_distorted > dff['retention_time'].max()] = dff['retention_time'].max()
  rt_distorted[rt_distorted < dff['retention_time'].min()] = dff['retention_time'].min()

  # initialize priors for the segmented linear regression
  # first element of vector is beta_0, or the intercept
  # second element is beta_1 and beta_2, the slopes of the two segments
  beta_init = np.array((
    np.repeat(10, num_experiments), 
    np.repeat(1, num_experiments)), dtype=float)

  logger.info('Optimizing priors with linear approximation for {} iterations.'.format(config['prior_iters']))

  mu_pred = np.zeros(num_peptides)
  # temporary data frame to quickly map over in the loop
  dft = pd.DataFrame(dict(
    stan_peptide_id=dff['stan_peptide_id'], 
    exp_id=dff['exp_id'], 
    pep=dff['pep'], 
    retention_time=mu_init[dff['stan_peptide_id']]))

  for i in range(0, config['prior_iters']):
    # for each experiment, fit a simple linear regression
    # between the distorted RTs and the initial canonical retention times
    for j in range(0, num_experiments):
        idx     = (dff['exp_id'] == j)
        rt_cur  = rt_distorted[idx]
        mu_cur  = mu_init[dff['stan_peptide_id'][idx]]
        pep_cur = dff['pep'][idx]

        # for this experiment, run a linear regression (1 degree polyfit)
        # of the mus and the distorted RTs. store the linear regression params
        m, c = np.polyfit(mu_cur, rt_cur, 1, w=(1 - pep_cur))
        beta_init[(0,1), j] = [c, m]

    # calculate new set of canonical RTs based on linear regression params
    mu_pred = (rt_distorted - beta_init[0][dff['exp_id']]) / beta_init[1][dff['exp_id']] 
    # make sure new canonical RTs are within same range as distorted RTs
    mu_pred[mu_pred <= config['mu_min']] = config['mu_min']
    #mu_pred[mu_pred <= 0] = rt_distorted.min()
    mu_pred[mu_pred >= rt_distorted.max()] = rt_distorted.max()
    dft['retention_time'] = np.copy(mu_pred)

    mu_prev = np.copy(mu_init)

    # new set of priors for canonical RTs based on weighted combination of
    # this set of predicted canonical RTs
    mu_init = dft.groupby('stan_peptide_id')[['pep', 'retention_time']].apply(get_mu).values

    logger.info('Iter {} | Avg. canonical RT shift: {:.5f}'.format(i + 1, 
      pow(np.sum(mu_prev - mu_init), 2) / len(mu_init)))


  # grab linear regression params
  beta_0 = beta_init[0]
  beta_1 = beta_init[1]

  # apply lower bound of (-1.0 * min(beta_1) * min(muInit)) to beta_0
  # where (-1.0 * min(beta_1) * min(muInit)) is the lowest possible intercept
  # given the lowest possible mu and lowest possible beta_1
  beta_0[beta_0 <= (-1 * beta_1.min() * mu_init.min())] = \
    (-1 * beta_1.min() * mu_init.min()) + 1e-3
  # apply lower bound of 0 to beta_1
  beta_1[beta_1 <= 0] = 1e-3
  
  # apply upper bound to prior canonical RTs
  mu_init[mu_init >= dff['retention_time'].max()] = 0.95 * dff['retention_time'].max()

  # create prior list for STAN
  init_list = {
    'mu': mu_init.tolist(),
    #'sigma': sigma_init,
    'beta_0': beta_0.tolist(),
    'beta_1': beta_1.tolist(),
    'sigma_slope': np.repeat(0.1, num_experiments).tolist(),
    'sigma_intercept': np.repeat(0.1, num_experiments).tolist()
  }

  return init_list

def generate_inits_two_piece_linear(dff, config):

  # this is essentially the same process as the linear fit
  init_list = generate_inits_linear(dff, config)

  # set beta_2 (slope of second segment) to the same as the slope of the 1st segment
  beta_2 = np.copy(init_list['beta_1'])
  # apply lower bound of 0 to beta_2
  beta_2[beta_2 <= 0] = 1e-3

  # set beta_2
  init_list['beta_2'] = beta_2.tolist()
  # set split point to be the median canonical RT
  init_list['split_point'] = np.repeat(np.median(init_list['mu']), 
                                       len(init_list['beta_0'])).tolist()

  return init_list

def muij_two_piece_linear(df, mu, params):
  muij = pd.Series(np.zeros(df.shape[0]))
  # if the mu is before the split point, only account for the first segment
  muij[mu <= df['split_point']] = \
    df['beta_0'] + (df['beta_1'] * mu)
  # if the mu is after the split point, account for both segments
  muij[mu  > df['split_point']] = \
    df['beta_0'] + (df['beta_1'] * df['split_point']) + \
    (df['beta_2'] * (mu-df['split_point']))

  return muij

def mu_two_piece_linear(df, mu, params):
  mu_pred = pd.Series(np.zeros(df.shape[0]))

  mu_pred[mu <= df['split_point']] = \
    (df['retention_time'] - df['beta_0']) / df['beta_1']

  mu_pred[mu  > df['split_point']] = \
    ((df['retention_time'] - df['beta_0'] - (df['beta_1']*df['split_point'])) \
    / df['beta_2']) + df['split_point']

  return mu_pred

def muij_linear(exp, exp_id, params):
  muij = params['exp']['beta_0'][exp_id] + \
    params['exp']['beta_1'][exp_id] * exp['mu']
  return muij

def mu_linear(df, params):
  mu_pred = (df['retention_time'] - df['beta_0']) / df['beta_1']
  return mu_pred

def sigmaij_linear_mu(df, params):
  # get sigmaij from the sigma_intercept and sigma_slope parameters for this experiment
  sigmaij = df['sigma_intercept'] + ((df['sigma_slope'] / 100) * df['mu'])
  return sigmaij

# uniform RT density from min(RT) to max(RT)
def uniform_null(exp):
  return 1 / (np.max(exp['retention_time']) - np.min(exp['retention_time']))

# normal RT density, truncated at min(RT) and max(RT)
def normal_null(exp):
  rt_mean = np.mean(exp['retention_time'])
  rt_std = np.std(exp['retention_time'])
  return norm.pdf(exp['retention_time'], loc=rt_mean, scale=rt_std)

def mixture_normal_normal(exp):
  # Fit3d - mixture between two normal densities
  
  rt_mean = np.mean(exp['retention_time'])
  rt_std = np.std(exp['retention_time'])

  comp1 = exp['pep'] * \
    norm.pdf(exp['retention_time'], loc=rt_mean, scale=rt_std)
  comp2 = (1.0 - exp['pep']) * \
    norm.pdf(exp['retention_time'], loc=exp['muij'], scale=exp['sigmaij'])

  return comp1 + comp2

def normal_drt(exp):
  return norm.pdf(exp['retention_time'], loc=exp['muij'], scale=exp['sigmaij'])

def laplace_drt(exp):
  return laplace.pdf(exp['retention_time'], loc=exp['muij'], scale=exp['sigmaij'])

models = {
  'linear': {
    'name': 'linear',
    # model names must be valid C++ class names. no dashes, etc.
    'model_name': 'FitLinear',
    # name of the executable file-folder in the models/ folder
    'model_binary': 'fit_RT_linear',
    # function to generate initial values
    # should return dict of initial values, where names of initial values
    # are the same as parameter names.
    # values should be in vanilla python list format, 
    # not numpy arrays or pandas series
    'init_func': generate_inits_linear,
    # when serializing the alignment results, which values to store in which file
    # these should be the same name as the STAN parameter names
    'exp_keys': ['beta_0', 'beta_1', 'sigma_intercept', 'sigma_slope'],
    'pair_keys': ['muij', 'sigma_ij'],
    'peptide_keys': ['mu'],
    # transform reference RT from reference space to absolute RT space
    'ref_to_rt': muij_linear,
    # inverse of above. transform absolute RTs into reference space
    'rt_to_ref': mu_linear,

    'sigmaij_func': sigmaij_linear_mu,
    # function for determining density over all RTs in experiment
    # P(RT|PSM-)
    'rt_minus_func': normal_null,
    # function for likelihood of RT given correct assignment
    # P(RT|PSM+)
    #'rt_plus_func': mixture_normal_normal
    'rt_plus_func': normal_drt
  },
  'two_piece_linear': {
    'name': 'two_piece_linear',
    'model_name': 'FitTwoPieceLinear',
    'model_binary': 'fit_RT3e',
    'init_func': generate_inits_two_piece_linear,
    'exp_keys': ['beta_0', 'beta_1', 'beta_2', 
      'split_point', 'sigma_intercept', 'sigma_slope'],
    'pair_keys': ['muij', 'sigma_ij'],
    'peptide_keys': ['mu'],
    'ref_to_rt': muij_two_piece_linear,
    'rt_to_ref': mu_two_piece_linear,
    'sigmaij_func': sigmaij_linear_mu,
    'rt_minus_func': normal_null,
    'rt_plus_func': normal_drt
  },
  'two_piece_linear_laplace': {
    'name': 'two_piece_linear_laplace',
    'model_name': 'FitTwoPieceLinearLaplace',
    'model_binary': 'fit_RT3e',
    'init_func': generate_inits_two_piece_linear,
    'exp_keys': ['beta_0', 'beta_1', 'beta_2', 
      'split_point', 'sigma_intercept', 'sigma_slope'],
    'pair_keys': ['muij', 'sigma_ij'],
    'peptide_keys': ['mu'],
    'ref_to_rt': muij_two_piece_linear,
    'rt_to_ref': mu_two_piece_linear,
    'sigmaij_func': sigmaij_linear_mu,
    'rt_minus_func': normal_null,
    'rt_plus_func': laplace_drt
  }
}

def get_model_from_config(config):
  model = 'two_piece_linear'
  if config['model'] is not None:
    if config['model'] in models:
      model = config['model']
    else:
      raise ConfigFileError('Model \"{}\" not found. Available choices are: {}'.format(
        model, models.keys()))
  else:
    logger.info('Alignment model not defined. Defaulting to \"two_piece_linear\" model')

  return models[model]

