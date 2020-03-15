#!/usr/bin/env python3
# coding: utf-8

import argparse
import json
import logging
import numpy as np
import os
import pandas as pd
import time

from dart_id.align import align
from dart_id.converter import process_files
from dart_id.exceptions import ConfigFileError
from dart_id.fido.BayesianNetwork import run_internal
from dart_id.helper import add_global_args, read_config_file, init_logger, load_params_from_file
from dart_id.models import models, get_model_from_config
from dart_id.report import generate_report
from scipy.stats import norm, lognorm, laplace, bernoulli, uniform

logger = logging.getLogger('root')

def update(dfa, params, config):
    dfa = dfa.reset_index(drop=True)

    #logger.info('{} / {} ({:.2%}) confident, alignable observations (PSMs) after filtering.'.format(dff.shape[0], dfa.shape[0], dff.shape[0] / dfa.shape[0]))

    # refactorize peptide id into stan_peptide_id, 
    # to preserve continuity when feeding data into STAN
    dfa['stan_peptide_id'] = dfa['sequence'].map({
        ind: val 
        for val, ind in enumerate(dfa['sequence'].unique())
    })

    num_experiments = dfa['exp_id'].max() + 1
    num_peptides = dfa['peptide_id'].max() + 1
    exp_names = np.sort(dfa['raw_file'].unique())
    pep_id_list = dfa['peptide_id'].unique()

    # validate parameters file. make sure it is from the same filters
    # or else the program will crash in the code below
    # check num_experiments, num_peptides
    if (
            params['exp'].shape[0] != num_experiments or 
            params['peptide'].shape[0] != (dfa['stan_peptide_id'].max() + 1)
        ):
        error_msg = 'Parameters files have different data than the input data provided. Ensure that both the input list and filters used to generate the alignment parameters and those provided to the current update are the __exact__ same.'
        raise ConfigFileError(error_msg)
    
    model = get_model_from_config(config)

    # mu from the STAN alignment
    dfa['mu'] = params['peptide']['mu'].values[dfa['stan_peptide_id']] 
    
    # Join transformation parameters (betas, sigmas)
    dfa = (dfa
        .join(params['exp'], on='exp_id', how='left', lsuffix='', rsuffix='_right')
        .drop(columns='exp_id_right')
    )

    # predict mus with RTs, and RTs with aligned mus
    dfa['mu_pred'] = model['rt_to_ref'](dfa, dfa['mu'], params)
    dfa['muij'] = model['ref_to_rt'](dfa, dfa['mu'], params)
    dfa['sigmaij'] = model['sigmaij_func'](dfa, params)
    # scaled sigma is the same ratio of muij / mu applied to sigmaij
    dfa['sigma_pred'] = dfa['sigmaij'] * dfa['mu_pred'] / dfa['muij']

    # get parameters for the null distributions for each experiment
    null_dists = dfa.groupby('exp_id')['retention_time'].agg([np.mean, np.std])
    #null_dists = np.array([norm(loc=null_dists.loc[i, 'mean'], scale=null_dists.loc[i, 'std']) for i in range(0, num_experiments)])
    # first column is mean, second is std
    null_dists = np.array([null_dists['mean'].values, null_dists['std'].values]).T

    # PEP ceiling at 1, otherwise will result in 
    # incorrect negative densities when plugging into Bayes' theorem
    dfa['pep'][dfa['pep'] > 1.0] = 1.0

    # output table
    df_new = pd.DataFrame()

    bootstrap_method = config['bootstrap_method'] if 'bootstrap_method' in config else None

    if bootstrap_method == 'none':
        bootstrap_method = None

    if bootstrap_method is None:
        logger.info('Bootstrap method not defined, using point estimates to update confidence instead.')
    else:
        logger.info('Using \"{}\" bootstrap method'.format(bootstrap_method))

    bootstrap_iters = 20 # default
    if 'bootstrap_iters' in config:
        bootstrap_iters = config['bootstrap_iters']
        if bootstrap_method is not None:
            logger.info('Using {} bootstrap iterations'.format(bootstrap_iters))

    # Calculate the number of observations per peptide
    # Used to determine how many draws we need from our distributions
    obs_per_peptide = (dfa
        .groupby('stan_peptide_id')
        .size()
    )
    max_obs_per_peptide = obs_per_peptide.max()

    laplace_pool = laplace.rvs(size=(max_obs_per_peptide * bootstrap_iters))
    uniform_pool = uniform.rvs(size=(max_obs_per_peptide * bootstrap_iters))
    null_sample_pool = norm.rvs(size=(max_obs_per_peptide * bootstrap_iters))

    # Group all data by peptide
    dfe = dfa.loc[:, ['stan_peptide_id', 'pep', 'mu_pred', 'mu', 'sigma_pred', 'exp_id']]
    dfe_group = dfe.groupby('stan_peptide_id')

    # Extract relevant values for each peptide
    all_mu_preds = dfe_group['mu_pred'].apply(np.array)
    all_mus = dfe_group['mu'].apply(np.array)
    all_sigma_preds = dfe_group['sigma_pred'].apply(np.array)
    all_peps = dfe_group['pep'].apply(np.array)
    all_exp_ids = dfe_group['exp_id'].apply(np.array)

    logger.info('Updating PEPs...')
    for i, e in enumerate(np.sort(dfa['exp_id'].unique())):

        # Timing debugging
        time_init = 0
        time_loo = 0
        time_draw_laplace = 0
        time_scale_laplace = 0
        time_draw_norm_and_uniform = 0
        time_scale_norm_and_uniform = 0
        time_sampling_with_replacement = 0
        time_medians = 0
        time_dist_building = 0
        time_bayes = 0
        time_append = 0

        _time = time.time()

        exp_name = exp_names[i]

        exp = dfa[dfa['exp_id'] == e]
        exp = exp.reset_index(drop=True)

        exp_peptides = exp['stan_peptide_id'].unique()

        logger.info('Exp ({} / {}) - {} - ({} Peptides, {} PSMs)'.format(i + 1, num_experiments, exp_name, len(exp_peptides), exp.shape[0]))

        time_init += (time.time() - _time)

        # vector of P(RT|delta=1) for this experiment.
        rt_plus = pd.Series(np.zeros(exp.shape[0]))

        if bootstrap_method is not None:

            _time = time.time()

            # to avoid using this experiment's own data to update the confidence
            # of its own observations, recalculate the reference RTs (mu) without the
            # data from this experiment, by: 
            # 1) non-parametric bootstrapping over the median of the predicted mus.
            # OR
            # 2) parametric bootstrapping, using the RT distribution parameters
                
            # Extract relevant values for each peptide
            mu_preds = all_mu_preds[exp_peptides]
            mus = all_mus[exp_peptides]
            sigma_preds = all_sigma_preds[exp_peptides]
            peps = all_peps[exp_peptides]
            exp_ids = all_exp_ids[exp_peptides]

            num_peptides = exp_ids.shape[0]
            
            # Leave out this experiment's observations
            leave_out = exp_ids.apply(lambda x: np.array(x == e))
            obs_per_pep = exp_ids.apply(len) - leave_out.apply(sum)

            def loo(x, y):
                return x[~y]

            mu_preds = mu_preds.combine(leave_out, loo)
            mus = mus.combine(leave_out, loo)
            sigma_preds = sigma_preds.combine(leave_out, loo)
            peps = peps.combine(leave_out, loo)
            exp_ids = exp_ids.combine(leave_out, loo)       

            # matrix of n by k estimated mus from the bootstrapping
            # will iterate over in the loop after the immediate one
            mu_k = np.zeros((num_peptides, bootstrap_iters))

            time_loo += (time.time() - _time)

            for j, stan_peptide_id in enumerate(exp_peptides):
                num_obs = obs_per_pep[stan_peptide_id]

                # Parametric bootstrap
                if (
                        bootstrap_method == 'parametric' or 
                        bootstrap_method == 'parametric_mixture' or 
                        bootstrap_method == 'parametric-mixture'
                    ):
                    
                    # Draw num_obs * bootstrap_iters samples
                    _time = time.time()
                    pos_samples = laplace_pool[0:(num_obs * bootstrap_iters)].reshape(bootstrap_iters, num_obs)
                    
                    time_draw_laplace += (time.time() - _time)

                    _time = time.time()
        
                    # Shift and scale sampled RTs by mu and sigma_pred, respectively
                    pos_samples = (pos_samples * sigma_preds[stan_peptide_id]) + mu_preds[stan_peptide_id]
                    
                    time_scale_laplace += (time.time() - _time)

                    if (
                            bootstrap_method == 'parametric_mixture' or 
                            bootstrap_method == 'parametric-mixture'
                        ):

                        _time = time.time()

                        null_samples = null_sample_pool[0:(num_obs * bootstrap_iters)].reshape(bootstrap_iters, num_obs)

                        coin_flips = uniform_pool[0:(num_obs * bootstrap_iters)].reshape(bootstrap_iters, num_obs)

                        time_draw_norm_and_uniform += (time.time() - _time)

                        _time = time.time()

                        # Shift and scale sampled RTs by mean and std of null dists
                        null_samples = (null_samples * null_dists[exp_ids[stan_peptide_id], 1]) + null_dists[exp_ids[stan_peptide_id], 0]

                        fp = coin_flips < np.repeat([peps[stan_peptide_id],], bootstrap_iters, axis=0)

                        # Overwrite original samples with samples from null distribution
                        pos_samples[fp] = null_samples[fp]

                        time_scale_norm_and_uniform += (time.time() - _time)
                

                # Non-parametric bootstrap
                elif (
                        bootstrap_method == 'non-parametric' or
                        bootstrap_method == 'non_parametric'
                    ):
                
                    # Pull random indices from the list of existing predicted mus
                    # To get random indices, just take N random variates from uniform_pool,
                    # (Otherwise used for coin flips in parametric bootstrap)
                    # and multiply by len, then floor, to get a list index
                    # This is just a cheap way to sample with replacement from the mu_preds

                    _time = time.time()

                    # Convert to a numpy array so we can do integer indexing
                    pos_samples = np.array(mu_preds[stan_peptide_id])[
                        np.floor(uniform_pool[0:(num_obs * bootstrap_iters)] * num_obs).astype(int)
                    ]
                    # Reshape into matrix
                    pos_samples = pos_samples.reshape(bootstrap_iters, num_obs)
                    
                    time_sampling_with_replacement += (time.time() - _time)

                _time = time.time()

                # Aggregate all sampled mus and store it in mu_k
                if config['mu_estimation'] == 'median':
                    mu_k[j] = np.median(pos_samples, axis=1)
                elif config['mu_estimation'] == 'mean':
                    mu_k[j] = np.mean(pos_samples, axis=1)
                elif config['mu_estimation'] == 'weighted_mean':
                    # or take the weighted mean
                    weights = ((1 - np.array(peps[stan_peptide_id])) - (1 - config['pep_threshold'])) / config['pep_threshold']
                    mu_k[j] = (np.sum(pos_samples * weights, axis=1) / np.sum(weights))

                time_medians += (time.time() - _time)

            _time = time.time()
            # map of stan_peptide_id onto 1:num_peptides
            pep_inds = {ind: var for var, ind in enumerate(exp_peptides)}
            pep_inds = exp['stan_peptide_id'].map(pep_inds)

            # for each bootstrap iteration:
            for k in range(0, bootstrap_iters):
                # evaluate the transformed RTs (predicted mus) on distributions
                # with the bootstrapped, estimated mus as the means.
                #rt_plus = rt_plus + laplace.pdf(exp['retention_time'], \
                #  loc=model['ref_to_rt'](exp, mu_k[:,j][pep_inds], params), \
                #  scale=exp['sigmaij'])

                rt_plus = rt_plus + laplace.pdf(exp['mu_pred'], 
                    loc=mu_k[:, k][pep_inds],
                    scale=exp['sigma_pred']
                )

            # divide total likelihood by # of iterations to normalize to area of 1
            rt_plus = rt_plus / bootstrap_iters

            time_dist_building += (time.time() - _time)

        else:
            _time = time.time()
            # not using bootstrap, but using adjusted mu as a point estimate
            # for updating the confidence
            rt_plus = model['rt_plus_func'](exp)
            time_dist_building += (time.time() - _time)

        _time = time.time()

        #                                         P(RT|delta=0)*P(delta=0)
        # PEP.new = P(delta=0|RT) =   ---------------------------------------------------
        #                             P(RT|delta=0)*P(delta=0) + P(RT|delta=1)*P(delta=1)
        #                         
        # delta=1 = Correct ID (true positive)
        # delta=0 = Incorrect (false positive)
        
        # P(RT|delta=0) = probability of peptides RT, given that PSM is incorrect
        #           estimate empirical density of RTs over the experiment
        
        rt_minus = model['rt_minus_func'](exp)

        # P(delta=0) = probability that PSM is incorrect (PEP)
        # P(delta=1) = probability that PSM is correct (1-PEP)
        
        # P(RT|delta=1) = probability that given the correct ID, the RT falls in the
        #           normal distribution of RTs for that peptide, for that experiment
          
        # delta=1 = Correct ID (true positive)
        # delta=0 = Incorrect (false positive)
        # 
        pep_new = (
            (rt_minus * exp['pep']) /
            ((rt_minus * exp['pep']) + (rt_plus * (1.0 - exp['pep'])))
        )

        time_bayes += (time.time() - _time)

        _time = time.time()

        # for PSMs for which we have alignment/update data
        exp_new = pd.DataFrame({
            'rt_minus': rt_minus.tolist(),
            'rt_plus': rt_plus.tolist(),
            'mu': exp['mu'].values.tolist(),
            'muij': exp['muij'].values.tolist(),
            'sigmaij': exp['sigmaij'].values.tolist(),
            'pep_new': pep_new.tolist(),

            'id': exp['id'].values,
            'exp_id': exp['exp_id'].values,
            'peptide_id': exp['peptide_id'].values,
            'stan_peptide_id': exp['stan_peptide_id'].values,
            'input_id': exp['input_id'].values,
            'exclude': exp['exclude'].values
        })
        # append to master DataFrame and continue
        df_new = df_new.append(exp_new)

        time_append += (time.time() - _time)

        logger.debug('time_init: {:.1f} ms'.format(time_init*1000))
        logger.debug('time_loo (bootstrap only): {:.1f} ms'.format(time_loo*1000))
        logger.debug('time_draw_laplace (parametric/mixture only): {:.1f} ms'.format(time_draw_laplace*1000))
        logger.debug('time_scale_laplace (parametric/mixture only): {:.1f} ms'.format(time_scale_laplace*1000))
        logger.debug('time_draw_norm_and_uniform (parametric/mixture only): {:.1f} ms'.format(time_draw_norm_and_uniform*1000))
        logger.debug('time_scale_norm_and_uniform (parametric-mixture only): {:.1f} ms'.format(time_scale_norm_and_uniform*1000))
        logger.debug('time_sampling_with_replacement (non-parametric only): {:.1f} ms'.format(time_sampling_with_replacement*1000))
        logger.debug('time_medians (bootstrap only): {:.1f} ms'.format(time_medians*1000))
        logger.debug('time_dist_building: {:.1f} ms'.format(time_dist_building*1000))
        logger.debug('time_bayes: {:.1f} ms'.format(time_bayes*1000))
        logger.debug('time_append: {:.1f} ms'.format(time_append*1000))
        

    # reorder by ID and reset the index
    df_new = df_new.sort_values('id')
    df_new = df_new.reset_index(drop=True)

    return df_new


def write_output(df, out_path, config):
    # remove diagnostic columns, unless they are specified to be kept
    if 'add_diagnostic_cols' not in config or config['add_diagnostic_cols'] == False:
        df_out = df.drop([
            'pep_new', 'participated', 'exclude', 'mu', 'muij', 
            'rt_minus', 'rt_plus', 'sigmaij', 'residual',
            'input_id', 'exp_id', 'peptide_id', 'stan_peptide_id'
        ], axis=1)

    # filter by PSM FDR?
    if 'psm_fdr_threshold' in config and type(config['psm_fdr_threshold']) == float:
        to_remove = (df_out['dart_qval'] > config['psm_fdr_threshold'])
        logger.info('{}/{} ({:.2%}) PSMs removed at a threshold of {:.2%} FDR.'.format(np.sum(to_remove), df_out.shape[0], np.sum(to_remove) / df_out.shape[0], config['psm_fdr_threshold']))
        df_out = df_out[~to_remove].reset_index(drop=True)

    # filter by protein FDR?
    if 'protein_fdr_threshold' in config and type(config['protein_fdr_threshold']) == float:
        if 'razor_protein_fdr' in df_out.columns:
            to_remove = ((df_out['razor_protein_fdr'] > config['protein_fdr_threshold']) | pd.isnull(df_out['razor_protein_fdr']))
            logger.info('{}/{} ({:.2%}) PSMs removed at a threshold of {:.2%} Protein FDR.'.format(np.sum(to_remove), df_out.shape[0], np.sum(to_remove) / df_out.shape[0], config['protein_fdr_threshold']))
            df_out = df_out[~to_remove].reset_index(drop=True)
        
    df.to_csv(out_path, sep='\t', index=False)


def write_parameters_df(df, out_path, config):
    '''Write a stripped down version of the input file, with diagnostic columns
    that the user might've chosen to not output in the main output file

    This parameters.txt file is mainly used to generate the HTML report

    Parameters
    ----------
    df: pandas.DataFrame
    out_path: str
    config: dict

    Returns
    -------
    None
    '''

    logger.info('Saving parameters.csv file')

    input_cols = list(config['col_names'].values())

    diagnostic_cols = [
        'dart_PEP', 'pep_new', 'participated', 'exclude', 
        'mu', 'muij', 'rt_minus', 'rt_plus', 'sigmaij', 'residual', 
        'input_id', 'exp_id', 'peptide_id', 'stan_peptide_id'
    ]

    df_out = df.loc[:, input_cols + diagnostic_cols]

    df_out.to_csv(out_path, sep='\t', index=False)


def main():
    start = time.time()

    # load command-line args
    parser = argparse.ArgumentParser()  
    add_global_args(parser)
    args = parser.parse_args()

    # load config file
    # this function also creates the output folder
    config = read_config_file(args)

    # initialize logger
    init_logger(config['verbose'], os.path.join(config['output'], 'dart.log'), config['log_file'])

    logger.info('Converting files and filtering PSMs')
    df, df_original = process_files(config)
    logger.info('Finished converting files and filtering PSMs.')

    # load params, either from a defined folder or from running the alignment
    params = {}
    if 'params_folder' in config and type(config['params_folder']) is str:
        params = load_params_from_file(config['params_folder'])
    else:
        logger.info('Beginning alignment procedure')
        params = align(df, config)
        logger.info('Alignment procedure finished')

    # now we have the params, run the update
    logger.info('Updating PEPs with alignment data...')
    df_new = update(df, params, config)

    # save the sparse combined input file?
    #df_new.to_csv(os.path.join(args.output, 'df_converted.txt'), sep='\t', index=False)

    # add new columns to original DF, and remove the duplicate ID column
    logger.info('Concatenating results to original data...')
    df_adjusted = pd.concat([
        df_original.loc[~df_original['remove']].reset_index(drop=True),
        df_new.drop(['id', 'input_id'], axis=1).reset_index(drop=True)
    ], axis=1)

    # add rows of PSMs originally removed from analysis
    if np.sum(df_original['remove']) > 0:
        logger.info('Reattaching {} PSMs excluded from initial filters'.format(df_original['remove'].sum()))
        # store a copy of the columns and their order for later
        df_cols = df_adjusted.columns
        # concatenate data frames
        df_adjusted = pd.concat([
                df_adjusted,
                df_original.loc[df_original['remove']]
            ], 
            axis=0, ignore_index=True, sort=True
        )
        # pd.concat reindexes the order of the columns, 
        # so just order it back to what it used to be
        df_adjusted = df_adjusted.reindex(df_cols, axis=1)

    # sort by ID, and reset index
    df_adjusted = df_adjusted.sort_values(['id'])
    df_adjusted = df_adjusted.reset_index(drop=True)

    # add residual RT (alignment error) column
    df_adjusted['residual'] = np.abs(
        df_adjusted[config['col_names']['retention_time']] - df_adjusted['muij']
    )

    # add dart_PEP column - which is pep_new, with the NaNs filled in
    # with the old PEPs.
    df_adjusted['dart_PEP'] = df_adjusted['pep_new']
    df_adjusted['dart_PEP'][pd.isnull(df_adjusted['pep_new'])] = (
        df_adjusted[config['col_names']['pep']][
            pd.isnull(df_adjusted['pep_new'])
        ]
    )
    # make sure that updated PEP does not exceed 1
    df_adjusted['dart_PEP'][df_adjusted['dart_PEP'] > 1] = 1

    # add q-value (FDR) column
    # rank-sorted, cumulative sum of PEPs is expected number of false positives
    # q-value is just that vector divided by # of observations, to get FDR
    logger.info('Calculating FDR (q-values)')
    
    # q-value, without fixing # of false positives to a discrete number
    #df_adjusted['q-value'] = \
    #  ( \
    #    np.cumsum(df_adjusted['dart_PEP'][np.argsort(df_adjusted['dart_PEP'])]) / \
    #    np.arange(1, df_adjusted.shape[0]+1) \
    #  )[np.argsort(np.argsort(df_adjusted['dart_PEP']))]
    
    # q-value, by fixing # of false positives to a discrete number
    # for now, set all null PEPs to 1. we'll remember the index and set them back to nan later
    null_peps = pd.isnull(df_adjusted['dart_PEP'])
    if null_peps.sum() > 0:
        df_adjusted['dart_PEP'][null_peps] = 1

    # get the index order of sorted PEPs
    pep_order = np.argsort(df_adjusted['dart_PEP'])
    # Take the ceiling of the cumulative sum of the sorted PEPs to get the pessimistic
    # estimate of the number of false positives when selecting at that level.
    # because using ceiling, PSMs with different PEPs but within the same relative interval
    # will get the same "num_fp" value.
    num_fp = np.ceil(np.cumsum(df_adjusted['dart_PEP'][pep_order])).astype(int)
    # count the number of occurrences of num_fp and sum them up to get the sample size for each
    # discrete false positive # threshold
    fp_counts = np.cumsum(num_fp.value_counts().sort_index()).values
    # divide # of false positivies by sample size to get q-value. sorting the index order brings
    # the order of values back to their original form
    df_adjusted['dart_qval'] = (num_fp / fp_counts[num_fp-1]).values[np.argsort(pep_order.values)]

    # set null PEPs and q-values back to nan
    if null_peps.sum() > 0:
        df_adjusted['dart_PEP'][null_peps] = np.nan
        df_adjusted['dart_qval'][null_peps] = np.nan

    # rename 'remove' column - which indicates whether or not the PSM participated in the
    # DART-ID alignment and update
    df_adjusted['participated'] = ~df_adjusted['remove']

    ## Run protein inference (fido)?
    if 'run_pi' in config and config['run_pi'] is True:
        logger.info('Running protein inference with Fido...')

        # build fido options into a dict (parameter_map)
        parameter_map = {
            'gamma': config['pi_gamma'] if 'pi_gamma' in config else None,
            'alpha': config['pi_alpha'] if 'pi_alpha' in config else None,
            'beta': config['pi_beta']  if 'pi_beta'  in config else None,

            'connected_protein_threshold': config['pi_connected_protein_thresh'],
            'omit_clean_peptide_name': ~config['pi_clean_peptide_name'],
            'all_psms': config['pi_use_all_psms'],
            'group_proteins': config['pi_group_proteins'],
            'prune_low_scores': config['pi_prune_low_scores'],
            'parameter_accuracy': config['pi_parameter_accuracy'],

            'proteins_column': config['col_names']['proteins'],
            'protein_delimiter': config['pi_protein_delimiter'],
            'leading_protein_column': config['col_names']['leading_protein'],
            'decoy_tag': config['pi_decoy_tag'],

            'sequence_column': config['col_names']['sequence'],
            #'error_prob_column': config['col_names']['pep']
            'error_prob_column': 'dart_PEP',

            # pass in output folder so fido can save some intermediate and output files
            'output': config['output']
        }
        logger.debug('parameter_map for fido:')
        logger.debug(json.dumps(parameter_map, indent=4))

        # run fido subroutine
        df_adjusted = run_internal(df_adjusted, parameter_map)

        logger.info('Fido finished')
        logger.info('FDR for PSM\'s razor protein, from protein inference, placed in \"razor_protein_fdr\" column')

    parameters_path = os.path.join(config['output'], 'parameters.txt')
    write_parameters_df(df_adjusted, parameters_path, config)

    # print figures?
    if config['print_figures']:
        logger.info('Generating IPython/HTML report')
        generate_report(config['output'], {
            'parameters': parameters_path,
            'config': config,
            # Parameters passed to papermill --> IPython must be serializable
            'exp_params': params['exp'].to_json(),
            'peptide_params': params['peptide'].to_json(),
            'pair_params': params['pair'].to_json()
        })

    # overwrite PEP?
    # if true, then store old PEP in "Spectra PEP" column,
    # and put the dart PEP in "PEP" column.
    # then drop the pep_new and dart_PEP columns
    if config['overwrite_pep']:
        logger.info('Overwriting PEP column with new PEP. Saving old PEP in \"Spectra PEP\" column.')
        df_adjusted['Spectra PEP'] = df_adjusted[config['col_names']['pep']]
        df_adjusted[config['col_names']['pep']] = df_adjusted['dart_PEP']
        df_adjusted = df_adjusted.drop(['pep_new', 'dart_PEP'], axis=1)

    # tell the user whether or not to expect diagnostic columns
    if config['add_diagnostic_cols']:
        logger.info('Adding diagnostic columns to output')

    # write to file
    if config['save_combined_output']:
        # if combining input files, then write to one combined file
        out_path = os.path.join(config['output'], config['combined_output_name'])
        logger.info('Combining input file(s) and writing adjusted data file to {} ...'.format(out_path))
        write_output(df_adjusted, out_path, config)
    
    if config['save_separate_output']:
        # if keeping input files separate, then use 'input_id' to retain the
        # order in which the input files were passed in
        logger.info('Saving output to separate files...')
        for i, f in enumerate(config['input']):

            # get output extension
            # default to the same extension as the input
            # if one in the config file exists, use that instead
            out_ext = os.path.splitext(os.path.basename(f))[1]
            if config['output_ext'] is not None:
                out_ext = config['output_ext']

            # contruct output path based on which input file it was
            out_path = os.path.join(
                config['output'], 
                (
                    os.path.splitext(os.path.basename(f))[0] + 
                    config['output_suffix'] + '_' + str(i) + out_ext
                )
            )

            # if saving back to the original input folder,
            # then base the output file name on the input file name instead.
            # no need to number these.
            if config['save_in_input_folder']:
                out_path = os.path.join(
                    os.path.dirname(f),
                    (
                        os.path.splitext(os.path.basename(f))[0] + 
                        config['output_suffix'] + out_ext
                    )
                )

            logger.info('Saving input file {} to {}'.format(i, out_path))
            df_a = df_adjusted.loc[df_adjusted['input_id'] == i]
            # save to file
            # other data formats might have a different separator, or have an index column
            write_output(df_a, out_path, config)

    print('Done! Process took {:.3f} seconds'.format(time.time() - start))

if __name__ == '__main__':
    main()
