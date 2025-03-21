---
layout: default
title: Configuration File
nav_order: 4
permalink: docs/config-file
---

# Example Configuration File

```yaml
### DART-ID configuration
### =========================

# List the input files here, or define them on the command line
# when running the tool.
# for example, dart_id -i /path/to/input1.txt /path/to/input2.txt
input:
  - /path/to/dat/FP061A/evidence.txt
  - /path/to/dat/FP062ABCD/evidence.txt
  - /path/to/dat/FP063/evidence.txt
  - /path/to/dat/FP064AG/evidence.txt


# Folder to output to. The updated evidence file, as well as the
# optional parameters files and figures will be deposited here
# this can also be specified on the command line. e.g.:
# -o /path/to/output/folder
#output: /path/to/output/folder

# Print diagnostic figures, as well as an HTML file that allows
# for quick browsing
print_figures: true

## Input Type Options
## ==========================

# column names of the input file
# as of now all input files have to be the same format
# change these as the input file changes,
# e.g., when a different search engine or search engine configuration is used
col_names:
  # These four columns are required. This program will not work without them.

  # Sequence can be the canonical amino-acid sequence, 
  # or the modified/annotated sequence, as provided by the search engine
  sequence: Modified sequence
  # The name of the raw/spectrum file, or a unique identifier for each
  # mass-spec run
  raw_file: Raw file
  # The retention/elution time of the ion, in minutes
  # This can also be in seconds, just make sure you update the priors in
  # model to reflect this change.
  retention_time: Retention time
  # The error probability of the peptide-spectrum-match. can be provided
  # by the search engine or by a separate program, e.g., Percolator
  pep: PEP

  # optional columns, that would be used for filtering or figure generation

  # Used to (optionally) append the ion charge state to the peptide sequence, 
  # so that peptides with different charge states are treated as different
  # peptide species.
  charge: Charge

  # Used to run the Fido protein inference algorithm
  leading_protein: Leading razor protein
  proteins: Proteins

  # The base peak width, i.e., the time range between when an ion 
  # first elutes to when it last elutes. Use this as a quality score 
  # in order to filter out poorly retained ions.
  retention_length: Retention length

  # Unused

  #intensity: Intensity
  #leading_gene: ~
  #genes: ~
  #exclude: ~
  #exp_id: ~
  #peptide_id: ~

## PSM Filters
## ======================

# Filters are used to exclude certain observations (PSMs) from
# the alignment process.
# Remove/comment-out filters from this list that you do not want to have.
filters:
  # Filter out entire raw files, especially if they are of a different run-time,
  # or if the LC for that experiment was problematic. The "expr" field is a regular
  # expression that will be checked against all raw files in the input.
  #- name: exclude_filename 
  #  expr: PS06[1-3][AB]|PS064F

  # Same as above, but as a whitelist instead of a blacklist
  #- name: include_filename
  #  expr: 2018A

  # Provide an exclusion list of UniProt IDs. Any PSM matching this
  # list will be filtered out
  # Either a file, with UniProt IDs separated by line breaks, can be
  # specified with the "file" field, or
  # a list of UniProt IDs can be provided in the "list" field
  #- name: uniprot_exclusion
  #  file: /path/to/list_of_uniprot_ids.txt
  #  list:
  #    - or_you_could
  #    - list_uniprot_ids_here
  #    - P36578
  #    - Q99797

  # Filter out contaminants as marked by the search engine
  # The "tag" option is the pattern used to filter out PSMs
  - name: contaminant
    tag: CON__

  # Filter out decoys as marked by the search engine
  # The "tag" option is the pattern used to filter out PSMs
  # - name: decoy
  #   tag: REV__

  # Filter out PSMs by the retention length, which is defined
  # by some search engines as the time at which this spectra is first
  # observed, to the time this spectra is last observed
  # 
  # If "dynamic" is set to true, then the threshold is a fraction of
  # the maximum RT for that raw file (i.e., the run-time). A value of 0.01
  # denotes that the threshold will be 1% of the total run-time of the experiment.
  - name: retention_length
    dynamic: true
    value: 0.01667

  # Filter out PSMs by their RT ranges in each experiment. This behavior is
  # similar but not exactly the same as the "retention_length" filter.
  # 
  # If "dynamic" is set to true, then the threshold is a fraction of
  # the maximum RT for that raw file (i.e., the run-time). A value of 0.01
  # denotes that the threshold will be 1% of the total run-time of the experiment.
  - name: smears
    dynamic: true
    value: 0.03333



### =======================
### !! ADVANCED SETTINGS !!
### =======================

# Only edit the following settings if you understand their effects
# Please refer to config_annotated.yaml for detailed descriptions for
# each configuration field

# Level of verbosity. Higher numbers = printing more information
# 0 = ERROR
# 1 = WARNING (default)
# 2 = INFO
# 3 = DEBUG
# verbose: 1


## Input
## ==========================

# Column delimiter of the input files. i.e., ',' for CSV, '\t' for tabular
# sep: \t

# The input data is loaded in with pandas, and it doesn't like
# some columns being mostly empty. This needs to be set to false 
# for input formats like MaxQuant.
# low_memory: false

# Instead of running a new STAN alignment, use a set of parameters
# from a previous run. The folder needs to include the three files
# outputted from a run with the "save_params" option on, and this run
# needs to be run with the exact same filters as that previous run.
# (exp_params.txt, peptide_params.txt, pair_params.txt)
# params_folder: /path/to/output_folder_from_prev_run

## Alignment Options
## ==========================

# Which alignment model to use
# Options: 'linear', 'two_piece_linear', 'two_piece_linear_laplace'
# model: 'two_piece_linear_laplace'

# add charge of ion onto the sequence, so that sequences ionized
# with different charge states will be aligned separately.
# 
# Sometimes peptide sequences will form chemical adducts on column
# that can reflect on the charge received by the peptide during the
# ioniziation process, and aligning differently charged peptides can
# account for these chromatographic changes
# add_charge_to_sequence: false

# Number of iterations to run when generating priors
# If the average error when generating priors is too high,
# or prohibitive for STAN, then increase these to get more accurate priors
# prior_iters: 10

# Number of iterations to run for STAN. If STAN is consistently hitting
# its iteration limit without reaching an optima it is happy with,
# then increase this number
# stan_iters: 20000

## Advanced Alignment Options

# Minimum value for mu, a canonical retention time (RT) for a peptide
# mu_min: 1

# Amount to distort RTs when calculating priors. If STAN is erroring out
# because the priors are already too close to the optima, then consider
# slowly increasing this value to give STAN more room to iterate.
# rt_distortion: 0


# Advanced STAN parameters (with cmdstan: https://mc-stan.org/users/interfaces/cmdstan),
# for the LBFGS optimization algorithm
# we recommend leaving these at their defaults.

# Line search step size for first iteration
# init_alpha: 0.001 

# Convergence tolerance on absolute changes in objective function value
# tol_obj: 1.e-12

# Convergence tolerance on relative changes in objective function value
# tol_rel_obj: 10000

# Convergence tolerance on the norm of the gradient
# tol_grad: 1.e-8

# Convergence tolerance on the relative norm of the gradient
# tol_rel_grad: 10000000

# Convergence tolerance on changes in parameter value
# tol_param: 1.e-8

# Amount of history to keep for L-BFGS
# history_size: 5

## Update Options
## ==========================

# DART-ID bootstraps the reference RT (mu), to account for uncertainty
# in the estimation and to penalize mu estimates derived from only a few data points (experiments)
# Ideally we would use the MCMC sampler (STAN) to sample the full posterior, but due to 
# technical/performance constraints we are doing this in python instead

# options -- parametric-mixture, parametric, non-parametric, none
# bootstrap_method: 'parametric_mixture'
# bootstrap_iters: 100

# How to aggregate bootstrapped samples
# The weighted mean uses the PEP of each PSM as the weights
# options -- mean, median, weighted_mean
# mu_estimation: 'median'

## Protein Inference Options
## ==========================

# Run protein inference on the newly updated PSMs with the Fido framework
# https://noble.gs.washington.edu/proj/fido
# Paper in J. Proteome Research: http://dx.doi.org/10.1021/pr100594k

# Most, if not all, parameters described below are also described in detail
# on the Fido website and by the helper tips for the command-line verison of Fido.

# To run protein inference, set this flag to true
# run_pi: true

# Parameters derived from a parameter search and optimizing over an objective
# that minimizes selecting false positives.
# Parameters listed below are the default for fido. Leave these, or specify them,
# to skip the parameter searching step.
# Comment out these three parameters to search for the best set of 3 parameters and
# then run protein inference with those.
# pi_gamma: 0.5
# pi_alpha: 0.1
# pi_beta:  0.01

# Log2 of maximum number of subgraph connected states. Graphs with more states
# than this threshold will be pruned. Increasing this number increases run-time,
# by a lot!
# pi_connected_protein_thresh: 14
# Clean up the peptide sequence string, by removing adjacent amino acids,
# modifications, and also switching isoleucine to leucine.
# pi_clean_peptide_name: false
# Default behavior is to cut all PSMs except for the highest scoring one,
# for each peptide, in order to simplify the graph. Set this to true to include
# all PSMs
# pi_use_all_psms: false
# Use protein group level inference
# pi_group_proteins: false
# Prune low-scoring PSMs from the graph before the main pruning procedure.
# The threshold in this case is 1e-2 (PEP > 0.99)
# pi_prune_low_scores: true
# Accuracy of the parameter selection. This will be ignored if pi_gamma, pi_alpha,
# and pi_beta are provided, as the selection will not be performed in the first place.
# 1 = best    / slower     (uses entire data file)
# 2 = relaxed / faster     (uses 300 observations)
# 3 = sloppy  / very fast  (uses 100 observations)
# pi_parameter_accuracy: 3

# Proteins in the "Proteins" column are assumed to be protein IDs in a string,
# separated by a delimiter, which is specified here:
# i.e., the delimiter is ';' if the "Proteins" string is:
#       "Protein1;Protein2;Protein3;Protein4"
# pi_protein_delimiter: ';'
# A substring that delineates decoy proteins. In the case of MaxQuant,
# all decoy proteins are prepended with the string "REV__"
# pi_decoy_tag: 'REV__'


## Output
## ==========================

# Save the parameters outputted by STAN into three text files.
# Use the "params_folder" option in a future run to use these
# parameters instead of running the alignment procedure again.
# save_params: true

# Default behavior is to only append two columns, the new PEP
# and the updated PEP. Set this to true to get many more columns
# added on.
# add_diagnostic_cols: false

# Overwrite the original PEP column with the updated PEP, and save
# the original PEP to the 'Spectra PEP' column.
# Useful for workflows that rely on the PEP column
# overwrite_pep: false

# Remove PSMs that have an FDR (q-value) below this value.
# 0.01 corresponds to selecting PSMs at an FDR of 1%
# psm_fdr_threshold: 0.01

# Remove PSMs that have an associated protein FDR (q-value) below this value.
# 0.01 corresponds to selecting proteins at an FDR of 1%
# protein_fdr_threshold: 0.01


# If providing multiple input files, combine them all into one
# tabular file and save it.
# save_combined_output: true
# The name of the combined output file.
# combined_output_name: ev_updated.txt

# If providing separate input files, then save the output files separately
# as well. This can be used in conjunction with 'save_combined_output'
# save_separate_output: false
# Save the separate output files into the same folder where they originally
# came from. **WARNING** this program does not check to see if it will overwrite
# an existing file. Please choose the options below carefully to avoid overwriting
# your original data!
# save_in_input_folder: false
# The suffix and extension of each of the separate output files. 
# For example, if one of the inputs was "evidence.txt", 
# the output would be "evidence_updated.txt"
# output_suffix: _updated
# output_ext: .txt

# Save logging messages to file?
# log_file: true

## Filters
## ==========================

# Lower threshold of PEP. PSMs with PEP higher than this value will not be
# considered for the alignment process
# These PSMs can still have their confidence updated, as long as there are 
# PSMs of the same sequence that have PEP below this value
# pep_threshold: 0.5

# Peptide sequences need to be observed in at least this number of experiments,
# at a PEP below the pep_threshold, in order to participate in the alignment process
# num_experiments: 3

# Minimum number of confident PSMs per experiment, in order to participate in RT alignment
# If an experiment has less than this number of confident PSMs, then all of its
# PSMs will be excluded from the RT alignment process
# min_psms_per_experiment: 50

```