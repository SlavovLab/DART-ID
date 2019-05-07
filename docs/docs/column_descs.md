---
layout: default
title: Column Descriptions
nav_order: 5
permalink: docs/column-descs
---

# Output Column Descriptions

Below are descriptions of all column added by DART-ID onto your input data.

Note, if you are providing peptide-level data (e.g., MaxQuant's ```evidence.txt``` instead of ```msms.txt```), then all values are for the best peptide and not all its PSMs. 

## Confidence Update

| ```dart_PEP``` | The updated posterior error probability (PEP) for a peptide after applying DART-ID. If this peptide did not participate in DART-ID (due to not enough observations, or decoy tag), then this value defaults to its original PEP as assigned by the search engine. |
| ```dart_qval``` | The q-value (false discovery rate, FDR) derived from the ```dart_PEP``` values. To filter at 1% FDR, for example, throw out all peptides with ```dart_qval``` > 0.01. |

## Protein Inference:

The Fido protein inference algorithm stores protein IDs and their associated probabilities in the ```protein_fdr.txt``` file in the output folder. In addition, we provide each peptide with the probability of its razor protein:

| ```razor_protein_fdr``` | The probability that the peptides razor protein (as assigned by the search engine) is present in the sample. |

## Diagnostic Columns

| ```pep_new``` | The same as ```dart_PEP```, except that for peptides that did not participate in DART-ID, this value is NA. |
| ```participated``` | True/False, whether or not this peptide participated in the DART-ID algorithm. |
| ```exclude``` | True/False, whether or not this peptide was excluded from alignment but still allowed to participate in the confidence update process (i.e., if True, this peptide was not used to infer retention times, but still had its error probability updated using the inferred RTs from the model) |
| ```mu``` | The reference RT for this peptide. |
| ```muij``` | The inferred retention time for this peptide i in experiment j, derived from mu and the experiment transforms beta_0, beta_1, beta_2, and split_point. |
| ```sigmaij``` | The standard deviation (sqrt(variance)) of RT for this peptide i in experiment j. Derived from mu and the experimental transforms, sigma_slope and sigma_intercept. |
| ```rt_minus``` | The P(RT \| ID is incorrect) portion of Bayes' Theorem. This is derived from evaluating the observed RT on the null distribution of RTs. |
| ```rt_plus``` | The P(RT \| ID is correct) portion of Bayes' Theorem. This is derived from evaluating the observed RT on the inferred RT distribution (built from muij and sigmaij). |
| ```residual``` | The difference between the observed RT and the inferred RT (muij). |
| ```input_id``` | ID number associated with the input data file (actual file, not raw file name). |
| ```exp_id``` | ID number associated with the raw file name. |
| ```peptide_id``` | ID number associated with the peptide identity. |
| ```stan_peptide_id``` | Re-mapped peptide_id that is fed into the optimization program STAN. |
