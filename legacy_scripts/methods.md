---
output:
  pdf_document: default
  html_document: default
---
# RTLib Project Methods

### Pairwise (Ali) Method

#### Experiment Alignment

* Select "reference" experiment, by scoring each experiment in the set by a certain criteria, and then choosing the experiment with the best score
* Align all experiments to reference experiment using robust linear regression, and only using a small subset of PSMs that are very confidently identified (PEP < 5e-4)
    * A good number of experiments get thrown out in this way - if they don't have 10 peptides both confidently identified and in common with the reference experiment (122 experiments -> 80 experiments)
* Store linear regression parameters in ```shift.coeffs.txt```, and store corrected RTs in ```RT.lib.txt```, where each row is a peptide, and each column is an experiment
    * If an experiment has more than one PSM for a peptide, then the median of those RTs is taken and put into the library (as the library cannot hold more than one observation per experiment, per peptide)
    
#### Bayesian Update

* For each experiment, make two subsets:
    * _FORW_ - very confidently identified PSMs (PEP < 0.02)
    * _REV_ - known misidentifications, tagged by MaxQuant with the ```REV_``` tag in the ```Leading razor protein``` column
* For each peptide in _FORW_
    * Calculate the ```dRT```, which is the difference between the corrected retention time (observed, modified by linear regression parameters) and the library retention time
    * ```{r}
    dRT = abs(RT.corrected - RT.lib)
    ```
* For each peptide in _REV_
    * Calculate the ```dRT``` between the corrected retention time, and a random retention time from the library
    * ```{r}
    RT.lib = sample(rt.lib$rt.median, nrow(ev.rev), replace=T)
    dRT = abs(RT.corrected - RT.lib)
    ```
* Create density functions for the _FORW_ and _REV_ distributions
    * _FORW_ should be very sharp at 0, and tail off sharply
    * _REV_ should be almost uniform
* Bayesian update performed by comparing _FORW_ and _REV_ density functions


### STAN Method

Use modeling language STAN to simultaneously align experiments and get peptide-specific retention time distributions

#### Experiment Alignment

Cleanse data, before passing into STAN function

* Only use PSMs with PEP < 0.05
* Remove ```CON``` Proteins (Contaminants)
    * Sample prep contaminants (Keratin, Trypsin)
    * LCMS Contaminants (Albumin - from past runs)
* Remove ```REV``` Reverse matches
* Only use Elite chromatography experiments, run on similar LC systems and methods
* Remove previously identified experiments with very poorly performing LC

Peptides are assigned an ID from MaxQuant, but we also assign them a ```Stan ID```, which is just the sequential order of the peptide sequences, sorted alphabetically. All PSMs also have an ```Observation ID``` that is unique to the PSM - in case experiments have more than one PSM per peptide.

Make sure to note the differences between these IDs to avoid indexing errors

##### Fit 1

Two piece linear regression:

* ```beta_0``` - y-intercept
* ```beta_1``` - slope of first segment
* ```beta_2``` - slope of second segment
* ```split_point``` - time (in minutes) of transition between first and second segments

Peptide-specific retention time distributions:

Retention time distribution modeled as normal distribution

* ```mu``` - canonical retention time for given peptide
* ```muijs``` - linear regression-corrected retention time for given peptide, and given experiment
* ```sigma``` - canonical standard deviation for given peptide
* ```sigma_global``` - global average standard deviation for all peptides

##### Fit 2

Same as Fit 1, except standard deviations modeled to increase with retention time

Peptide-specific retention time distributions:

Retention time distribution modeled as laplace distribution

* ```sigma_slope``` - slope of sigma with respect to retention time (minutes)
* ```sigma_slope_global``` - global average of sigma slope
* ```sigma_intercept``` - y-intercept

#### Bayesian Update - Experiment-Centric

Similar to Pairwise/Ali method

For each experiment:

* Find peptides in this experiment that also have alignment data from STAN
* Adjust canonical retention times for these peptides
    * Use ```muijs``` from STAN output
    * OR
    * Adjust ```mu``` manually using ```beta_0```, ```beta_1```, ```beta_2```, and ```split_point```
* Calculate dRT
    * Adjust canonical retention time based on linear regression parameters for this experiment
    * For experiment-centric alignment, peptide-specific standard deviation plays no role
    ```{r}
    muijs = apply.regression.params(mu)
    ```
    * dRT = abs(adjusted canonical RT - observed RT)
    ```{r}
    dRT = abs(muijs - RT)
    ```
* _FORW_ density function is just the density of dRT
* _REV_ density function:
    * Same as _FORW_, but shuffle observed RTs before calculating dRT, in order to randomly match and create the null distribution
    * Would ideally use ```REV``` PSMs, but these are not fed into the STAN function in order to increase the accuracy of the linear regressions
* Bayesian update compares the _FORW_ and _REV_ density functions

#### Bayesian Update - Peptide-Centric

Build global null density function of retention times (empirical density function), which represents the probability of selecting a certain retention time randomly

For each experiment:

* Find peptides in this experiment that also have alignment data from STAN
* Adjust canonical retention times for these peptides
    * Use ```muijs``` from STAN output
    * OR
    * Adjust ```mu``` manually using ```beta_0```, ```beta_1```, ```beta_2```, and ```split_point```
* For Fit2, adjust peptide standard deviations using the linear regression parameters ```sigma_slope``` and ```sigma_intercept```

For each peptide in this experiment:

* Build density function for the retention time of this peptide, using ```muijs``` and ```sigma```
    * Fit1: Normal distribution
    * Fit2: Laplace distribution (double exponential)
* Compute likelihood of observed RT with respect to the density function
* Compute likelihood of observed RT with respect to global RT density function (null distribution)
* Bayesian update compares correct likelihood and incorrect likelihood



