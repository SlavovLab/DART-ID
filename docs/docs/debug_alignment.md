---
layout: default
title: Debug Alignment Errors
nav_order: 6
permalink: docs/debug-alignment-errors
---

# Debugging Alignment Errors

Occassionally STAN will fail to align the provided data. You might receive one such error:

```
2019-05-05 13:14:57,744 [MainThread  ] [INFO ]  Running STAN with command: fit_RT3e optimize iter=20000 algorithm=lbfgs init_alpha=0.001 tol_obj=1e-12 tol_rel_obj=10000 tol_grad=1e-08 tol_rel_grad=10000000 tol_param=1e-08 history_size=5 init=path/to/stan_init_list.json data file=path/to/stan_input.json output file=path/to/stan_output.csv

...

Initial log joint probability = -145484
    Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
      99      -79470.2     0.0116717   1.03451e+06           1           1      117   
    Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
     199      -70672.1    0.00631696        233576           1           1      220   
...
     899      -53134.8   0.000872197        193198      0.9749      0.9749      959   
    Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
     908      -52983.6   0.000539274        283494   1.081e-12       0.001     1048  LS failed, Hessian reset 
Error evaluating model log probability: Non-finite function evaluation.
Error evaluating model log probability: Non-finite function evaluation.

Optimization terminated with error: 
  Line search failed to achieve a sufficient decrease, no more progress can be made
2019-05-05 13:17:20,874 [MainThread  ] [INFO ]  STAN Model Finished. Run time: 143.130 seconds.
2019-05-05 13:17:20,876 [MainThread  ] [WARNI]  Stan model exited with error 70
```

We will continue with the last set of parameters despite this error, but if you are getting alignment errors it is worth looking into further to ensure that the alignment reached a sufficient optima. Set  ```print_figures: true``` in the configuration file to generate the ```figures.html``` file, and then inspect each experiment to check whether or not the alignment performed as expected. In some cases, STAN will error out even when it hits a good optima, and so these errors can be ignored.

One potential cause for these errors is prematurely getting stuck at a local minima in the optimization space. To prevent this, reduce the number of initial value iterations (```prior_iters``` in the config file), or artificially induce error in the initial values by randomly distorting retention times (```rt_distortion``` in the config file).

Some times errors can be caused by STAN starting in a bad place in the optimization space. If you think this is the problem, simply increase ```stan_attempts``` and hope that you get a lucky seed one attempt.

Finally, one solution that we have found works well is to increase experiment-specific variance by including more experiments -- try even including LC-MS/MS experiments with different conditions (different C18 column, different pump, run months apart, etc.). Adding more variance in this dimension increases the success rate, in our experience.

## It still doesn't work

If you've tried the above and still have issues, then you can run STAN standalone and interpret its debugging output. We have to run the standalone version of stan [```cmdstan```](https://mc-stan.org/users/interfaces/cmdstan). See [this getting started document](https://github.com/stan-dev/cmdstan/wiki/Getting-Started-with-CmdStan) to install ```cmdstan```. 

Run your alignment normally, and in your output folder you should find the ```stan_input.json``` and ```stan_init_list.json``` files that are now generated and saved before alignment begins.

Find the model executable in ```dart_id/models/fit_RT3e/<os>```. For Windows, for example, this file is: ```dart_id/models/fit_RT3e/windows/fit_RT3e.exe```. 

Alternatively, you can build your own executables. Clone the ```cmdstan``` repository (see above getting started link), and our DART-ID repository. Build our model by running:

```bash
cd /path/to/cmdstan
make /path/to/DART-ID/dart_id/models/fit_RT3e
```

### Running the stand-alone alignment model

Now run our compiled model:

```bash
./fit_RT3e optimize iter=20000 save_iterations=1 init='/path/to/stan_init_list.json' data file='/path/to/stan_input.json' output file='output.csv' diagnostic_file='diagnostic.txt'
```

In ```output.csv``` each line of the file corresponds to an iteration, and the columns are the number of parameters. For large alignments this file can be several gigabytes big so we don't recommend trying to open those in R or Excel. Instead extract parameters from each iteration by using the ```tail``` or ```sed``` commands:

```bash
# extract last line (parameters from last iteration, before crash)
tail -n 1 output.csv > output_last_line.csv

# extract range of lines (parameters from range of iterations)
# see: https://stackoverflow.com/q/83329
sed -n '16224,16482p;16483q' output.csv > output_16224_16482.csv
```

If you find a reason why our model is failing to align experiments, please contact the authors! We are actively looking to improve the robustness of our model and welcome any and all feedback.
