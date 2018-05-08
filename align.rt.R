align.rt <- function(ev, 
                     pep_thresh=0.5,
                     rtl_filter=5,
                     n_exp=3,
                     rt.distortion=10,
                     prior.iters=10,
                     STAN.model='dat/fit_RT3c.RData',
                     STAN.iters=2e4,
                     pars.out=NULL,
                     figs.out=NULL) {
  
  library(scales)
  library(readr)
  library(tidyverse)
  library(splines)
  library(rstan)
  library(rmutil)
  
  cat('Subsetting confident, alignable observations...\n')
  
  ev.f <- subset.exp(ev, pep_thresh=pep_thresh, n_exp=n_exp, rtl_filter=rtl_filter)
  
  cat(nrow(ev.f), '/', nrow(ev), '(',format(nrow(ev.f)/nrow(ev)*100, digits=2),'%) alignable observations after filtering.\n')
  
  num_experiments <- max(ev.f$`exp_id`)
  
  ## "true peptide id" matches peptide id in evidence file
  ## "peptide id" is index in 1:num_peptides for stan
  raw_peptide_id <- ev.f$`Peptide ID`
  stan_peptide_id <- as.numeric(as.factor(ev.f$`Peptide ID`))
  ev.f$`Stan ID` <- stan_peptide_id
  
  num_total_observations <- nrow(ev.f)
  num_peptides <- length(unique(stan_peptide_id))
  
  retention_times <- ev.f$`Retention time`
  exp_id <- ev.f$exp_id
  exp_names <- levels(as.factor(ev.f$`Raw file`))
  
  cat('Building peptide-experiment pairs...\n')
  
  pep_exp_all <- paste(stan_peptide_id, exp_id, sep=" - ")
  pep_exp_pairs <- unique(pep_exp_all)
  num_pep_exp_pairs <- length(pep_exp_pairs)
  
  ## map peptides to stan peptides in stan
  muij_map <- match(pep_exp_all, pep_exp_pairs)
  splt <- strsplit(pep_exp_pairs, " - ")
  muij_to_pep <- as.numeric(sapply(splt, function(x) x[1]))
  muij_to_exp <- as.numeric(sapply(splt, function(x) x[2]))
  
  cat(num_pep_exp_pairs, 'peptide-experiment pairs.\n')
  
  ## Remove zeros
  pep <- ev.f$PEP
  pep <- pep + .Machine$double.eps
  
  ## Gather all data required by the stan model into a list
  data <- list(num_experiments=num_experiments, 
               num_peptides=num_peptides,
               num_pep_exp_pairs=num_pep_exp_pairs, 
               muij_map=muij_map,
               muij_to_exp=muij_to_exp, 
               muij_to_pep=muij_to_pep,
               num_total_observations=num_total_observations,
               experiment_id=exp_id, 
               peptide_id=stan_peptide_id,
               retention_times=retention_times,
               mean_log_rt=mean(log(retention_times)),
               sd_log_rt=sd(log(retention_times)),
               mean_rt=mean(retention_times),
               sd_rt=sd(retention_times),
               pep=pep,
               max_retention_time=max(retention_times))
  
  # fix randomness
  set.seed(0)
  
  cat('Initializing fit priors...\n')
  
  # Initialize canonical retention time priors for each peptide
  muInit <- sapply(unique(stan_peptide_id), function(pid) {
    # mean of peptide retention times, weighted by PEP
    weights <- ((1 - pep[stan_peptide_id == pid]) - (1 - pep_thresh)) / pep_thresh
    sum(retention_times[stan_peptide_id == pid] * weights) / sum(weights) + rnorm(1, 0, rt.distortion)
  })
  # negative or very low retention times not allowed. floor at 1 minute
  mu.min <- 1
  muInit[muInit <= mu.min] <- mu.min
  # canonical retention time shouldn't be bigger than largest real RT
  mu.max <- max(retention_times)
  muInit[muInit > mu.max] <- mu.max
  
  # take retention times and distort by +- 10 mins
  rt_distorted <- retention_times + rnorm(num_total_observations, 0, rt.distortion)
  # make sure distorted retention times stay within bounds of real ones
  rt_distorted[rt_distorted > max(retention_times)] <- max(retention_times)
  rt_distorted[rt_distorted < min(retention_times)] <- min(retention_times)
  
  # initialize priors for the segmented linear regression
  # first element of vector is beta_0, or the intercept
  # second element is beta_1 and beta_2, the slopes of the two segments
  beta_init <- rbind(rep(10, num_experiments), rep(1, num_experiments))
  
  cat('Optimizing priors with linear approximation for', prior.iters, 'iterations.\n')
  
  for(i in 1:prior.iters) {
    # for each experiment, fit a simple linear regression
    # between the distorted RTs and the initial canonical retention times
    lm_coefs <- sapply(sort(unique(exp_id)), function(id) {
      rt_cur <- rt_distorted[exp_id == id] 
      mu_cur <- muInit[match(stan_peptide_id[exp_id == id], unique(stan_peptide_id))]
      pep_cur <- pep[exp_id == id]
      lm_cur <- lm(rt_cur ~ mu_cur, weights = 1-pep_cur)
      
      lm_cur$coefficients
    })
    # store linear regression parameters
    beta_init[1,] <- lm_coefs[1,]
    beta_init[2,] <- lm_coefs[2,]
    # store this set of priors before iterating again
    muInitPrev <- muInit
    
    # calculate new set of canonical RTs based on linear regression params
    mu_pred <- (rt_distorted - beta_init[1, exp_id]) / beta_init[2, exp_id] 
    # make sure new canonical RTs are within same range as distorted RTs
    mu_pred[mu_pred <= 0] <- min(rt_distorted)
    mu_pred[mu_pred >= max(rt_distorted)] <- max(rt_distorted)
    
    muPrev <- muInit
    
    # new set of priors for canonical RTs based on weighted combination of
    # this set of predicted canonical RTs
    muInit <- sapply(unique(stan_peptide_id), function(pid) {
      weights <- ((1 - pep[stan_peptide_id == pid]) - (1 - pep_thresh)) / pep_thresh
      sum(mu_pred[stan_peptide_id == pid] * weights) / sum(weights)
    })
    
    cat('Iter', i, '| Avg. error:', sum(muPrev - muInit)^2 / length(muInit), '\n' )
  }
  
  beta_0 <- beta_init[1, ]
  beta_1 <- beta_2 <- beta_init[2, ]
  # apply lower bound of (-1.0 * min(beta_1) * min(muInit)) to beta_0
  # where (-1.0 * min(beta_1) * min(muInit)) is the lowest possible intercept
  # given the lowest possible mu and lowest possible beta_1
  beta_0[beta_0 <= -1.0 * min(beta_1) * min(muInit)] <- (-1.0 * min(beta_1) * min(muInit)) + 1e-3
  
  # apply upper bound to prior canonical RTs
  muInit[muInit >= max(retention_times)] <- 0.99 * max(retention_times)
  
  # create prior list for STAN
  initList <- list(mu=muInit,
                   beta_0 = beta_0,
                   beta_1 = beta_1,
                   beta_2 = beta_2,
                   sigma_slope=rep(0.1, num_experiments),
                   sigma_intercept=rep(0.1, num_experiments),
                   split_point=rep(median(muInit), num_experiments))
  #print(initList)
  #return(NULL)
  
  #sm <- stan_model(file='fit_RT3c.stan')
  load(STAN.model)
  
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  cat('Running STAN Model...\n')
  
  start <- Sys.time()
  pars <- optimizing(sm, data=data, init=initList, iter=STAN.iters, verbose=TRUE)$par
  cat('STAN Model Finished. Run time:', Sys.time() - start, 'seconds. \n')
  
  if(!is.null(pars.out)) {
    cat('Saving parameters to', pars.out, '\n')
    save(pars, file=pars.out)
  }
  
  if(is.null(figs.out)) return()
  if(!dir.exists(figs.out)) {
    dir.create(figs.out)
  }
  
  cat('Generating fit summaries. Printing in', figs.out, '\n')
  
  muij_fit <- pars[grep("muij", names(pars))]
  sigma_ij <- pars[grep("sigma_ij", names(pars))]
  
  residual_vec <- c()
  exp_vec <- c()
  col_vec <- c()
  
  pep_col_code <- cut(pep, breaks=10)
  
  for(exp in 1:num_experiments) {
    cat('\r', 'Generating summary for exp', exp, '/', num_experiments, '                  ')
    flush.console();
    
    betas <- c(pars[paste0("beta_0[", exp, "]")], 
               pars[paste0("beta_1[", exp, "]")], 
               pars[paste0("beta_2[", exp, "]")])
    
    ## Get canonical retention times
    mu <-  pars[grep("mu\\[", names(pars))]
    
    ## get experiment indices
    exp_indices <- muij_to_exp==exp
    
    ## split points
    split <-  pars[grep(paste0("split_point\\[", exp, "\\]"), names(pars))]
    
    ## estimated muij for experiment
    exp_fit <- muij_fit[exp_indices]
    sigmas <- pars["sigma_intercept[20]"] + pars["sigma_slope[20]"]/100 * sort(exp_fit)
    
    ## predicted values (duplicated because multiple observations)
    predicted <- (muij_fit[muij_map])[exp_indices[muij_map]]
    ## predicted sd (duplicated because multiple observations)
    predicted_sd <- (sigma_ij[muij_map])[exp_indices[muij_map]]
    
    mus <- mu[(muij_to_pep[muij_map])[exp_indices[muij_map]]]
    observed <- retention_times[exp_indices[muij_map]]
    obs_peps <- pep[exp_indices[muij_map]]
    obs_code <- pep_col_code[exp_indices[muij_map]]
    residual <- observed - predicted
    
    residual_vec <- c(residual_vec, residual)
    col_vec <- c(col_vec, obs_code)
    exp_vec <- c(exp_vec, rep(exp, length(residual)))
    
    pdf(file=paste0(figs.out, '/exp', exp, '_', exp_names[exp], '.pdf'), width=8, height=10)
    par(mfrow=c(2,2), oma=c(0, 0, 2, 0))
    
    plot(mus, observed, pch=19, cex=0.2, 
         xlab='Canonical RT', ylab='Observed RT',
         main='Segmented Fit')
    abline(v=split, col="blue", lty=2)
    segments(x0=0, y0=betas[1], x1=split, y1=betas[1]+betas[2]*split, col="green", lwd=1.5)
    segments(x0=split, y0=betas[1]+betas[2]*split, x1=400, y1=betas[1] + betas[2]*split +betas[3]*(400-split), col="red", lwd=1.5)
    
    plot(predicted, observed, pch=19, cex=0.2,
         xlab='Predicted RT', ylab='Observed RT')
    abline(a=0, b=1, col="blue")
    abline(v=split, col="blue")
    
    cols <- rev(heat.colors(10))
    
    plot(predicted, observed-predicted, pch=19, cex=0.2, col=cols[obs_code],
         xlab='Predicted RT', ylab='Residual RT',
         main='Residuals of fit')
    lines(predicted[order(predicted)],
          sapply(predicted_sd, function(s) qlaplace(.025, 0, s))[order(predicted)],
          col="red")
    lines(predicted[order(predicted)],
          sapply(predicted_sd, function(s) qlaplace(.975, 0, s))[order(predicted)],
          col="red")
    abline(v=split, col="blue")
    
    mtext(paste0(exp_names[exp]), outer = TRUE, cex = 1.5)
    
    dev.off()
  }
}

## ----------

adjust.pep <- function(
  ev=NULL,
  pars=NULL,
  pep_thresh=0.5,
  rtl_filter=5,
  n_exp=3,
  out.path=NULL
) {
  # load libraries
  library(tidyverse)
  library(rmutil)
  
  if(is.null(ev)) {
    stop('Please provide an experiment file to operate on')
  } else if (typeof(ev) == 'character') {
    ev <- read_tsv(ev)
  } else if (typeof(ev) == 'list') {
  } else {
    stop('Format of experiment file not recognized. Please provide a data frame or a path to the evidence file')
  }
  
  # check pars file
  if(typeof(pars) == 'character') {
    load(pars)
  } else if (typeof(pars) == 'double') {}
  else {
    stop('Invalid pars argument. Either enter path to params file, or pass the params object itself')
  }
  
  # check output file
  if(typeof(out.path) == 'character') {}
  else {
    stop('Invalid out path. Please enter a valid path to save the modified evidence file')
  }
  
  # subset experiment - in exactly same way it was done in the alignment part
  ev.f <- subset.exp(ev, pep_thresh=pep_thresh, n_exp=n_exp, rtl_filter=rtl_filter)
  
  experiment_factors <- as.factor(ev.f$`Raw file`)
  experiment_ids <- ev.f$exp_id
  num_exps <- length(unique(ev.f$exp_id))
  
  ## "true peptide id" matches peptide id in evidence file
  ## "peptide id" is index in 1:num_peptides for stan
  raw_peptide_id <- ev.f$`Peptide ID`
  pep.id.list <- unique(raw_peptide_id)
  stan_peptide_id <- as.numeric(as.factor(ev.f$`Peptide ID`))
  ev.f$`Stan ID` <- stan_peptide_id
  
  num_total_observations <- nrow(ev.f)
  num_peptides <- length(unique(stan_peptide_id))
  
  retention_times <- ev.f$`Retention time`
  mean.log.rt <- mean(log(retention_times))
  sd.log.rt <- sd(log(retention_times))
  max.rt <- max(retention_times) # global maximum retention time
  
  # parse linear regression params from STAN output
  beta0 <- pars[sprintf('beta_0[%i]', seq(1, num_exps))]
  beta1 <- pars[sprintf('beta_1[%i]', seq(1, num_exps))]
  beta2 <- pars[sprintf('beta_2[%i]', seq(1, num_exps))]
  split.point = pars[sprintf('split_point[%i]', seq(1, num_exps))]
  sigma.slope = pars[sprintf('sigma_slope[%i]', seq(1, num_exps))]
  sigma.intercept = pars[sprintf('sigma_intercept[%i]', seq(1, num_exps))]
  sigma.slope.global = pars['sigma_slope_global']
  mus <- pars[grep('mu\\[', names(pars))]
  
  muij_fit <- pars[grep("muij", names(pars))]
  sigma_ijs <- pars[grep('sigma_ij', names(pars))]
  
  # output table
  ev.new <- data.frame(
    Raw.file=character(),
    `Best MS/MS`=numeric(),
    rt.minus=numeric(),
    rt.plus=numeric(),
    muijs=numeric(),
    sigmas=numeric(),
    PEP.new=numeric()
  )
  
  print('Processing Experiments...')
  for(i in 1:num_exps) {
    # counter also doubles as experiment_id
    exp_id <- i
    exp_name <- levels(experiment_factors)[exp_id]
    # get exp subset of ev
    exp <- subset(ev, ev$`Raw file`==exp_name, 
                  c('Raw file', 'Sequence', 'PEP', 'Retention time', 
                    'Best MS/MS', 'Sequence ID', 'Peptide ID'))
    
    cat('\r', i, '/', num_exps, exp_name, nrow(exp), '                           ')
    flush.console()
    
    # not all peptides in this experiment have data from the model
    # we can only update those that have that data. others will not be touched
    exp.matches <- exp$`Peptide ID` %in% pep.id.list
    exp.f <- subset(exp, exp.matches)
    
    # convert true_peptide_id to stan_peptide_id
    exp.peptide.map <- match(exp.f$`Peptide ID`, pep.id.list)
    exp.peptides <- unique(exp.peptide.map)
    
    # apply segmented linear regression parameters for this experiment
    exp.mus <- mus[exp.peptides]
    exp.mus <- sapply(exp.mus, function(mu) {
      if(mu < split.point[exp_id]) {
        return(beta0[exp_id] + (beta1[exp_id] * mu))
      } else {
        return(beta0[exp_id] + (beta1[exp_id] * split.point[exp_id]) + 
                 (beta2[exp_id] * (mu - split.point[exp_id])))
      }
    })
    # get sigmas (standard deviations) for peptides
    exp.sigmas <- sigma.intercept[exp_id] + 
      sigma.slope[exp_id] / 100 * mus[exp.peptides]
    
    # PEP.new = P(-|RT) = P(RT|-)*P(-) / (P(RT|-)*P(-) + P(RT|+)*P(+)
    # + <- PSM=Correct
    # - <- PSM=Incorrect
    
    # P(RT|-) = probability of peptides RT, given that PSM is incorrect
    #           calculated from the uniform density from 0 to max(RT)
    #exp.rt.minus <- 1 / max(exp.f$`Retention time`) #experiment-specific
    
    # Fit3b+c, lognormal density over all retention times
    exp.rt.minus <- dlnorm(exp.f$`Retention time`, meanlog=mean.log.rt, sdlog=sd.log.rt)
    
    # P(-) = probability that PSM is incorrect (PEP)
    # P(+) = probability that PSM is correct (1-PEP)
    
    # P(RT|+) = probability that given the correct ID, the RT falls in the
    #           lognormal distribution of RTs for that peptide, for that experiment
    #
    # this is defined in fit_RT3.stan as a mixture between the laplace and uniform distribution
    # where the laplace distribution is weighted by 1-PEP
    # and the uniform distribution is weighted by PEP
    # -- summing to a total density of 1
    #
    exp.rt.plus <- unlist(sapply(exp.peptides, function(id) {
      rts <- exp.f$`Retention time`[exp.peptide.map==id]
      muij <- as.numeric(exp.mus[exp.peptides==id])
      sigma <- as.numeric(exp.sigmas[exp.peptides==id])
      pep <- exp.f$PEP[exp.peptide.map==id]
      
      # ensure that pep does not exceed 1
      # will result in incorrect negative densities when applying mixture model
      pep[pep > 1] <- 1
      
      # mixture between laplace and uniform distribution
      # Fit3c - lognormal density + normal
      ((pep) * dlnorm(rts, meanlog=mean.log.rt, sdlog=sd.log.rt)) +
        ((1-pep) * dnorm(rts, mean=muij, sd=sigma))
    }))
    # sometimes rt.plus will go so low that it will round to 0
    # just round this back up to the smallest number R will handle
    exp.rt.plus[exp.rt.plus == 0] = .Machine$double.xmin
    
    exp.PEP <- exp.f$PEP
    # sometimes MQ will output PEP > 1, which makes no sense, and will
    # result in negative values for our adjusted PEP
    # set all PEP > 1 to PEP = 1
    exp.PEP[exp.PEP > 1] <- 1
    
    # now we can update the PEP
    # PEP.new = P(-|RT) = P(RT|-)*P(-) / (P(RT|-)*P(-) + P(RT|+)*P(+)
    # + <- PSM=Correct
    # - <- PSM=Incorrect
    PEP.new <- (exp.rt.minus*exp.PEP) / 
      ((exp.rt.minus*exp.PEP) + (exp.rt.plus*(1-exp.PEP)))
    
    # output table for this experiment, for updated data
    exp.new <- data.frame(
      Raw.file=as.character(exp.f$`Raw file`),
      `Best MS/MS`=as.numeric(exp.f$`Best MS/MS`),
      rt.minus=as.numeric(exp.rt.minus),
      rt.plus=as.numeric(exp.rt.plus),
      muijs=as.numeric(exp.mus[match(exp.peptide.map, exp.peptides)]),
      sigmas=as.numeric(exp.sigmas[match(exp.peptide.map, exp.peptides)]),
      PEP.new=as.numeric(PEP.new)
    )
    # append non matched data to output table
    exp.new <- rbind(exp.new, data.frame(
      Raw.file=as.character(exp$`Raw file`[!exp.matches]),
      `Best MS/MS`=as.numeric(exp$`Best MS/MS`[!exp.matches]),
      rt.minus=NA,
      rt.plus=NA,
      muijs=NA,
      sigmas=NA,
      PEP.new=NA
    ))
    
    # append to master output table
    ev.new <- rbind(ev.new, exp.new)
  }
  
  # reorder ev.new in the same fashion as the original ev
  ev.new <- ev.new[order(ev.new$Best.MS.MS),]
  rownames(ev.new) <- NULL
  
  # combine ev and ev.new
  ev.adjusted <- cbind(ev, ev.new)
  # remove some columns we dont need
  ev.adjusted <- ev.adjusted[,!(names(ev.adjusted) %in% 
                                  c('Raw.file'))]
  
  ev.ff <- ev.adjusted[, c('Sequence', 'Proteins', 'Leading razor protein', 'Raw file', 
                           'Retention time', 'Retention length', 'PIF', 'PEP', 'Intensity',
                           'Reporter intensity corrected 0', 'Reporter intensity corrected 1', 
                           'Reporter intensity corrected 2', 'Reporter intensity corrected 3', 
                           'Reporter intensity corrected 4', 'Reporter intensity corrected 5', 
                           'Reporter intensity corrected 6', 'Reporter intensity corrected 7', 
                           'Reporter intensity corrected 8', 'Reporter intensity corrected 9',
                           'Reverse', 'Peptide ID', 'Sequence ID', 'Modifications',
                           'rt.minus', 'rt.plus', 'muijs', 'sigmas', 'PEP.new', 'Best.MS.MS', 'id')] %>%
    rename(`Best MS/MS`='Best.MS.MS')
  
  write_tsv(ev.ff, out.path)
}


subset.exp <- function(ev, pep_thresh=0.5, n_exp=3, rtl_filter=5) {
  # load exclusion list - common contaminants + keratins
  exclude <- read_lines('pd_exclude.txt')
  
  # ev.f - subset of confident, alignable observations
  ev.f <- ev %>% 
    filter(PEP < pep_thresh) %>%
    filter(!grepl('REV*', `Leading razor protein`)) %>%
    filter(!grepl('CON*',`Leading razor protein`))  %>%
    filter(!grepl(paste(exclude, collapse='|'), Proteins)) %>%
    filter(`Retention length` < rtl_filter) %>%
    #mutate(`Peptide ID`=as.numeric(as.factor(`Modified sequence`))) %>%
    select("Peptide ID", "Raw file", "Retention time", "PEP") %>%
    ## Add factor indices
    mutate(exp_id=as.numeric(as.factor(`Raw file`)))
  
  ## Remove peptides that occur in less than n experiments
  # by default n is 3
  toRemove <- ev.f %>%
    group_by(`Peptide ID`) %>%
    summarise(num.exps=length(unique(`exp_id`))) %>%
    filter(num.exps <= n_exp) %>%
    pull(`Peptide ID`)
  ev.f <- ev.f %>%
    filter(!(`Peptide ID` %in% toRemove))
  
  return(ev.f)
}
