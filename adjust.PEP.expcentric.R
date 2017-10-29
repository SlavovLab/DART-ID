# adjust.pep.expcentric.R
# experiment-centric method of adjusting PEPs, using STAN parameters
# i.e., compare distributions of correct vs. incorrect IDs

adjust.pep.expcentric <- function(ev, par.file='dat/params.Fit2.RData', path.out=NULL) {
  library(tidyverse)
  source('lib.R')
  
  ## load evidence
  #ev <- read_tsv('dat/evidence.txt')
  
  
  ## load STAN params
  #load(par.file)
  #load('dat/params_all.RData')
  load.ev(ev, par.file=par.file)
  
  # remove abnormal LC experiments
  # load experiments from correlation testing in similar.lc.R
  exps.lc <- unlist(read_csv('dat/exps.corr.txt')[,2])
  names(exps.lc) <- NULL
  
  ev <- ev %>%
    filter(grepl('[0-9]{6}A', `Raw file`)) %>% # Only use Elite experiments
    filter(`Raw file` %in% exps.lc) # Remove abnormal LC experiments
  
  dens.forw <- list()
  dens.rev <- list()
  dRT.forw <- list()
  dRT.rev <- list()
  
  # Updating PEPs
  ev.new <- data.frame()
  for(i in 1:num_exps) {
    # counter also doubles as experiment_id
    exp_id <- i
    exp_name <- levels(experiment_factors)[exp_id]
    # get exp subset of ev
    exp <- subset(ev, ev$`Raw file`==exp_name, 
                  c('Raw file', 'Sequence', 'PEP', 'Retention time', 
                    'Best MS/MS', 'Peptide ID', 'Leading razor protein'))
    
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
    
    exp.dRT <- exp.mus[match(exp.peptide.map, exp.peptides)] - exp.f$`Retention time`
    
    set.seed(1)
    rev <- exp.mus[match(exp.peptide.map, exp.peptides)] - sample(exp.f$`Retention time`)
    forw <- exp.dRT[exp.f$PEP < 0.02]
    # take the absolute value of the forward and reverse distribution
    forw <- abs(forw)
    rev <- abs(rev)
    
    den.forw <- density(forw, na.rm=TRUE)
    den.rev <- density(rev, na.rm=TRUE)
    
    dens.forw[[i]] <- den.forw
    dens.rev[[i]] <- den.rev
    #dRT.forw[[i]] <- forw
    #dRT.rev[[i]] <- rev
    
    # P(dRT | Correct) = probability that given the correct ID
    #           dRT belongs to distribution of correct IDs
    exp.f$rt.plus <- approx(den.forw$x, den.forw$y, xout=abs(exp.dRT))$y
    
    # P(dRT|-) = probability of peptides dRT, given that PSM is incorrect
    #           belongs to distribution of incorrect IDs
    # P(dRT | Incorrect) = 1/(N-1) sum_{j not equal to i} P(dRT for j | j is correct)
    sum.coef <- 1 / (nrow(exp.f) - 1)
    exp.f$rt.minus <- rep(NA, nrow(exp.f))
    for(j in 1:nrow(exp.f)) {
      exp.f$rt.minus[j] <- sum(exp.f$rt.plus[-j], na.rm=TRUE) 
    }
    exp.f$rt.minus <- sum.coef * exp.f$rt.minus
    #exp.f$rt.minus <- approx(den.rev$x, den.rev$y, xout=abs(exp.dRT))$y
    
    # fix missing data - distributions can be very sparse sometimes
    exp.f$rt.plus[exp.f$rt.plus == 0 | is.na(exp.f$rt.plus)] = .Machine$double.xmin
    exp.f$rt.minus[exp.f$rt.minus == 0 | is.na(exp.f$rt.minus)] = .Machine$double.xmin
    
    #plot(density(exp.f$rt.plus, na.rm=TRUE), col='blue')
    #lines(density(exp.f$rt.minus, na.rm=TRUE), col='red')
    
    # Bayes Theorem
    # PEP.new = P(-|RT) = P(RT|-)*P(-) / (P(RT|-)*P(-) + P(RT|+)*P(+)
    # + <- PSM=Correct
    # - <- PSM=Incorrect
    exp.f$PEP.new <- (exp.f$PEP * exp.f$rt.minus) / 
      ((exp.f$PEP * exp.f$rt.minus) + ((1 - exp.f$PEP) * exp.f$rt.plus))
    
    # output table for this experiment, for updated data
    exp.new <- data.frame(
      Raw.file=as.character(exp.f$`Raw file`),
      Obs.ID=as.numeric(exp.f$`Best MS/MS`),
      rt.minus=as.numeric(exp.f$rt.minus),
      rt.plus=as.numeric(exp.f$rt.plus),
      muijs=as.numeric(exp.mus[match(exp.peptide.map, exp.peptides)]),
      sigmas=as.numeric(exp.sigmas[match(exp.peptide.map, exp.peptides)]),
      PEP.new=as.numeric(exp.f$PEP.new)
    )
    # append non matched data to output table
    exp.new <- rbind(exp.new, data.frame(
      Raw.file=as.character(exp$`Raw file`[!exp.matches]),
      Obs.ID=as.numeric(exp$`Best MS/MS`[!exp.matches]),
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
  ev.new.f <- ev.new[order(ev.new$Obs.ID),]
  rownames(ev.new.f) <- NULL
  
  # combine ev and ev.new
  ev.adjusted <- cbind(ev, ev.new.f)
  # remove some columns we dont need
  ev.adjusted <- ev.adjusted[,!(names(ev.adjusted) %in% 
                                  c('Obs.ID', 'Raw.file'))]
  
  ev.ff <- ev.adjusted[, c('Sequence', 'Proteins', 'Leading razor protein', 'Raw file', 
                           'Retention time', 'Retention length', 'PIF', 'PEP', 'Intensity',
                           'Reporter intensity corrected 0', 'Reporter intensity corrected 1', 
                           'Reporter intensity corrected 2', 'Reporter intensity corrected 3', 
                           'Reporter intensity corrected 4', 'Reporter intensity corrected 5', 
                           'Reporter intensity corrected 6', 'Reporter intensity corrected 7', 
                           'Reporter intensity corrected 8', 'Reporter intensity corrected 9',
                           'Reverse', 'Peptide ID', 
                           'rt.minus', 'rt.plus', 'muijs', 'sigmas', 'PEP.new', 'id')]
  
  if(is.null(path.out)) {
    
  } else {
    write.table(ev.ff, path.out, sep='\t', row.names=FALSE, quote=FALSE)
  }
  return(ev.ff)
}


#par(mfrow=c(2,2))
#for(i in sample.int(num_exps, size=4)) {
#  plot(dens.forw[[i]], col='blue', xlab='dRT (min)', 
#       main=paste0('Correct (blue) vs. Incorrect (red)\n', 'Exp: ', i))
#  lines(dens.rev[[i]], col='red')
#}
