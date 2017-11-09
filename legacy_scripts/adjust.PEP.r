## init
library(readr)

## load STAN params
#load('dat/params.RData')
load('dat/params_all.RData')

## load evidence
ev <- read_tsv('dat/evidence.txt')
#ev.f <- ev[ev$PEP < 0.05, c('Sequence', 
ev.f <- ev[, c('Sequence',
                            'Proteins', 'Leading razor protein',
                            'Raw file', 'Retention time', 'PEP',
                            'Peptide ID', 'Best MS/MS')]
colnames(ev.f) <- c('Sequence','Proteins','Razor.protein',
                    'Raw.file','Retention.time','PEP', 'Peptide.ID', 'Obs.ID')

# get experiment ids
experiment_factors <- as.factor(ev.f$Raw.file)
experiment_id <- as.numeric(experiment_factors)
num_exps = max(experiment_id)
# get peptide ids
true_peptide_id <- ev.f$Peptide.ID
peptide_id <- as.numeric(as.factor(true_peptide_id))
# observation IDs
obs_id <- as.numeric(ev.f$Obs.ID)

retention_times <- ev.f$Retention.time

# map muij to experiments and peptides
pep_exp_all <- paste(peptide_id, experiment_id, sep=" - ")
pep_exp_pairs <- unique(pep_exp_all)
num_pep_exp_pairs <- length(pep_exp_pairs)
muij_map <- match(paste(peptide_id, experiment_id, sep=" - "), pep_exp_pairs)
splt <- strsplit(pep_exp_pairs, " - ")
muij_to_pep <- as.numeric(sapply(splt, function(x) x[1]))
muij_to_exp <- as.numeric(sapply(splt, function(x) x[2]))

# output table
ev.new <- data.frame(
  Obs.ID=numeric(),
  rt.minus=numeric(),
  rt.plus=numeric(),
  muijs=numeric(),
  sigmas=numeric(),
  PEP.new=numeric()
)

print('Processing...')
for (i in 1:num_exps) {
  # counter also doubles as experiment_id
  exp_id <- i
  # get exp subset of ev.f
  exp <- subset(ev.f, experiment_id==exp_id)
  exp_name <- exp$Raw.file[1]
  
  cat('\r', i, '/', num_exps, exp_name, '                           ')
  flush.console()
  
  # some experiments have very few rows??? skip these
  if(dim(exp)[1] < 10) next
  
  # unique peptide_ids of peptides in this experiment
  exp.peptide.map = peptide_id[experiment_id==exp_id]
  exp.peptides = unique(exp.peptide.map)
  # also get true_peptide_id so we can fetch more data from
  # the original data frame
  exp.true.peptide.map = true_peptide_id[experiment_id==exp_id]
  
  # retention times for this experiment
  exp.rts <- retention_times[experiment_id==exp_id]
  # adjusted mean retention times, for peptide i in experiment j
  muijs <- pars[sprintf('muij[%i]', which(muij_to_exp==exp_id))]
  # NOTE: need to disable scientific notation since 
  #       weird "e" tricks will break indexing later
  names(muijs) <- format(exp.peptides, trim=TRUE, scientific=FALSE)
  # adjusted retention time standard deviations (not experiment specific)
  sigmas <- pars[sprintf('sigma[%i]', exp.peptides)]
  names(sigmas) <- gsub('[^0-9]', '', names(sigmas))
  
  # PEP.new = P(-|RT) = P(RT|-)*P(-) / (P(RT|-)*P(-) + P(RT|+)*P(+)
  # + <- PSM=Correct
  # - <- PSM=Incorrect
  
  # P(RT|-) = probability of peptides RT, given that PSM is incorrect
  #           calculated from the empirical distribution of all RTs (rt.edf)
  exp.rt.edf <- density(exp.rts)
  exp.rt.minus <- approx(exp.rt.edf$x, exp.rt.edf$y, xout=exp.rts)$y
  
  # P(-) = probability that PSM is incorrect (PEP)
  # P(+) = probability that PSM is correct (1-PEP)
  
  # P(RT|+) = probability that given the correct ID, the RT falls in the
  #           lognormal distribution of RTs for that peptide, for that experiment
  exp.rt.plus <- unlist(sapply(exp.peptides, function(id) {
    # NOTE: more non-scientific formatting as described above
    pid <- format(id, trim=TRUE, scientific=FALSE)
    
    rts <- exp.rts[exp.peptide.map==id]
    muij <- muijs[as.character(pid)]
    sigma <- sigmas[as.character(pid)]
    dnorm(rts, mean=muij, sd=sigma)
  }))
  
  exp.PEP <- exp$PEP
  # sometimes MQ will output PEP > 1, which makes no sense, and will
  # result in negative values for our adjusted PEP
  # set all PEP > 1 to PEP = 1
  exp.PEP[exp.PEP > 1] <- 1
  
  # now we can update the PEP
  PEP.new <- (exp.rt.minus*exp.PEP) / 
    ((exp.rt.minus*exp.PEP) + (exp.rt.plus*(1-exp.PEP)))
  
  #hist(PEP.new[PEP.new<0.1], breaks=seq(0, 0.1, by=0.001), main=i)
  
  # output table for this experiment
  exp.new <- data.frame(
    Obs.ID=as.numeric(exp$Obs.ID),
    rt.minus=as.numeric(exp.rt.minus),
    rt.plus=as.numeric(exp.rt.plus),
    muijs=as.numeric(muijs[match(exp.peptide.map, as.numeric(names(muijs)))]),
    sigmas=as.numeric(sigmas[match(exp.peptide.map, as.numeric(names(sigmas)))]),
    PEP.new=as.numeric(PEP.new)
  )
  # append to master output table
  ev.new <- rbind(ev.new, exp.new)
}

# remove experiments with few observations
# (7/4/17 quick fix until STAN can be run again with better params)
ev <- ev[-which(experiment_id==61),]

# reorder ev.new in the same fashion as the original ev
ev.new.f <- ev.new[order(ev.new$Obs.ID),]
rownames(ev.new.f) <- NULL

# combine ev and ev.new
ev.adjusted <- cbind(ev, ev.new.f)
# remove some columns we dont need
ev.adjusted <- ev.adjusted[,!(names(ev.adjusted) %in% 
                                c('Obs.ID'))]

write.table(ev.adjusted, 'dat/ev.adjusted.txt', sep='\t', row.names=FALSE, quote=FALSE)

ev.ff <- ev.adjusted[, c('Sequence', 'Proteins', 'Leading razor protein', 'Raw file', 
                'Retention time', 'Retention length', 'PIF', 'PEP', 'Intensity',
                'Reporter intensity corrected 0', 'Reporter intensity corrected 1', 
                'Reporter intensity corrected 2', 'Reporter intensity corrected 3', 
                'Reporter intensity corrected 4', 'Reporter intensity corrected 5', 
                'Reporter intensity corrected 6', 'Reporter intensity corrected 7', 
                'Reporter intensity corrected 8', 'Reporter intensity corrected 9',
                'Reverse', 'Peptide ID', 
                'rt.minus', 'rt.plus', 'muijs', 'sigmas', 'PEP.new')]
write.table(ev.ff, 'dat/ev.adj.txt', sep='\t', row.names=FALSE, quote=FALSE)
