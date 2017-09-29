## init
library(readr)

## load STAN params
#load('dat/params_elite.RData')
#load('dat/params_all.RData')
load('dat/params.corr.RData')

## load evidence
ev <- read_tsv('dat/evidence.txt')
# only keep experiments done on Elite system
ev <- ev[grep('[0-9]{6}A', ev$`Raw file`),]
# remove abnormal LC experiments
# load experiments from correlation testing in similar.lc.R
exps.lc <- unlist(read_csv('dat/exps.corr.txt')[,2])
names(exps.lc) <- NULL
ev <- ev[ev$`Raw file` %in% exps.lc,]


ev.f <- ev[ev$PEP < 0.05, c('Sequence', 
                            'Proteins', 'Leading razor protein',
                            'Raw file', 'Retention time', 'PEP',
                            'Peptide ID', 'Best MS/MS')]
colnames(ev.f) <- c('Sequence','Proteins','Razor.protein',
                    'Raw.file','Retention.time','PEP', 'Peptide.ID', 'Obs.ID')

# if we don't have regression data for some experiments that have ALL PEPs > the
# threshold of 0.05, then remove them
ev <- ev[(ev$`Raw file` %in% ev.f$Raw.file),]

# get experiment ids
experiment_factors <- as.factor(ev$`Raw file`)
experiment_id <- as.numeric(experiment_factors)
num_exps = max(experiment_id)
# get peptide ids
true_peptide_id <- ev.f$Peptide.ID
pep.id.list <- unique(true_peptide_id)
peptide_id <- as.numeric(as.factor(true_peptide_id))
# get observation IDs - for mapping back to original data format
obs_id <- as.numeric(ev.f$Obs.ID)

# parse linear regression params from STAN output
beta0 <- pars[sprintf('beta_0[%i]', seq(1, num_exps))]
beta1 <- pars[sprintf('beta_1[%i]', seq(1, num_exps))]
beta2 <- pars[sprintf('beta_2[%i]', seq(1, num_exps))]
split.point = pars[sprintf('split_point[%i]', seq(1, num_exps))]
# parse canonical retention times from STAN output
mus <- pars[grep('mu\\[', names(pars))]
# parse RT standard deviations from STAN output
sigmas <- pars[grep('sigma\\[', names(pars))]

# build global false positive density function
ev.rt.edf <- density(ev$`Retention time`)

# output table
ev.new <- data.frame(
  Raw.file=character(),
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
  exp_name <- levels(experiment_factors)[exp_id]
  # get exp subset of ev
  exp <- subset(ev, ev$`Raw file`==exp_name, 
                c('Raw file', 'Sequence', 'PEP', 'Retention time', 
                  'Best MS/MS', 'Peptide ID'))

  cat('\r', i, '/', num_exps, exp_name, nrow(exp), '                           ')
  flush.console()
  
  # some experiments have very few rows??? skip these
  #if(dim(exp)[1] < 10) next
  
  # not all peptides in this experiment have data from the model
  # we can only update those that have that data. others will not be touched
  exp.matches <- exp$`Peptide ID` %in% pep.id.list
  exp.f <- subset(exp, exp.matches)
  
  # convert true_peptide_id to peptide_id
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
  exp.sigmas <- sigmas[exp.peptides]
  
  # PEP.new = P(-|RT) = P(RT|-)*P(-) / (P(RT|-)*P(-) + P(RT|+)*P(+)
  # + <- PSM=Correct
  # - <- PSM=Incorrect
  
  # P(RT|-) = probability of peptides RT, given that PSM is incorrect
  #           calculated from the empirical distribution of all RTs in this experiment
  #exp.rt.edf <- density(exp$`Retention time`)
  exp.rt.minus <- approx(ev.rt.edf$x, ev.rt.edf$y, 
                         xout=exp.f$`Retention time`)$y
  
  # P(-) = probability that PSM is incorrect (PEP)
  # P(+) = probability that PSM is correct (1-PEP)
  
  # P(RT|+) = probability that given the correct ID, the RT falls in the
  #           lognormal distribution of RTs for that peptide, for that experiment
  exp.rt.plus <- unlist(sapply(exp.peptides, function(id) {
    rts <- exp.f$`Retention time`[exp.peptide.map==id]
    muij <- exp.mus[exp.peptides==id]
    sigma <- sigmas[id]
    dnorm(rts, mean=muij, sd=sigma)
  }))
  # sometimes rt.plus will go so low that it will round to 0
  # just round this back up to 1e-100
  exp.rt.plus[exp.rt.plus == 0] = 1e-100
  
  exp.PEP <- exp.f$PEP
  # sometimes MQ will output PEP > 1, which makes no sense, and will
  # result in negative values for our adjusted PEP
  # set all PEP > 1 to PEP = 1
  exp.PEP[exp.PEP > 1] <- 1
  
  # now we can update the PEP
  PEP.new <- (exp.rt.minus*exp.PEP) / 
    ((exp.rt.minus*exp.PEP) + (exp.rt.plus*(1-exp.PEP)))
  
  #hist(PEP.new[PEP.new<0.1], breaks=seq(0, 0.1, by=0.001), main=i)
  
  # output table for this experiment, for updated data
  exp.new <- data.frame(
    Raw.file=as.character(exp.f$`Raw file`),
    Obs.ID=as.numeric(exp.f$`Best MS/MS`),
    rt.minus=as.numeric(exp.rt.minus),
    rt.plus=as.numeric(exp.rt.plus),
    muijs=as.numeric(exp.mus[match(exp.peptide.map, exp.peptides)]),
    sigmas=as.numeric(exp.sigmas[match(exp.peptide.map, exp.peptides)]),
    PEP.new=as.numeric(PEP.new)
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

write.table(ev.adjusted, 'dat/ev.adjusted.elite.txt', sep='\t', row.names=FALSE, quote=FALSE)

ev.ff <- ev.adjusted[, c('Sequence', 'Proteins', 'Leading razor protein', 'Raw file', 
                'Retention time', 'Retention length', 'PIF', 'PEP', 'Intensity',
                'Reporter intensity corrected 0', 'Reporter intensity corrected 1', 
                'Reporter intensity corrected 2', 'Reporter intensity corrected 3', 
                'Reporter intensity corrected 4', 'Reporter intensity corrected 5', 
                'Reporter intensity corrected 6', 'Reporter intensity corrected 7', 
                'Reporter intensity corrected 8', 'Reporter intensity corrected 9',
                'Reverse', 'Peptide ID', 
                'rt.minus', 'rt.plus', 'muijs', 'sigmas', 'PEP.new')]
write.table(ev.ff, 'dat/ev.adj.corr.txt', sep='\t', row.names=FALSE, quote=FALSE)
#write.table(ev.ff, 'dat/ev.adj.elite_sigmas10.txt', sep='\t', row.names=FALSE, quote=FALSE)
