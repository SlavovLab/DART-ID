# init
library(tidyverse)
library(rmutil)

#load('dat/params.Fit3.RData')
#load('dat/params.Fit3b.RData')
load('dat/params.Fit3c.RData')

ev <- read_tsv('dat/evidence_elite.txt')

# alias the modified peptide ID as the peptide ID
ev <- ev %>%
  rename(`Sequence ID`=`Peptide ID`) %>%
  rename(`Peptide ID`=`Mod. peptide ID`) 

## Filter of PEP < pep_thresh
## !!!!!!!
## make sure the pep_thresh is the same that the parameters were built with
## !!!!!!!
pep_thresh <- 0.5
ev.f <- ev %>% 
  filter(PEP < pep_thresh) %>%
  filter(!grepl('REV*', `Leading razor protein`)) %>% # Remove Reverse matches
  filter(!grepl('CON*',`Leading razor protein`))  %>% # Remove Contaminants
  select("Peptide ID", "Raw file", "Retention time", "PEP") %>%
  mutate(exp_id=`Raw file`) %>%  # new column - exp_id = numeric version of experiment file
  mutate_at("exp_id", funs(as.numeric(as.factor(.))))

## Remove peptides that occur in less than n experiments
toRemove <- ev.f %>% 
  group_by(`Peptide ID`) %>% 
  summarise(n=length(unique(`exp_id`))) %>% 
  filter(n <= 10) %>%
  pull(`Peptide ID`)
ev.f <- ev.f %>% 
  filter(!(`Peptide ID` %in% toRemove))

experiment_factors <- as.factor(ev.f$`Raw file`)
experiment_ids <- ev.f$exp_id
num_exps <- length(unique(ev.f[["exp_id"]]))

## "true peptide id" matches peptide id in evidence file
## "peptide id" is index in 1:num_peptides for stan
raw_peptide_id <- ev.f$`Peptide ID`
pep.id.list <- unique(raw_peptide_id)
stan_peptide_id <- as.numeric(as.factor(ev.f[["Peptide ID"]]))
ev.f[["Stan ID"]] <- stan_peptide_id

num_total_observations <- nrow(ev.f)
num_peptides <- length(unique(stan_peptide_id))

retention_times <- ev.f[["Retention time"]]
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
  Obs.ID=numeric(),
  rt.minus=numeric(),
  rt.plus=numeric(),
  muijs=numeric(),
  sigmas=numeric(),
  PEP.new=numeric()
)

print('Processing...')
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
  
  # get muijs for this experiment
  #exp.f$mus <- muij_fit[exp.peptide.map]
  # get sigma_ijs for this experiment
  #exp.f$sigmas <- sigma_ijs[exp.peptide.map]
  
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
  
  # global - Fit3
  #exp.rt.minus <- 1 / max(retention_times) 
  
  # Fit3b+c, lognormal density over all retention times
  exp.rt.minus <- dlnorm(exp.f$`Retention time`, meanlog=mean.log.rt, sdlog=sd.log.rt)
  
  # P(RT|-) = probability of peptides RT, given that PSM is incorrect
  #           calculated from the empirical distribution of all RTs in this experiment
  #exp.rt.edf <- density(exp$`Retention time`, adjust=2) # smooth it out a bit
  #exp.rt.minus <- approx(exp.rt.edf$x, exp.rt.edf$y, 
  #                       xout=exp.f$`Retention time`)$y
  
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
    
    # Fit3 - uniform density + laplace
    #((pep) * dunif(rts, min=0, max=max.rt)) + 
    # ((1-pep) * dlaplace(rts, m=muij, s=sigma))
    
    # Fit3b - lognormal density + laplace
    #((pep) * dlnorm(rts, meanlog=mean.log.rt, sdlog=sd.log.rt)) +
    #((1-pep) * dlaplace(rts, m=muij, s=sigma))
    
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

#write.table(ev.adjusted, 'dat/ev.adjusted.elite.txt', sep='\t', row.names=FALSE, quote=FALSE)

ev.ff <- ev.adjusted[, c('Sequence', 'Proteins', 'Leading razor protein', 'Raw file', 
                         'Retention time', 'Retention length', 'PIF', 'PEP', 'Intensity',
                         'Reporter intensity corrected 0', 'Reporter intensity corrected 1', 
                         'Reporter intensity corrected 2', 'Reporter intensity corrected 3', 
                         'Reporter intensity corrected 4', 'Reporter intensity corrected 5', 
                         'Reporter intensity corrected 6', 'Reporter intensity corrected 7', 
                         'Reporter intensity corrected 8', 'Reporter intensity corrected 9',
                         'Reverse', 'Peptide ID', 'Sequence ID', 'Modifications',
                         'rt.minus', 'rt.plus', 'muijs', 'sigmas', 'PEP.new', 'id')]

#write.table(ev.ff, 'dat/ev.adj.Fit3.txt', sep='\t', row.names=FALSE, quote=FALSE)
write_tsv(ev.ff, 'dat/ev.adj.Fit3c.txt')



## analyze the results -----

ev.fit2 <- parse.ev.adj('dat/ev.adj.Fit2.txt')
ev.fit2 <- read_tsv('dat/ev.adj.Fit2.txt')

#ev.a <- ev.fit2 %>%
ev.a <- ev.ff %>%
  mutate(exp_id=`Raw file`) %>%  # new column - exp_id = numeric version of experiment file
  mutate_at("exp_id", funs(as.numeric(as.factor(.)))) %>%
  group_by(`Peptide ID`, `exp_id`) %>%
  select('Peptide ID', 'exp_id', 'rt.minus', 'rt.plus', 'Retention time', 
         'muijs', 'sigmas', 'PEP', 'PEP.new', 'id')

i <- sample.int(n=nrow(ev.a), size=1e5)

ev.a[i,] %>%
  filter(!is.na(PEP.new)) %>%
  #filter(`exp_id`==1) %>%
  ggplot(aes(x=PEP, y=PEP.new)) +
    #geom_point(alpha=0.1)+
    stat_density2d(aes(fill=..level..), geom='polygon', n=50) +
    geom_abline(intercept=0, slope=1, color='red') +
    scale_x_log10(limits=c(1e-8, 1)) + scale_y_log10(limits=c(1e-8, 1)) +
    labs(title='RTLib: Fit2 (updated sigmas)')
