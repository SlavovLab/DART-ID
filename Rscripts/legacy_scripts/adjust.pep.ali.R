# adjust.pep.ali.R
# experiment-centric method of adjusting PEPs, using STAN parameters
# i.e., compare distributions of correct vs. incorrect IDs

adjust.pep.ali <- function(path.in='dat/evidence.txt',
                           path.out='dat/ev+dRT.elite.txt') {

library(readr)

## load STAN params
load('dat/params_elite.RData')
#load('dat/params_all.RData')

## load evidence
ev <- read_tsv(path.in)
# only keep experiments done on Elite system
ev <- ev[grep('[0-9]{6}A', ev$`Raw file`),]

ev <- ev[, c('Raw file', 'Sequence', 'PEP', 'Retention time', "PIF", 
             "Reporter intensity corrected 0", "Reporter intensity corrected 1", 
             "Reporter intensity corrected 2", "Reporter intensity corrected 3",
             "Reporter intensity corrected 4", "Reporter intensity corrected 5", 
             "Reporter intensity corrected 6", "Reporter intensity corrected 7", 
             "Reporter intensity corrected 8", "Reporter intensity corrected 9",
             "id", "Leading razor protein", 'Best MS/MS', 'Peptide ID')]
ev.f <- ev[ev$PEP < 0.05,]

# if we don't have regression data for some experiments that have ALL PEPs > the
# threshold of 0.05, then remove them
ev <- ev[(ev$`Raw file` %in% ev.f$`Raw file`),]
ev$prot.label <- substr(ev$`Leading razor protein`, 1, 3)

# get experiment ids
experiment_factors <- as.factor(ev$`Raw file`)
experiment_id <- as.numeric(experiment_factors)
num_exps = max(experiment_id)
# get peptide ids
true_peptide_id <- ev.f$`Peptide ID`
pep.id.list <- unique(true_peptide_id)
peptide_id <- as.numeric(as.factor(true_peptide_id))
# get observation IDs - for mapping back to original data format
obs_id <- as.numeric(ev.f$`Best MS/MS`)

# parse linear regression params from STAN output
beta0 <- pars[sprintf('beta_0[%i]', seq(1, num_exps))]
beta1 <- pars[sprintf('beta_1[%i]', seq(1, num_exps))]
beta2 <- pars[sprintf('beta_2[%i]', seq(1, num_exps))]
split.point = pars[sprintf('split_point[%i]', seq(1, num_exps))]
# parse canonical retention times from STAN output
mus <- pars[grep('mu\\[', names(pars))]
# parse RT standard deviations from STAN output
sigmas <- pars[grep('sigma\\[', names(pars))]

# convert peptide ID to STAN peptide id
ev$ID <- match(ev$`Peptide ID`, pep.id.list)
# attach the adjusted median RT to these peptides
ev$mu <- as.vector(mus[sprintf('mu[%d]', ev$ID)])

## adjust original RT using linear regression params
ev.exp.map <- match(ev$`Raw file`, levels(experiment_factors))
ev$beta0 <- beta0[sprintf('beta_0[%d]', ev.exp.map)]
ev$beta1 <- beta1[sprintf('beta_1[%d]', ev.exp.map)]
ev$beta2 <- beta2[sprintf('beta_2[%d]', ev.exp.map)]
ev$split.point <- split.point[sprintf('split_point[%d]', ev.exp.map)]

# apply segmented linear regression parameters for this experiment
ev$RT.corrected <- NA
# 2 cases - whether it lies on the first line or the second
ev.split <- (ev$mu - ev$split.point) < 0
ev.split[is.na(ev.split)] <- FALSE
# 1st line
ev$RT.corrected[ev.split] <- 
  ev$beta0[ev.split] +
  (ev$beta1[ev.split] * ev$mu[ev.split])

ev.split <- (ev$mu - ev$split.point) >= 0
ev.split[is.na(ev.split)] <- FALSE
# 2nd line
ev$RT.corrected[ev.split] <-
  ev$beta0[ev.split] +
  (ev$beta1[ev.split] * ev$split.point[ev.split]) +
  (ev$beta2[ev.split] * 
     (ev$mu[ev.split] - ev$split.point[ev.split]))
# difference in adjusted canonical vs. actual retention time
ev$dRT <- abs(ev$RT.corrected-ev$`Retention time`)

# Updating PEPs
ex.out <- data.frame()
exps <- unique(ev.f$`Raw file`)
counter <- 0
for (i in exps) {
  counter <- counter + 1
  
  ex <- subset(ev, ev$`Raw file`==i & ev$PEP <= 1)
  row.names(ex) <- NULL
  
  # incorrect IDs
  rev <- subset(ex, ex$prot.label=='REV')
  # correct IDs
  forw <- subset(ex, ex$prot.label=='sp|' & ex$PEP < 0.02)
  
  # build density functions. ignore missing data
  den.for <- density(forw$dRT, na.rm=TRUE)
  den.rev <- density(rev$dRT, na.rm=TRUE)
  
  # TRUE - belongs to distribution of correct IDs
  ex$Tr <- approx(den.for$x, den.for$y, xout=ex$dRT)$y
  # FALSE - belongs to distribution of incorrect IDs
  ex$Fa <- approx(den.rev$x, den.rev$y, xout=ex$dRT)$y
  # fix missing data - TRUE distribution can be very sparse sometimes
  ex[is.na(ex$Tr) & !is.na(ex$dRT), "Tr"] <- .Machine$double.xmin
  # Bayes Theorem
  ex$PEP.new <- (ex$PEP*ex$Fa)/((ex$PEP*ex$Fa)+((1-ex$PEP)*ex$Tr))
  
  ex.out <- rbind(ex.out, ex)
  
  cat('\r', 'Processing ', counter, '/', length(exps), ' ', i, nrow(forw),
      '                                          ')
  flush.console()
}

ex.out$PEP.updated <- ex.out$PEP.new
ex.out$PEP.updated[is.na(ex.out$PEP.updated)] <-
  ex.out$PEP[is.na(ex.out$PEP.updated)]

row.names(ex.out) <- NULL
#ex.out$Tr <- NULL
#ex.out$Fa <- NULL
#ex.out$dRT <- NULL

if(is.null(path.out)) {
  #source('parse.ev.adj.R')
  #return(parse.ev.adj(ex.out, type='Ali'));
  return(ex.out)
}

write.table(ex.out, path.out, sep = '\t', row.names = FALSE, quote = FALSE)

}


par(mfrow=c(2,2))
for(i in sample.int(num_exps, size=4)) {
  ex <- subset(ev, ev$`Raw file`==exps[i] & ev$PEP <= 1)
  
  # incorrect IDs
  rev <- subset(ex, ex$prot.label=='REV')
  # correct IDs
  forw <- subset(ex, ex$prot.label=='sp|' & ex$PEP < 0.02)
  
  # build density functions. ignore missing data
  den.for <- density(forw$dRT, na.rm=TRUE)
  den.rev <- density(rev$dRT, na.rm=TRUE)
  
  plot(den.for, col='blue', main='Correct (blue) V Incorrect (red)', xlab='dRT (min)')
  lines(den.rev, col='red')
}
