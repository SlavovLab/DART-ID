# analyze contaminated PSMs ---------

# init -------
library('readr')

# load msms.txt -----
path.data <- '/Volumes/salamanca/GoogleDrive/SingleCell_Data/'
ev <- read_tsv(paste(c(path.data, 'All_human/MQ_exp_37-51/msms.txt'), collapse=''),
                 col_types=cols_only(
                   `Raw file`=col_character(),
                   Sequence=col_character(),
                   Proteins=col_character(),
                   PEP=col_number(),
                   `Retention time`=col_number(),
                   PIF=col_number(),
                   `Precursor Intensity`=col_number()
                 ))
names(ev)[names(ev) == 'Raw file'] <- 'Raw.file'
# psm/pep ------- 
# Number of PSMs per contaminant peptide and per non-contaminant peptide

contam.inds <- grep('CON_', ev$Proteins)

contam.ev <- ev[contam.inds,]
uncontam.ev <- ev[-contam.inds,]

# contaminants:
contam.freq <- as.data.frame(table(contam.ev$Sequence))
with(contam.freq, hist(Freq[Freq < 50], 
    breaks=c(seq(0,50,by=5)), 
    main=paste('PSMs/Contaminated Peptide\n', '50+ PSMs/pep:', sum(Freq>50)),
    xlab='PSMs/Peptide'
))

# non contaminants
uncontam.freq <- as.data.frame(table(uncontam.ev$Sequence))
with(uncontam.freq, hist(Freq[Freq < 50], 
    breaks=c(seq(0,50,by=5)), 
    main=paste('PSMs/Uncontaminated Peptide\n', '50+ PSMs/pep:', sum(Freq>50)),
    xlab='PSMs/Peptide'
))

# dist ------
# A distribution of fraction of PSMs spent on contaminants. That is,
# for each file compute the fraction of PSMs mapped to a 
# CON__ peptide and plot these fractions as distributions 
# for at several PEP thresholds, i.e., 0.01, 0.05, 0.1, 1

exps <- unique(ev$Raw.file)
pep.thresh <- c(0.005, 0.01, 0.05, 0.1, 1)
contam.frac <- matrix(data=NA, nrow=length(exps), ncol=length(pep.thresh))
counter <- 0
for (i in exps) {
  counter <- counter + 1
  new.exp <- subset(ev, ev$Raw.file==i)
  
  counter.j <- 0
  # for each PEP thresh
  for (j in pep.thresh) {
    counter.j <- counter.j + 1  
    # apply PEP threshold
    new.exp.pep = subset(new.exp, new.exp$PEP < j)
    # compute contaminated fraction
    contam.frac[counter,counter.j] <- 
      length((grep('CON_', new.exp.pep$Proteins))) / nrow(new.exp.pep)
  }
}


