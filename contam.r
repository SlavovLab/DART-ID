# analyze contaminated PSMs ---------

# init -------
library('readr')
library('ggplot2')
library('reshape2')

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
names(ev)[names(ev) == 'Retention time'] <- 'Retention.time'
# psm/pep ------- 
# Number of PSMs per contaminant peptide and per non-contaminant peptide

contam.inds <- grep('CON_', ev$Proteins)

contam.ev <- ev[contam.inds,]
uncontam.ev <- ev[-contam.inds,]

# contaminants:
contam.freq <- as.data.frame(table(contam.ev$Sequence))
with(contam.freq, hist(Freq[Freq < 50], 
    breaks=c(seq(0,50,by=1)), 
    main=paste('PSMs/Contaminated Peptide\n', '50+ PSMs/pep:', sum(Freq>50)),
    xlab='PSMs/Peptide'
))

# non contaminants
uncontam.freq <- as.data.frame(table(uncontam.ev$Sequence))
with(uncontam.freq, hist(Freq[Freq < 50], 
    breaks=c(seq(0,50,by=1)), 
    main=paste('PSMs/Uncontaminated Peptide\n', '50+ PSMs/pep:', sum(Freq>50)),
    xlab='PSMs/Peptide'
))

# get list of peptides/proteins with over 30 PSMs
pep.freq <- as.data.frame(table(ev$Sequence))
red.peps.list <- as.character(pep.freq[pep.freq$Freq>30, 'Var1'])
ev.f = subset(ev, ev$Sequence %in% red.peps.list)

row.names(ev.f) <- NULL
write.table(ev.f, '~/git/RTLib/dat/psm30.txt', sep = '\t', 
            row.names = FALSE, quote = FALSE)

# condense this data
ev.f2 <- ev.f[match(red.peps.list, ev.f$Sequence), c('Sequence', 'Proteins')] 
ev.f2 <- cbind(ev.f2, as.data.frame(table(ev.f$Sequence))$Freq )
names(ev.f2)[dim(ev.f2)[2]] <- 'PSMs'

row.names(ev.f2) <- NULL
write.table(ev.f2, '~/git/RTLib/dat/psm30_condensed.txt', sep = '\t', 
            row.names = FALSE, quote = FALSE)


# the variance for the RTs of contaminant and non-contaminant peptides. ------

exps = unique(ev$Raw.file)
all.rt.var = data.frame(
  Raw.file=character(),
  Sequence=character(),
  RTVar=numeric(),
  Contam=logical(),
  stringsAsFactors=FALSE
)
c <- 0
for (i in exps) {
  c <- c + 1
  new.exp = subset(ev, ev$Raw.file==i)
  new.exp.peps = unique(new.exp$Sequence)
  
  rt.var = data.frame(
    Raw.file=as.character(rep(i, length(new.exp.peps))),
    Sequence=as.character(new.exp.peps),
    RTVar=as.numeric(rep(NA, length(new.exp.peps))),
    Contam=as.logical(rep(FALSE, length(new.exp.peps))),
    stringsAsFactors=FALSE
  )
  pep.inds <- match(new.exp.peps, new.exp$Sequence)
  rt.var$Contam[grep('CON_', new.exp[pep.inds,c('Proteins')]$Proteins)] <- TRUE
  rt.var$RTVar <- aggregate(Retention.time ~ Sequence, data=new.exp, FUN=var)$Retention.time
  
  all.rt.var <- rbind(all.rt.var, rt.var)
  
  cat('\r', 'Processing:', c, '/', length(exps), i, '                         ');
  flush.console()
}

row.names(all.rt.var) <- NULL
write.table(all.rt.var, '~/git/RTLib/dat/all.rt.var.txt', sep = '\t', 
            row.names = FALSE, quote = FALSE)

## plot contam/uncontam variances ------

# experiment to use
exp <- exps[10]

new.exp = subset(all.rt.var, ev$Raw.file==exp)
par(mfrow=c(2,1)) 
with(new.exp, hist(RTVar[Contam==TRUE & RTVar < 1], 
 breaks=seq(0,1,by=0.01),
 main=paste('RT Variance for Contaminated Peptides\n', 
            'Peps w/ Variance > 1:', sum(Contam==TRUE & RTVar > 1, na.rm=TRUE),
            '(', formatC(sum(Contam==TRUE & RTVar > 1, na.rm=TRUE) / 
                           sum(Contam==TRUE, na.rm=TRUE) * 100, digits=4), '%)',
            '\n', exp),
 xlab='Retention Time Variance'
))
with(new.exp, hist(RTVar[Contam==FALSE & RTVar < 1], 
 breaks=seq(0,1,by=0.01),
 main=paste('RT Variance for Non-Contaminated Peptides\n', 
            'Peps w/ Variance > 1:', sum(Contam==FALSE & RTVar > 1, na.rm=TRUE),
            '(', formatC(sum(Contam==FALSE & RTVar > 1, na.rm=TRUE) / 
                           sum(Contam==FALSE, na.rm=TRUE) * 100, digits=4), '%)',
            '\n', exp),
 xlab='Retention Time Variance'
))


# dist ------
# A distribution of fraction of PSMs spent on contaminants. That is,
# for each file compute the fraction of PSMs mapped to a 
# CON__ peptide and plot these fractions as distributions 
# for at several PEP thresholds, i.e., 0.01, 0.05, 0.1, 1

exps <- unique(ev$Raw.file)
pep.thresh <- c(0.005, 0.01, 0.05, 0.1, 1)
contam.frac <- as.data.frame(matrix(data=NA, nrow=length(exps), ncol=length(pep.thresh)))
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
colnames(contam.frac) <- paste('PEP <', as.character(pep.thresh))

ggplot(melt(contam.frac), aes(x=value, fill=variable)) + 
  geom_histogram(binwidth=0.05) +
  facet_grid(variable~.) + 
  labs(x='Contaminated PSM Fraction', title='Contaminated PSM Fractions')

