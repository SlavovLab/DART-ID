## compare RT variance between multiple alignment with STAN, 
# and reference experiment method

library('readr')
library('MASS')

# reference experiment method -------

fix.evidence.v2('dat/evidence.txt', 'dat/ev.ed.txt')
find.ref.exp('dat/ev.ed.txt', 'dat/ref.picks.txt')
ref.exp <- "160712A_NC_19B_targeted"
make.lib('dat/ev.ed.txt', 'dat/RT.lib.txt','dat/shift.coeffs.txt', ref.exp)
dRT('dat/RT.lib.txt','dat/shift.coeffs.txt','dat/evidence.txt','dat/evidence+dRT.txt')

# reference experiment method - elite only -------

fix.evidence.elite('dat/evidence.txt', 'dat/ev.ed.elite.txt')
find.ref.exp('dat/ev.ed.elite.txt', 'dat/ref.picks.elite.txt')
ref.exp.elite <- "160712A_NC_19B_targeted"
make.lib('dat/ev.ed.elite.txt', 'dat/RT.lib.elite.txt',
         'dat/shift.coeffs.elite.txt', ref.exp.elite)
dRT('dat/RT.lib.elite.txt','dat/shift.coeffs.elite.txt',
    'dat/evidence.txt','dat/evidence+dRT.elite.txt')

## analyze ref experiment method ------- 

ev <- read_tsv('dat/evidence+dRT.txt')
rt.var <- aggregate(RT.corrected ~ Sequence, data=ev, FUN=var)
with(rt.var, hist(RT.corrected[RT.corrected<10], 
                  breaks=seq(0,10,by=0.1),
                  main=paste('dRT Ref Exp Method - RT Variance\n',
                             '# Peptides > 10:',sum(RT.corrected>10, na.rm=TRUE),
                             '(',formatC(sum(RT.corrected>10, na.rm=TRUE)/
                                           length(RT.corrected)*100, digits=4),'% )'),
                  xlab='RT Variance'
))

# elite
ev.elite <- read_tsv('dat/evidence+dRT.elite.txt')
rt.var.elite <- aggregate(RT.corrected ~ Sequence, data=ev.elite, FUN=var)
with(rt.var.elite, hist(RT.corrected[RT.corrected<10], 
                  breaks=seq(0,10,by=0.1),
                  main=paste('dRT Ref Exp Method - RT Variance - Elite Only\n',
                             '# Peptides > 10:',sum(RT.corrected>10, na.rm=TRUE),
                             '(',formatC(sum(RT.corrected>10, na.rm=TRUE)/
                                           length(RT.corrected)*100, digits=4),'% )'),
                  xlab='RT Variance'
))

## analyze multiple alignment method -----

# load evidence, and extract peptide_ids and experiment_ids
evidence <- read_tsv("dat/evidence.txt")
subEvidence <- evidence[evidence$PEP < 0.05, c("Peptide ID", "Raw file", "Retention time")]

true_peptide_id <- subEvidence[["Peptide ID"]]
peptide_id <- as.numeric(as.factor(subEvidence[["Peptide ID"]]))

peptideCounts <- as.data.frame(sort(table(peptide_id), decreasing=TRUE))

# load STAN params
params <- load('dat/params.RData')

# load standard deviations
ma.rt.var <- pars[grep('sigma', names(pars))]
hist(as.numeric(ma.rt.var[ma.rt.var < 10]), 
     breaks=seq(0,10,by=0.1),
     main=paste('Multiple Alignment (STAN) - RT Variance\n',
                '# Peptides > 10:', sum(ma.rt.var>10), '(',
                formatC(sum(ma.rt.var>10)/length(ma.rt.var)*100, digits=4),'% )'),
     xlab='RT Variance')

# STAN w/ Elite only
subEvidence <- evidence[evidence$PEP < 0.05, c("Peptide ID", "Raw file", "Retention time")]
subEvidence <- subEvidence[grep('[0-9]{6}?A', subEvidence$`Raw file`),]

params <- load('dat/params_elite.RData')
ma.rt.var.elite <- pars[grep('sigma', names(pars))]
hist(as.numeric(ma.rt.var.elite[ma.rt.var.elite < 10]), 
     breaks=seq(0,10,by=0.1),
     main=paste('Multiple Alignment (STAN) - RT Variance\nElite Only\n',
                '# Peptides > 10:', sum(ma.rt.var.elite>10), '(',
                formatC(sum(ma.rt.var.elite>10)/length(ma.rt.var.elite)*100, digits=4),'% )'),
     xlab='RT Variance')
