library(tidyverse)
source('Rscripts/lib.R')

ev <- read_tsv("/gd/Slavov_Lab/Albert/RTLib_Alignments/SQC_20180621_2/ev_updated.txt")

## add percolator data

source('Rscripts/add_percolator.R')

## run protein quant -------

source('Rscripts/protein_quant.R')

## ---------
boxs <- list(Spectra=apply(dmat, 2, sum), 
             Percolator=apply(dmat_perc, 2, sum),
             DART=apply(dmat_new, 2, sum))
## ---------

pdf(file='manuscript/Figs/peps_per_exp.pdf', width=1.5, height=2.75)

par(mar=c(2.5,1.75,1.25,0.5),
    pty='m', las=1, cex.axis=0.75)

# plot(0, 0, type='n',
#      xlim=c(0, 3), ylim=c(0, 1000),
#      xlab=NA, ylab=NA,
#      xaxt='n', yaxt='n')

boxplot(boxs, col=c(av[1], av[3], av[2]), 
        xaxt='n', yaxt='n',
        outwex=1, outcex=0.75, outpch='x', outcol=rgb(0, 0, 0, 0.5))
axis(1, at=c(1, 2, 3), tck=-0.02, labels=NA,
     #labels=c('Spectra', 'Percolator', 'DART-ID'),
     srt=45)
axis(2, at=seq(0, 2000, 500), labels=seq(0, 2000, 500), tck=-0.02, 
     mgp=c(0, 0.1, 0), las=3)
text(c(1, 2, 3), -125, srt=45, adj=c(1, 0.5), xpd=T,
     labels=c('Spectra', 'Percolator', 'DART-ID'), cex=0.7)

mtext('Peptides Quantified per Experiment', 2, line=1, las=3, cex=0.75)

dev.off()
