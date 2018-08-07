library(tidyverse)
source('Rscripts/lib.R')

#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180621_2/ev_updated.txt")
ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180724_3/ev_updated.txt")

## add percolator data

source('Rscripts/add_percolator.R')

## run protein quant -------

source('Rscripts/protein_quant.R')

## ---------
no_conf_ids <- apply(dmat, 1, sum) == 0
boxs <- list(Spectra=apply(dmat, 2, sum), 
             Percolator=apply(dmat_perc, 2, sum),
             Percolator_prev=apply(dmat_perc[!no_conf_ids,], 2, sum),
             DART=apply(dmat_new, 2, sum),
             DART_prev=apply(dmat_new[!no_conf_ids,], 2, sum))
## ---------

pdf(file='manuscript/Figs/peps_per_exp_v6.pdf', width=1.25, height=3)

par(mar=c(2.5,2,0.5,0.5),
    pty='m', las=1, cex.axis=0.75)

# plot(0, 0, type='n',
#      xlim=c(0, 3), ylim=c(0, 1000),
#      xlab=NA, ylab=NA,
#      xaxt='n', yaxt='n')

boxplot(boxs, 
        col=c(av[1], av[3], paste0(av[3], '66'), av[2], paste0(av[2], '66')), 
        xaxt='n', yaxt='n', ylim=c(0, 2600),
        outwex=1, outcex=0.75, outpch='x', outcol=rgb(0, 0, 0, 0.5))
axis(1, at=c(1, 2.5, 4.5), tck=-0.02, labels=NA,
     #labels=c('Spectra', 'Percolator', 'DART-ID'),
     srt=45)
axis(2, at=seq(0, 2500, by=500), labels=seq(0, 2500, by=500), tck=-0.02, 
     mgp=c(0, 0.1, 0), las=3)
text(c(1, 2.5, 4.5), par('usr')[3]-75, srt=45, adj=c(1, 0.5), xpd=T,
     labels=c('Spectra', 'Percolator', 'DART-ID'), cex=0.7)
#text(c())

mtext('Peptides ID\'d per Experiment   ', 2, line=0.95, las=3, cex=1)

dev.off()
