library(tidyverse)
source('Rscripts/lib.R')

#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180621_2/ev_updated.txt")
#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180724_3/ev_updated.txt")
#ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180813_with_PI/ev_updated.txt')
ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180815_2/ev_updated.txt')

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

#save(boxs, file='dat/peptide_ids_20180816.rds')

## ---------

load('dat/peptide_ids_20180816.rds')

# horizontal boxplot ------------------------------------------------------

pdf(file='manuscript/Figs/peps_per_exp_v8.pdf', width=1.75, height=1.5)

par(mar=c(1,2.5,0.1,0.25),
    oma=c(0,0.25,1.15,0),
    pty='m', las=1, cex.axis=0.6)

boxplot(rev(boxs), horizontal=T,
        col=rev(c(cb[1], cb[3], paste0(cb[3], '44'), cb[2], paste0(cb[2], '88'))), 
        xaxt='n', yaxt='n', ylim=c(50, 2550),
        outwex=1, outcex=0.75, outpch='x', outcol=rgb(0, 0, 0, 0))
axis(1, at=seq(0, 2500, by=500), labels=c(0, NA, 1000, NA, 2000, NA), tck=-0.02, 
     mgp=c(0, -0.1, 0))
axis(2, at=c(1.5, 3.5, 5), labels=c('DART-ID', 'Percolator', 'Spectra'),
     tck=-0.02, mgp=c(0, 0.3, 0), las=1)

mtext('    Peptides ID\'d per Experiment', side=3, line=0, las=1, font=2, cex=0.7, outer=T)

dev.off()

## ---------

pdf(file='manuscript/Figs/peps_per_exp_v7.pdf', width=1.75, height=1.5)

par(mar=c(2,1.75,0.25,0.25),
    pty='m', las=1, cex.axis=0.65)

# plot(0, 0, type='n',
#      xlim=c(0, 3), ylim=c(0, 1000),
#      xlab=NA, ylab=NA,
#      xaxt='n', yaxt='n')

boxplot(boxs, 
        col=c(av[1], av[3], paste0(av[3], '66'), av[2], paste0(av[2], '66')), 
        xaxt='n', yaxt='n', ylim=c(75, 2525),
        outwex=1, outcex=0.75, outpch='x', outcol=rgb(0, 0, 0, 0.5))
axis(1, at=c(1, 2.5, 4.5), tck=-0.02, labels=NA,
     #labels=c('Spectra', 'Percolator', 'DART-ID'),
     srt=45)
axis(2, at=seq(0, 2500, by=500), labels=seq(0, 2500, by=500), tck=-0.02, 
     mgp=c(0, 0.1, 0), las=3)
text(c(1, 2.5, 4.5), par('usr')[3]-150, srt=45, adj=c(1, 0.5), xpd=T,
     labels=c('Spectra', 'Percolator', 'DART-ID'), cex=0.6)
#text(c())

mtext('Peptides ID\'d per Experiment           ', 2, line=1, las=3, cex=0.7)

dev.off()
