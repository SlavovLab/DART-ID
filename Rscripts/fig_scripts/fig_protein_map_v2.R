library(tidyverse)
source('Rscripts/lib.R')

ev <- read_tsv("/gd/Slavov_Lab/Albert/RTLib_Alignments/SQC_20180621_2/ev_updated.txt")

## add percolator data

source('Rscripts/add_percolator.R')

## run protein quant -------

source('Rscripts/protein_quant.R')

## ---------

pdf(file='manuscript/Figs/protein_map_v2.pdf', width=3.5, height=2.75)

par(mar=c(1.75, 1.75, 1.25, 0.25),
    cex.axis=0.75, cex.lab=0.75,
    las=1, pty='m')

image(t(dmat_c), useRaster=T, col=c("#FFFFFF", "#000000", "#FF0000", "#0000FF"),
      xlab=NA, ylab=NA,
      xaxt='n', yaxt='n')

axis(1, tck=-0.02, 
     at=seq(0, 180, by=30) / ncol(dmat_c),
     labels=seq(0, 180, by=30),
     mgp=c(0, 0.1, 0))
axis(2, tck=-0.02, 
     at=seq(0, 8000, by=1000) / nrow(dmat_c),
     labels=seq(0, 8, 1),
     mgp=c(0, 0.5, 0))

mtext('Experiment', 1, line=0.8, cex=0.75)
mtext('Peptide Sequence * 1000', 2, line=1, cex=0.75, las=3)
mtext('Peptide Coverage Increase', 3, line=0.1, cex=0.75, font=2)

dev.off()

## ----------

# do this with the lower level raster functions

dmat_cc <- dmat_c
dmat_cc[dmat_c==0] <- '#FFFFFF'
dmat_cc[dmat_c==1] <- '#000000'
dmat_cc[dmat_c==2] <- '#FF0000'
dmat_cc[dmat_c==3] <- '#0000FF'

dmat_cc <- as.raster(dmat_cc)

## ---------

pdf(file='manuscript/Figs/protein_map_v3.pdf', width=3.5, height=2.75)

par(mar=c(3.25, 1.75, 1.25, 0.25),
    cex.axis=0.75, cex.lab=0.75,
    las=1, pty='m')

plot(0, 0, type='n', 
     xlim=c(0, ncol(dmat_c)), ylim=c(0, nrow(dmat_c)),
     xaxt='n', yaxt='n', xaxs='i', yaxs='i',
     xlab=NA, ylab=NA)
rasterImage(dmat_cc, 0, 0, ncol(dmat_c), nrow(dmat_c), interpolate=F)

axis(1, tck=-0.02, 
     at=seq(0, 180, by=30),
     labels=seq(0, 180, by=30),
     mgp=c(0, 0.05, 0))
axis(2, tck=-0.02, 
     at=seq(0, 8000, by=1000),
     labels=seq(0, 8, 1),
     mgp=c(0, 0.5, 0))

legend(x=0, y=-1250, xjust=0, yjust=1,
       c('Not Quantified', 'Spectra', 'DART-ID (Upgraded)', 'DART-ID (Downgraded)'),
       #col=c('#FFFFFF', '#000000', '#FF0000', '#0000FF'),
       pch=22, pt.cex=1.5, pt.lwd=1, 
       pt.bg=c('#FFFFFF', '#000000', '#FF0000', '#0000FF'), 
       ncol=2, bty='n', cex=0.7,
       x.intersp=0.75, y.intersp=1, text.width=50, xpd=T)

mtext('Experiment', 1, line=0.7, cex=0.75)
mtext('Peptide Sequence * 1000', 2, line=1, cex=0.75, las=3)
mtext('Quantified Peptide Coverage Increase, FDR < 1%', 3, line=0.1, cex=0.75, font=2)

dev.off()

