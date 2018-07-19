library(tidyverse)
source('Rscripts/lib.R')

ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180621_2/ev_updated.txt")

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

# only take first n experiments
dmat_cc <- dmat_cc[,1:50]
# remove rows of all 0s
dmat_cc <- dmat_cc[!apply(dmat_cc, 1, sum) == 0,]

# reorder the peptides/proteins so that the 1s are stacked on the bottom
#odr <- rev(order(apply(dmat, 1, sum) + apply(dmat_new, 1, sum)))
odr <- rev(order(apply(dmat_cc == 1 | dmat_cc == 3, 1, sum)))
dmat_cc <- dmat_cc[odr,]

dmat_cc[dmat_cc==0] <- '#FFFFFF'
dmat_cc[dmat_cc==1] <- '#000000'
dmat_cc[dmat_cc==2] <- '#FF0000'
#dmat_cc[dmat_cc==3] <- '#0000FF'
dmat_cc[dmat_cc==3] <- '#000000'

## ---------

#pdf(file='manuscript/Figs/protein_map_v4.pdf', width=3.5, height=6)
png(file='manuscript/Figs/protein_map_v6.png', width=3.5, height=6, units='in', res=250)

par(mar=c(3.75, 2, 1.25, 0.5),
    cex.axis=0.75, cex.lab=0.85,
    las=1, pty='m')

plot(0, 0, type='n', 
     xlim=c(0, ncol(dmat_cc)), ylim=c(0, nrow(dmat_cc)),
     xaxt='n', yaxt='n', xaxs='i', yaxs='i',
     xlab=NA, ylab=NA)
rasterImage(as.raster(dmat_cc[nrow(dmat_cc):1,]), 0, 0, ncol(dmat_cc), nrow(dmat_cc), interpolate=F)

# draw rect around new peptide IDs
num_new_ids <- nrow(dmat_cc) - sum(apply(dmat_cc, 1, function(x) { any(x == '#000000') }))
rect(xleft=0, xright=ncol(dmat_cc), ybottom=nrow(dmat_cc)-num_new_ids, ytop=nrow(dmat_cc),
     border=NA, col=rgb(1, 0, 0, 0.1))
#text(15, 13500, 'New Peptide IDs', adj=c(0, 0.5), col='red', cex=1)

axis(1, tck=-0.02, 
     at=seq(0, ncol(dmat_cc), by=10),
     #labels=seq(0, 240, by=30),
     mgp=c(0, 0.1, 0))
axis(2, tck=-0.02, 
     at=seq(0, nrow(dmat_cc)+1000, by=1000),
     labels=seq(0, nrow(dmat_cc)+1000, by=1000) / 1000,
     #labels=seq(0, 26, 2),
     mgp=c(0, 0.4, 0))

# legend(x=0, y=-550, xjust=0, yjust=1,
#        c('Not Identified', 'Spectra', 'DART-ID (Upgraded)', 'DART-ID (Downgraded)'),
#        #col=c('#FFFFFF', '#000000', '#FF0000', '#0000FF'),
#        pch=22, pt.cex=1.5, pt.lwd=1, 
#        pt.bg=c('#FFFFFF', '#000000', '#FF0000', '#0000FF'), 
#        ncol=2, bty='n', cex=0.7,
#        x.intersp=1.1, y.intersp=1, text.width=13, xpd=T)
legend(x=0, y=-550, xjust=0, yjust=1,
       c('Not Identified', 'Spectra', 'DART-ID (Upgraded)'),
       #col=c('#FFFFFF', '#000000', '#FF0000', '#0000FF'),
       pch=22, pt.cex=1.5, pt.lwd=1, 
       pt.bg=c('#FFFFFF', '#000000', '#FF0000'), 
       ncol=2, bty='n', cex=0.7,
       x.intersp=1.1, y.intersp=1, text.width=13, xpd=T)


mtext('Experiment', 1, line=0.9, cex=1)
mtext('# Distinct Peptides * 1000', 2, line=1.15, cex=1, las=3)
mtext('Peptide Coverage Increase, FDR < 1%', 3, line=0.1, cex=1, font=2)

dev.off()

