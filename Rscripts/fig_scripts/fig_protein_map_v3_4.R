library(tidyverse)
source('Rscripts/lib.R')

#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180621_2/ev_updated.txt")
ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180724_3/ev_updated.txt")

## add percolator data

source('Rscripts/add_percolator.R')

## run protein quant -------

source('Rscripts/protein_quant.R')

# take first 50 experiments
dmat_cc <- dmat_c[,1:50]

# remove rows of all 0s
dmat_cc <- dmat_cc[!apply(dmat_cc, 1, sum)==0,]

# reorder
odr <- rev(order(apply(dmat_cc==1, 1, sum)))
dmat_cc <- dmat_cc[odr,]

## ----------

# do this with the lower level raster functions

dmat_cc[dmat_cc==0] <- '#FFFFFF'
dmat_cc[dmat_cc==1] <- '#000000'
dmat_cc[dmat_cc==2] <- '#FF0000'
#dmat_cc[dmat_cc==3] <- '#0000FF'
dmat_cc[dmat_cc==3] <- '#000000'

dmat_cc <- as.raster(dmat_cc)

## ---------

#pdf(file='manuscript/Figs/protein_map_v4.pdf', width=3.5, height=6)
png(file='manuscript/Figs/protein_map_v7.png', width=3.5, height=6, units='in', res=250)

par(mar=c(4, 2, 1.25, 0.6),
    cex.axis=0.75, cex.lab=0.75,
    las=1, pty='m')

plot(0, 0, type='n', 
     xlim=c(0, ncol(dmat_cc)), ylim=c(0, nrow(dmat_cc)),
     xaxt='n', yaxt='n', xaxs='i', yaxs='i',
     xlab=NA, ylab=NA)
rasterImage(dmat_cc[nrow(dmat_cc):1,], 0, 0, ncol(dmat_cc), nrow(dmat_cc), interpolate=F)

# draw rectangle around new peptides

# bottom edge of the rect
ry0 <- sum(apply(dmat_c[,1:50], 1, function(x) { any(x == 1)}))
ry1 <- nrow(dmat_cc)
rect(xleft=0, xright=ncol(dmat_cc), ybottom=ry0, ytop=ry1, 
     col=rgb(1,0,0,0.1), border=NA)
#text(15, 6500, 'New Peptide IDs', adj=c(0, 0.5), cex=1, col='red')

axis(1, tck=-0.02, 
     at=seq(0, 50, by=10),
     mgp=c(0, 0.1, 0))
axis(2, tck=-0.02, 
     at=seq(0, 11000, by=1000),
     labels=seq(0, 11, 1),
     mgp=c(0, 0.5, 0))

# legend(x=0, y=-750, xjust=0, yjust=1,
#        c('Not Quantified', 'Spectra', 'DART-ID (Upgraded)', 'DART-ID (Downgraded)'),
#        #col=c('#FFFFFF', '#000000', '#FF0000', '#0000FF'),
#        pch=22, pt.cex=1.5, pt.lwd=1, 
#        pt.bg=c('#FFFFFF', '#000000', '#FF0000', '#0000FF'), 
#        ncol=2, bty='n', cex=0.85,
#        x.intersp=1.1, y.intersp=1, text.width=75, xpd=T)
legend(x=0, y=-500, xjust=0, yjust=1,
       c('Not Quantified', 'Spectra', 'DART-ID (Upgraded)'),
       #col=c('#FFFFFF', '#000000', '#FF0000'),
       pch=22, pt.cex=1.5, pt.lwd=1, 
       pt.bg=c('#FFFFFF', '#000000', '#FF0000'), 
       ncol=2, bty='n', cex=0.85,
       x.intersp=1.1, y.intersp=1, text.width=18, xpd=T)

mtext('Experiment', 1, line=1, cex=1)
mtext('# Distinct Peptides * 1000', 2, line=1, cex=1, las=3)
#mtext('Quantified Peptide Coverage Increase, FDR < 1%', 3, line=0.1, cex=0.75, font=2)
mtext('Peptide Coverage Increase, FDR < 1%', 3, line=0.1, cex=1, font=2)

dev.off()

