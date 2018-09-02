library(tidyverse)
source('Rscripts/lib.R')

#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180621_2/ev_updated.txt")
#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180724_3/ev_updated.txt")
#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180808_5exp_parametric/ev_updated.txt")
#ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC2_20180812_1/ev_updated.txt')
#ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180813_with_PI/ev_updated.txt')
ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180815_2/ev_updated.txt')

## add percolator data

source('Rscripts/add_percolator.R')

## run protein quant -------

source('Rscripts/protein_quant.R')

# reorder and raster ------------------------------------------------------

# take first 50 experiments
dmat_cc <- dmat_c[,1:50]

# remove rows of all 0s
dmat_cc <- dmat_cc[!apply(dmat_cc, 1, sum)==0,]

# reorder
odr <- rev(order(apply((dmat_cc==1 | dmat_cc==3), 1, sum)))
dmat_cc <- dmat_cc[odr,]

# do this with the lower level raster functions

dmat_cc[dmat_cc==0] <- NA
dmat_cc[dmat_cc==1] <- '#000000'
dmat_cc[dmat_cc==2] <- '#FF0000'
#dmat_cc[dmat_cc==3] <- '#0000FF'
dmat_cc[dmat_cc==3] <- '#000000'

dmat_cc <- as.raster(dmat_cc)

## ---------

#pdf(file='manuscript/Figs/protein_map_v4.pdf', width=3.5, height=6)
png(file='manuscript/Figs/protein_map_v9.png', width=3.5, height=6, units='in', res=250)

par(mar=c(4.5, 2, 1.25, 2.5),
    cex.axis=0.75, cex.lab=0.75,
    las=1, pty='m')

plot(0, 0, type='n', 
     xlim=c(0, ncol(dmat_cc)), ylim=c(0, nrow(dmat_cc)),
     xaxt='n', yaxt='n', xaxs='i', yaxs='i',
     xlab=NA, ylab=NA)

# draw rectangle around new peptides

# bottom edge of the rect
ry0 <- sum(apply(dmat_c[,1:50], 1, function(x) { any(x == 1)}))
ry1 <- nrow(dmat_cc)
rect(xleft=0, xright=ncol(dmat_cc), ybottom=ry0, ytop=ry1, 
     col=rgb(1, 0, 0, 0.15), border=NA)

rasterImage(dmat_cc[nrow(dmat_cc):1,], 0, 0, ncol(dmat_cc), nrow(dmat_cc), interpolate=F)


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
legend(x=0, y=-650, xjust=0, yjust=1,
       c('Not Quantified', 'Spectra', 'DART-ID (Upgraded)'),
       #col=c('#FFFFFF', '#000000', '#FF0000'),
       pch=22, pt.cex=1.5, pt.lwd=1, 
       pt.bg=c('#FFFFFF', '#000000', '#FF0000'), 
       ncol=2, bty='n', cex=0.85,
       x.intersp=1.1, y.intersp=1, text.width=18, xpd=T)

mtext('Experiment', 1, line=1, cex=1)
mtext('Distinct Peptide ID * 1000', 2, line=1, cex=1, las=3)
#mtext('Quantified Peptide Coverage Increase, FDR < 1%', 3, line=0.1, cex=0.75, font=2)
mtext('Peptide Coverage Increase, FDR < 1%', 3, line=0.1, cex=1, font=2)


text(56.5, 5.85e3, expression('DART-ID'[2]), xpd=T, srt=270, adj=c(0.5, 0.5), cex=1)
text(52.5, 5.85e3, 'Including peptides without confident spectra',
     xpd=T, srt=270, adj=c(0.5, 0.5), cex=0.75)
text(56.5, 1650, expression('DART-ID'[1]), xpd=T, srt=270, adj=c(0.5, 0.5), cex=1)
text(52.5, 1650, 'Peptides with confident spectra',
     xpd=T, srt=270, adj=c(0.5, 0.5), cex=0.75)


dev.off()


# two-panel ---------------------------------------------------------------

ry0 <- sum(apply(dmat_c[,1:50], 1, function(x) { any(x == 1)}))
ry1 <- nrow(dmat_cc)

png(file='manuscript/Figs/protein_map_v10.png', width=3.5, height=6, units='in', res=250)

layout(c(1, 2), heights=c((ry1-ry0)/ry0), 1)

par(oma=c(0,3,0.5,0.75))

par(mar=c(4.5, 0, 0, 0))

plot(0, 0, type='n', 
     xlim=c(0, ncol(dmat_cc)), ylim=c(0, ry1-ry0),
     xaxt='n', yaxt='n', xaxs='i', yaxs='i',
     xlab=NA, ylab=NA)

rasterImage(dmat_cc[nrow(dmat_cc):ry0,], 0, 0, ncol(dmat_cc), ry1-ry0, interpolate=F)

#axis(1, tck=-0.02, at=seq(0, 50, by=10), label=NA, mgp=c(0, 0.1, 0))
axis(2, tck=-0.02, at=seq(1000-(ry0-3000), 11000, by=1000), 
     labels=seq(4, 14, 1), mgp=c(0, 0.5, 0), las=1)

mtext(expression('Distinct Peptides'%*%1000), side=2, cex=1, line=1.5, at=-500)

mtext(expression(''%up%' DART-ID'[2]), side=1, line=0.25, cex=1)
mtext(expression('Peptides without confident spectra'), side=1, line=1.15, cex=1)

par(mar=c(4.5, 0, 0, 0))

plot(0, 0, type='n', 
     xlim=c(0, ncol(dmat_cc)), ylim=c(0, ry0),
     xaxt='n', yaxt='n', xaxs='i', yaxs='i',
     xlab=NA, ylab=NA)

rasterImage(dmat_cc[ry0:1,], 0, 0, ncol(dmat_cc), ry0, interpolate=F)

axis(1, tck=-0.02, at=seq(0, 50, by=10), mgp=c(0, 0.1, 0))
axis(2, tck=-0.02, at=seq(0, 11000, by=1000), labels=seq(0, 11, 1), mgp=c(0, 0.5, 0))

mtext('Experiment', 1, line=1, cex=1)
mtext(expression(''%down%' DART-ID'[1]), side=3, line=0.9, cex=1)
mtext(expression('Peptides with confident spectra'), side=3, line=-0.1, cex=1)

legend(x=0, y=-850, xjust=0, yjust=1,
       c('Not Quantified', 'Spectra', 'DART-ID (Upgraded)'),
       #col=c('#FFFFFF', '#000000', '#FF0000'),
       pch=22, pt.cex=1.2, pt.lwd=1, 
       pt.bg=c('#FFFFFF', '#000000', '#FF0000'), 
       ncol=2, bty='n', cex=0.85,
       x.intersp=0.9, y.intersp=1, text.width=17, xpd=T)

dev.off()

