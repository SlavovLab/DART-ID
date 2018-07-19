library(tidyverse)
source('Rscripts/lib.R')

ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180621_2/ev_updated.txt')

## load cormat data -----

source('Rscripts/validation_cormats.R')

# exclude carriers from cor mats
cor_type_a <- cor_type_a[-c(1, 5), -c(1, 5)]
cor_type_b <- cor_type_b[-c(1, 5), -c(1, 5)]

## ------

pdf(file='manuscript/Figs/cormats_type_v2.pdf', width=4.5, height=3.5)

# 1, 2 = cell type cor mats, 3 = colorbar
#layout(rbind(c(1, 1, 1, 3),
#             c(2, 2, 2, 3)))
layout(t(c(1, 2, 3)), widths=c(1, 1, 0.25))

par(oma=c(5, 2.5, 5, 0))

# cell type cormats

colfunc <- colorRampPalette(c('blue', 'white', 'red'))
ncols <- 50

cell_labels = parse(text=c('U937[3]', 'U937[2]', 'U937[1]',
                           'Jurkat[3]', 'Jurkat[2]', 'Jurkat[1]'))

par(mar=c(0,0,0,0.5), pty='s',
    cex.axis=0.75, cex.lab=0.75)

image(cor_type_a[nrow(cor_type_a):1,], zlim=c(-1, 1), col=colfunc(ncols),
      xlab=NA, ylab=NA, xaxt='n', yaxt='n')

# axis(1, at=seq(0, 1, length.out=nrow(cor_type_a)), labels=NA,
#      tck=-0.02)
axis(1, at=seq(0, 1, length.out=nrow(cor_type_a)),
     #labels=c('J', 'J', 'J', '50J', 'U', 'U', 'U', '50U'), 
     labels=rev(cell_labels),
     tck=-0.01, las=3, mgp=c(0, 0.25, 0))
axis(2, at=seq(0, 1, length.out=nrow(cor_type_a)),
     #labels=c('50U', 'U', 'U', 'U', '50J', 'J', 'J', 'J'), 
     labels=cell_labels,
     tck=-0.01, las=1, mgp=c(0, 0.25, 0))

#mtext(paste0('Spectra - ', length(old_prots), ' proteins'), 2, line=3, cex=0.75, las=3, font=1)
#mtext(paste0('Spectra'), side=2, line=2.5, cex=0.75, las=3, font=1)
mtext(paste0('Spectra'), side=1, line=2.5, cex=0.75, font=1)
#mtext('J vs. U Correlation', side=3, line=0.5, cex=0.8, font=2, las=1)

par(mar=c(0,0.5,0,0))
image(cor_type_b[nrow(cor_type_b):1,], zlim=c(-1, 1), col=colfunc(ncols),
      xlab=NA, ylab=NA, xaxt='n', yaxt='n')

# axis(1, at=seq(0, 1, length.out=nrow(cor_type_a)), labels=NA,
#      tck=-0.02)
axis(1, at=seq(0, 1, length.out=nrow(cor_type_a)),
     #labels=c('J', 'J', 'J', '50J', 'U', 'U', 'U', '50U'),
     labels=rev(cell_labels),
     tck=-0.01, las=3, mgp=c(0, 0.25, 0))
axis(2, at=seq(0, 1, length.out=nrow(cor_type_a)), labels=NA,
     tck=-0.02)
# axis(2, at=seq(0, 1, length.out=nrow(cor_type_a)),
#      #labels=c('50U', 'U', 'U', 'U', '50J', 'J', 'J', 'J'), 
#      labels=cell_labels,
#      tck=-0.01, las=1, mgp=c(0, 0.25, 0))

#mtext(paste0('DART-ID - ', length(new_prots), ' proteins'), 2, line=3, cex=0.75, las=3, font=1)
#mtext(paste0('DART-ID'), side=2, line=2.5, cex=0.75, las=3, font=1)
mtext(paste0('DART-ID'), side=1, line=2.5, cex=0.75, font=1)

par(mar=c(1.5,0.75,1.5,1.75), pty='m')
image(matrix(seq(-1, 1, length.out=ncols), ncol=ncols), col=colfunc(ncols),
      xlab=NA, ylab=NA, xaxt='n', yaxt='n')
axis(4, at=seq(0, 1, length.out=11), labels=seq(-1, 1, length.out=11), 
     tck=-0.1, las=1, mgp=c(0, 0.4, 0))

mtext('Jurkat vs. U937 Correlation, Separate Proteins    ', side=3, line=-1, cex=1, font=2, las=1, outer=T)

dev.off()


# ratios ------------------------------------------------------------------

pdf(file='manuscript/Figs/cormats_ratio.pdf', width=2, height=3)

# 1, 2 = j/u ratio cormats, 3 = colorbar
layout(rbind(c(1, 1, 1, 3),
             c(2, 2, 2, 3)))

par(oma=c(2, 0, 2, 0))

colfunc <- colorRampPalette(c('white', 'red'))
ratio_expressions <- c(
  'J[1]/U[1]', 'J[1]/U[2]', 'J[1]/U[3]',
  'J[2]/U[1]', 'J[2]/U[2]', 'J[2]/U[3]',
  'J[3]/U[1]', 'J[3]/U[2]', 'J[3]/U[3]')
ratio_expressions <- parse(text=ratio_expressions)

par(mar=c(0,3.5,0,0), pty='s',
    cex.axis=0.75, cex.lab=0.75)

image(cor_ratio_a[nrow(cor_ratio_a):1,], zlim=c(0.5, 1), col=colfunc(ncols),
      xlab=NA, ylab=NA, xaxt='n', yaxt='n')
axis(1, at=seq(0, 1, length.out=nrow(cor_ratio_a)),
     labels=NA, tck=-0.02)
axis(2, at=seq(0, 1, length.out=nrow(cor_ratio_a)),
     tck=-0.01, labels=rev(ratio_expressions), mgp=c(0, 0.4, 0), las=1)

mtext('Log2 J/U Ratio Correlation', 3, line=0.5, cex=0.8, font=2, las=1)

par(mar=c(0,3.5,0,0))

image(cor_ratio_b[nrow(cor_ratio_b):1,], zlim=c(0.5, 1), col=colfunc(ncols),
      xlab=NA, ylab=NA, xaxt='n', yaxt='n')
axis(1, at=seq(0, 1, length.out=nrow(cor_ratio_a)),
     tck=-0.01, labels=ratio_expressions, mgp=c(0, 0.3, 0), las=3)
axis(2, at=seq(0, 1, length.out=nrow(cor_ratio_a)),
     tck=-0.01, labels=rev(ratio_expressions), mgp=c(0, 0.4, 0), las=1)

par(mar=c(2, 0.75, 2, 2), pty='m')
image(matrix(seq(-1, 1, length.out=ncols), ncol=ncols), col=colfunc(ncols),
      xlab=NA, ylab=NA, xaxt='n', yaxt='n')
axis(4, at=seq(0, 1, length.out=6), labels=seq(0.5, 1, length.out=6), 
     tck=-0.1, las=1, mgp=c(0, 0.3, 0))

dev.off()