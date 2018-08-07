library(tidyverse)
library(RColorBrewer)
source('Rscripts/lib.R')

lighten <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col*factor
  col[col > 255] <- 255
  col <- rgb(t(col), maxColorValue=255)
  col
}
unsaturate <- function(color, factor=1.4) {
  col <- col2rgb(color)
  colhsv <- rgb2hsv(col)
  colhsv[2,] <- colhsv[2,]/factor
  hsv(colhsv[1,], colhsv[2,], colhsv[3,])
}

## ----------

pdf(file='manuscript/Figs/alignment_demo_v7.pdf', width=3.5, height=2.5)
par(oma=c(0,0,0,0),
    mar=c(2,2.5,1.25,5.5),
    cex.lab=0.85, cex.axis=0.85, cex.main=1,
    pty='m')

set.seed(15)

exps <- c(1.25, 1.03, 0.65)
exp_intercepts <- c(-3.4, 0, 7.1)

#mus <- c(23.13, 24.19, 25.51, 26.94, 28.68, 29.79)
#mus <- c(23.53, 25.51, 26.94, 28.68)
mus <- c(20.5, 21.16, 21.7, 22.4)
points_per_exp <- 1
musx <- rep(rep(mus, each=points_per_exp),length(exps))
rts <- musx * rep(exps, each=points_per_exp*length(mus)) + 
  rep(exp_intercepts, each=points_per_exp*length(mus))
error <- rnorm(length(rts), 0, 0.4)
rts <- rts + error

pep_cols <- brewer.pal(6, 'Set1')
# swap 1 and 2
pep_cols[1] <- brewer.pal(6, 'Set1')[2]
pep_cols[2] <- brewer.pal(6, 'Set1')[1]

# randomly set some rts to 0 to simulate missing data
#num_to_remove <- 0.2
#rts[sample.int(length(rts), size=ceiling(num_to_remove*length(rts)))] <- NA

# remove peptide beta from exp A and B
rts[c(2, 6)] <- NA
# remove beta from B
#rts[6] <- NA
# remove gamma from A
rts[3] <- NA
# remove alpha from C
rts[9] <- NA

rts[9] <- rts[9] + 0.4
rts[8] <- rts[8] - 0.7


plot(0, 0, type='n', xlim=c(19.5,22.75), ylim=c(19.75,25.5),
     xaxt='n',yaxt='n', xlab=NA, ylab=NA)

# lines for mus
segments(x0=mus, x1=mus, y0=0, 
         #y1=mus*rep(exps[1], length(mus)), 
         y1=33.5,
         #col='grey90', 
         #col=paste0(pep_cols[1:4], '88'),
         col=lighten(unsaturate(pep_cols[1:4], 2.5), 1.2),
         lty=1, lwd=3)

points(musx, rts, 
       #col=rep(cols, each=length(mus)),
       col=rep(pep_cols[1:4], length(exps)),
       pch=4, 
       lwd=3, cex=1.5)

# plot missing data
points(x=musx[is.na(rts)], 
       y=(rep(exps, each=length(mus))[is.na(rts)] * musx[is.na(rts)]) + rep(exp_intercepts, each=length(mus))[is.na(rts)], 
       pch=1, 
       #col=paste0(rep(pep_cols[1:4], length(exps))[is.na(rts)], '66'),
       col=lighten(unsaturate(rep(pep_cols[1:4], length(exps))[is.na(rts)], 2.5), 1.2),
       lwd=3, cex=1.5)

# lines for experiments
#cols <- brewer.pal(3, 'Dark2')
cols <- rep(rgb(0,0,0), 3)
for(i in 1:length(exps)) {
  abline(a=exp_intercepts[i], b=exps[i], 
         col=paste0(cols[i],'FF'), 
         lty=1, 
         #lty=c(4, 2, 1)[i],
         lwd=2)
}

#axis(side=1, tck=-0.02, at=seq(22, 32, by=2), 
#     labels=seq(12, 22, by=2), mgp=c(0, 0.1, 0))
#axis(side=2, tck=-0.02, at=seq(22, 34, by=2), 
#     labels=seq(12, 24, by=2), mgp=c(0, 0.4, 0), las=1)
axis(side=1, at=seq(0, 30, by=1),
     tck=-0.02, mgp=c(0, 0.1, 0))
axis(side=2, at=seq(0, 30, by=1),
     tck=-0.02, mgp=c(0, 0.4, 0), las=1)

mtext('Reference RT (min)', side=1, line=1)
mtext('Observed RT (min)', side=2, line=1.4, las=3)

mtext('Experiment A', side=4, at=25.30, las=1, line=0.25, col=cols[1])
mtext('Experiment B', side=4, at=23.60, las=1, line=0.25, col=cols[2])
mtext('Experiment C', side=4, at=22.00, las=1, line=0.25, col=cols[3])

text(mus, par('usr')[4]+0.4, labels=parse(text=c('beta', 'alpha', 'gamma', 'theta')),
     col=pep_cols[1:4], xpd=T, cex=1.25, font=2)
text(20.2, par('usr')[4]+0.4, 'Peptide:', adj=c(1, 0.5), xpd=T, cex=1)

#mtext('       Aligning Exps with Canonical RTs', side=3, outer=T, cex=1, line=0, font=2, adj=0)

dev.off()