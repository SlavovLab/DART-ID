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

pdf(file='manuscript/Figs/alignment_demo_v5.pdf', width=3.5, height=2.5)
par(oma=c(0,0,1,0),
    mar=c(2,0,1.25,0),
    cex.lab=0.85, cex.axis=0.85, cex.main=1,
    pty='s')

set.seed(15)

exps <- c(1.15, 1.03, 0.85)
exp_intercepts <- c(-1, 0, 2.5)

#mus <- c(23.13, 24.19, 25.51, 26.94, 28.68, 29.79)
mus <- c(23.13, 25.51, 26.94, 28.68)
points_per_exp <- 1
musx <- rep(rep(mus, each=points_per_exp),length(exps))
rts <- musx * rep(exps, each=points_per_exp*length(mus)) + 
  rep(exp_intercepts, each=points_per_exp*length(mus))
error <- rnorm(length(rts), 0, 0.4)
rts <- rts + error

pep_cols <- brewer.pal(6, 'Set1')

# randomly set some rts to 0 to simulate missing data
#num_to_remove <- 0.2
#rts[sample.int(length(rts), size=ceiling(num_to_remove*length(rts)))] <- NA

# remove peptide alpha from exp A and B
rts[c(1, 5)] <- NA
# remove beta from B
#rts[6] <- NA
# remove gamma from A
rts[3] <- NA
# remove beta from C
rts[10] <- NA

rts[9] <- rts[9] + 0.4
rts[8] <- rts[8] - 0.7


plot(0, 0, type='n', xlim=c(22,30), ylim=c(22,33),
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

axis(side=1, tck=-0.02, at=seq(22, 32, by=2), 
     labels=seq(12, 22, by=2), mgp=c(0, 0.1, 0))
axis(side=2, tck=-0.02, at=seq(22, 34, by=2), 
     labels=seq(12, 24, by=2), mgp=c(0, 0.4, 0), las=1)

mtext('Canonical RT (min)', side=1, line=1)
mtext('Observed RT (min)', side=2, line=1.4, las=3)

mtext('Exp A', side=4, at=33.55, las=1, line=0.25, col=cols[1])
mtext('Exp B', side=4, at=31.30, las=1, line=0.25, col=cols[2])
mtext('Exp C', side=4, at=28.30, las=1, line=0.25, col=cols[3])

text(mus, par('usr')[4]+0.75, labels=parse(text=c('alpha', 'beta', 'gamma', 'delta')),
     col=pep_cols[1:4], xpd=T, cex=1.25, font=2)
text(21.5, par('usr')[4]+0.75, 'Peptide:', adj=c(1, 0.5), xpd=T, cex=1)

mtext('       Aligning Exps with Canonical RTs', side=3, outer=T, cex=1, line=0, font=2, adj=0)

dev.off()