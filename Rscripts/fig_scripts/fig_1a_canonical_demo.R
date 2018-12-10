library(tidyverse)
library(RColorBrewer)
source('Rscripts/lib.R')

# ----------------------------------------------------------------

pdf(file='manuscript/Figs/canonical_demo_v4.pdf', width=3.5, height=2.5)

par(oma=c(0,0,0.25,0),
    mar=c(2,3,0.25,6), cex.axis=0.85)

set.seed(9)

mus <- c(20.5, 21.16)
num_exps <- 6

# each experiment has a random shift
#exps <- runif(num_exps, -0.3, 0.3)
exps <- c(0.13, 0.04, -0.14, -0.2, 0, 0.2)

# and added onto those shifts is noise
points_per_exp <- 1
rts <- rep(mus, each=num_exps*points_per_exp) + rep(rep(exps, each=points_per_exp), length(mus))
# add noise
rts <- rts + rnorm(length(rts), mean=0, sd=0.05)
# randomly remove some data
#remove_fraction <- 0.7
#rts[sample.int(length(rts), ceiling(remove_fraction*length(rts)))] <- NA

cols <- brewer.pal(6, 'Set1')

plot(0, 0, type='n',
     ylim=c(20, 21.3), xlim=c(0.5, num_exps+0.5),
     xaxt='n', yaxt='n', xlab=NA, ylab=NA)
abline(v=1:num_exps, col=rgb(0,0,0,0.1), lty=1, lwd=2)
abline(h=mus, col=cols[c(2, 1)], lty=1, lwd=3)

x <- rep(rep(seq(1,num_exps), 2), each=points_per_exp)
# remove some rts
rts[c(3, 7, 8, 12)] <- NA

points(x, rts, 
       pch=4,
       #pch=rep(c(1,4), each=num_exps*points_per_exp), 
       col=rep(cols[c(2, 1)], each=num_exps*points_per_exp),
       cex=1.25, lwd=2)

points(c(1, 2, 6), rep(mus[2], 3), pch=1, col=paste0(cols[1],'88'), cex=1.25, lwd=2)
points(c(3), rep(mus[1], 1), pch=1, col=paste0(cols[2], '88'), cex=1.25, lwd=2)

# l <- legend(21.4, 4, c('PSMs:', 'Canonical RT:'), cex=0.9,
#             y.intersp=1.5,
#             xpd=T, xjust=0, yjust=0.5, bty='n')
# 
# points(rep(l$text$x[1], 2) + l$rect$w - c(0.70, 0.50), 
#        rep(l$text$y[1], 2), 
#        pch=c(1, 4), lwd=2, cex=1.25, xpd=T)
# segments(x0=l$text$x[2]+l$rect$w-0.25, x1=l$text$x[2]+l$rect$w-0.25, 
#          y0=l$text$y[2]-0.35, y1=l$text$y[2]+0.35,
#          col='red', lty=1, lwd=2, xpd=T)

#par(lheight=0.85)
text(7.3, 20.125, 'Peptide:', cex=1, adj=c(0, 0.5), xpd=T)
legend(6.75, 19.9, c('Observed', 'Missing'), xjust=0, yjust=0.5, 
       pch=c(4, 1), pt.cex=1.25, lwd=2, lty=NA, text.width=0.3,
       col=c(rgb(0,0,0,1), rgb(0,0,0,0.3)),
       cex=0.85, x.intersp=0.4, y.intersp=1.25, bty='n', xpd=T)
par(lheight=1)

axis(2, at=seq(19, 23, by=0.5), tck=-0.02, mgp=c(0, 0.4, 0), las=1)
axis(1, at=1:num_exps, tck=-0.02, mgp=c(0, 0.1, 0),
     labels=LETTERS[1:num_exps])

mtext('Reference RT', side=4, at=mus[1]+0.08, 
      line=0.3, las=1, col=cols[2])
mtext(bquote(.('Peptide ')*beta), side=4, at=mus[1]-0.08, 
      line=0.3, las=1, col=cols[2])

mtext('Reference RT', side=4, at=mus[2]+0.08, 
      line=0.3, las=1, col=cols[1])
mtext(bquote(.('Peptide ')*alpha), side=4, at=mus[2]-0.08, 
      line=0.3, las=1, col=cols[1])

mtext('Experiment', side=1, line=1, cex=1)
mtext('Retention time (min)', side=2, line=2.1, cex=1)
#mtext('      Selecting Canonical RT from PSM RTs', side=3, outer=T, line=0, cex=1, font=2, adj=0)

dev.off()