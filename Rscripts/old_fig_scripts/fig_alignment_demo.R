library(tidyverse)
library(RColorBrewer)
source('Rscripts/lib.R')

## ----------

pdf(file='manuscript/Figs/alignment_demo_v4.pdf', width=3.5, height=2.5)
par(oma=c(0,0,0,1.5),
    mar=c(2,1.5,0.5,0.1),
    cex.lab=0.85, cex.axis=0.85, cex.main=1,
    pty='s')

set.seed(15)

exps <- c(1.06, 1.01, 0.95)
exp_intercepts <- c(-0.1, 0, 0.1)

#mus <- runif(20, min=20,max=30)
# mus <- c(20.30265, 20.63603, 20.92263, 22.13038, 22.43004, 23.09638, 23.49773,
#          23.96849, 25.01534, 25.65994, 25.95323, 26.40581, 26.94700, 27.44276,
#          27.61737, 28.42095, 28.68403, 29.08916, 29.23397, 29.79338)
mus <- c(23.13, 24.19, 25.51, 26.94, 28.68, 29.79)
points_per_exp <- 1
musx <- rep(rep(mus, each=points_per_exp),length(exps))
rts <- musx * rep(exps, each=points_per_exp*length(mus))
error <- rnorm(length(rts), 0, 0.5)
rts <- rts + error
#rts[49] <- rts[49] - 0.6

# randomly set some rts to 0 to simulate missing data
num_to_remove <- 0.2
rts[sample.int(length(rts), size=ceiling(num_to_remove*length(rts)))] <- NA

#cols <- rep(rgb(0,0,0,0.3), length(musx))
#cols[-(1:(2*length(mus)))] <- rgb(0, 0, 1, 0.5)
#cols[49] <- 'blue'

plot(0, 0, type='n',
     xlim=c(22,30), ylim=c(22,32),
     xaxt='n',yaxt='n',
     xlab=NA, ylab=NA)
#points(musx[-49], rts[-49], pch='+', cex=1.25, col=rgb(0,0,0,0.4))
#points(musx[49], rts[49], pch='x', cex=1.25, col=rgb(0,0,1))

# lines for mus
#abline(v=mus[-9], col=rgb(0,0,0,0.1), lty=1)
segments(x0=mus, x1=mus,
         y0=0, 
         #y1=mus*rep(exps[1], length(mus)), 
         y1=32.25,
         #col=rgb(0,0,0,0.2), 
         #col=c('green', 'blue', 'purple'),
         col='grey90',
         lty=1, lwd=2)

# lines for experiments
#abline(a=0, b=1, col='black')
#abline(a=0, b=exps[3], col='blue')
ltys <- rev(c(1, 4, 5))
#cols <- c('green', 'blue', 'purple')
cols <- brewer.pal(3, 'Dark2')
for(i in 1:length(exps)) {
  abline(a=0, b=exps[i], 
         #col=rgb(0,0,1,1), 
         col=paste0(cols[i],'FF'),
         #lty=ltys[i], 
         lty=1,
         lwd=2)
}

points(musx, rts, 
       #col=rgb(0,0,0,1),
       col=rep(cols, each=length(mus)),
       #pch=rep(1:6, length(exps)), 
       pch=3, lwd=3, cex=1.5)

# 
# #abline(v=mus[9], col=rgb(0,0,0,0.4), lty=2, lwd=1.5)
# segments(x0=mus[9], x1=mus[9], y0=0, y1=rts[9], 
#          col=rgb(0,0,0,0.4), lty=1, lwd=1.5)




# # draw the sideways distribution
# polygon(mus[9] - (dnorm(seq(22, 26, length.out=200), mean=(mus[9] * exps[3]), sd=0.4)) * 0.75,
#         seq(22, 26, length.out=200), 
#         col=rgb(0, 0, 1, 0.3), border=NA)
# # draw mean of sideways distribution
# segments(x0=mus[9], x1=mus[9] - (dnorm(mus[9] * exps[3], mean=(mus[9] * exps[3]), sd=0.4) * 0.75), 
#          y0=mus[9] * exps[3], y1=mus[9] * exps[3], 
#          col='blue', lty=2)
# text(x=mus[9]+0.3, y=rts[49]+0.5, labels="Inferred RT Distribution\nPeptide X, Experiment A", 
#      cex=0.7, adj=c(0,1))


axis(side=1, tck=-0.02, at=seq(22, 32, by=2), mgp=c(0, 0.1, 0))
axis(side=2, tck=-0.02, at=seq(22, 34, by=2), mgp=c(0, 0.4, 0), las=1)

mtext('Canonical RT (min)', side=1, line=1)
mtext('Observed RT (min)', side=2, line=1.4, las=3)
#mtext('Global RT Alignment', side=3, cex=1, line=0.25, font=2)

mtext('Exp A', side=4, at=32.15, las=1, line=0.25, col=cols[1])
mtext('Exp B', side=4, at=30.75, las=1, line=0.25, col=cols[2])
mtext('Exp C', side=4, at=28.85, las=1, line=0.25, col=cols[3])

# legend('topleft', c('Aligned Canonical RT', 'PSM'), 
#        pch=c(NA,'x'), col=c(rgb(0,0,0,0.5), rgb(0,0,0,0.5)), 
#        lty=c(1, NA), lwd=2, pt.cex=1, cex=0.85, 
#        bty='n', inset=c(0.03,0), y.intersp=1.1)
#legend('topleft', c(''))

dev.off()