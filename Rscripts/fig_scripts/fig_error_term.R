library(tidyverse)
source('Rscripts/lib.R')

## ------------

pdf(file='manuscript/Figs/error_term_v2.pdf', width=4, height=3)

layout(rbind(c(1,1,3,3),
             c(2,2,3,3)))

par(mgp=c(0.2, 0.5, 0),
    #pin=c(2,2), asp=1,
    pty='s', las=1,
    cex.lab=1, cex.axis=1.25, cex.main=1,
    xaxs='r', yaxs='r')

set.seed(8)

# exp A
par(mar=c(2, 1, 2, 1))

mus <- runif(15, min=20, max=30)
exp <- 0.95
rts <- mus * exp
rts <- rts + rnorm(length(rts), 0, 0.75)

plot(mus, rts,
     xlim=c(20, 28), ylim=c(20, 28),
     xlab=NA, ylab=NA,
     xaxt='n', yaxt='n',
     pch='x', cex=1.5)

abline(a=0, b=exp, lty=1, lwd=1.5)

segments(x0=mus, x1=mus, y0=rts, y1=mus*exp, 
         col='red', lwd=2, lty=1)

text(20, 27.25, labels='Experiment A', cex=1, adj=c(0, 0))
mtext('Observed RT (min)', 2, line=2.1, las=3, at=18, cex=1)
mtext('...', 1, cex=2, line=0.5, adj=0.5)

#axis(side=1, at=seq(20, 28, 2), labels=c(20, NA, 24, NA, 28), 
#     tck=-0.02, padj=-0.1)
axis(side=2, at=seq(20, 28, 2), labels=c(20, NA, 24, NA, 28), 
     tck=-0.02, padj=0.2)

# exp B
par(mar=c(4, 1, 0, 1))

#mus <- runif(15, min=20, max=28)
exp <- 0.95
rts <- mus * exp
rts <- rts + rnorm(length(rts), 0, 0.75)

plot(mus, rts,
     xlim=c(20, 28), ylim=c(20, 28),
     xlab=NA, ylab=NA,
     xaxt='n', yaxt='n',
     pch='x', cex=1.5)

abline(a=0, b=exp, lty=1, lwd=1.5)

segments(x0=mus, x1=mus, y0=rts, y1=mus*exp, 
         col='red', lwd=2, lty=1)

axis(side=1, at=seq(20, 28, 2), labels=c(20, NA, 24, NA, 28), 
     tck=-0.02, padj=-0.1)
axis(side=2, at=seq(20, 28, 2), labels=c(20, NA, 24, NA, 28), 
     tck=-0.02, padj=0.2)

text(20, 27.25, labels='Experiment Z', cex=1, adj=c(0, 0))
mtext('Canonical RT (min)', 1, line=2, adj=0.5, cex=1)

# error bars
par(pty='m', 
    mar=c(3, 0.5, 2.25, 1),
    cex.main=1.5,
    mgp=c(0, 0.5, 0))

plot(0, 0, type='n',
     xlim=c(-0.25, 10.5), ylim=c(0, 1),
     xlab=NA, ylab=NA,
     xaxs='i', yaxs='i',
     xaxt='n', yaxt='n')

bars <- c(0.05, 0.04, 0.06, 0, 0, 0, 0.05, 0.9)
barplot(bars, width=1, space=0.25, add=T,
        col=c(rep(rgb(1, 0, 0, 0.3), 7), rgb(1, 0, 0, 1)),
        xaxt='n', yaxt='n')

text(x=5.75, y=0.03, labels='...', cex=3, adj=c(0.5, 0))

axis(side=1, tck=-0.02, padj=-0.1, 
     #at=seq(0.75, 9.5, by=1.25), 
     at=c(0.75, 2, 3.25, 8.25, 9.5),
     labels=c('A', 'B', 'C', 'Z', 'All'),
     cex.axis=1)
axis(side=2, tck=-0.02, padj=0.2,
     at=seq(0, 1, 0.2),
     labels=seq(0, 5, length.out=6))

mtext('Single Error Term', 3, line=0.5, cex=1, font=2)
mtext('Experiment', 1, line=1.75, cex=1)
mtext('Total Residual RT (min * 1000)', 2, line=1.75, las=3, cex=1)

dev.off()
