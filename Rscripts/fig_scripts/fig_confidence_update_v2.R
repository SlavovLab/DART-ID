library(tidyverse)
#library(gridExtra)
#library(grid)
#library(gtable)
source('Rscripts/lib.R')

## Peptide Update ------

pdf(file='manuscript/Figs/confidence_update_v5.pdf', width=3.5, height=3)

layout(rbind(c(1, 2),
             c(1, 3)),
       widths=c(3, 1))

par(oma=c(0,0,0,0.5),
    mar=c(2.5,3,0.5,0.75),
    xaxs='i', yaxs='i', pty='m', cex.axis=1)

plot(0, 0, type='n',
     xlim=c(21, 25), ylim=c(-0.1, 1.85),
     xaxt='n', yaxt='n',
     #xlab="Retention Time (min)", ylab="Density",
     #main='Confidence Update',
     xlab=NA, ylab=NA)

mu <- 25.01534
exp <- 0.95

denx <- seq(21, 25, length.out=1000)
dist_sd <- 0.275
polygon(denx, dnorm(denx, mean=mu*exp, sd=dist_sd), 
        col=rgb(0, 0, 1, 0.3))
lines(denx[400:1000], dnorm(denx[400:1000], mean=mu*exp, sd=dist_sd), col='blue', lwd=2)

segments(x0=mu*exp, x1=mu*exp,
         y0=0, y1=dnorm(mu*exp, mean=mu*exp, sd=dist_sd),
         col='blue', lwd=2, lty=2)
abline(h=0, col='black', lwd=2)

# null distribution
polygon(c(21,denx,25), c(0,(dlnorm(denx, meanlog=3.1, sdlog=0.07) * 1),0), 
        col=rgb(1, 0, 0, 0.3), border=NA)
lines(denx, dlnorm(denx, meanlog=3.1, sdlog=0.07) * 1, col='red', lwd=2)

# points
x1 <- mu*exp-1.6
x2 <- mu*exp-0.3
points(c(x1, x2), rep(0,2), pch=c(25, 24), bg='black', col='black', cex=2)

#text(x=mu*exp,y=1.65, labels="Peptide X\nExperiment A", font=1, cex=1.2)

axis(side=1, tck=-0.02, mgp=c(0, 0.3, 0))
axis(side=2, tck=-0.02, las=1, mgp=c(0, 0.4, 0))

mtext('Retention time (min)', side=1, line=1.5, cex=1)
mtext('Probability density', side=2, line=1.85, cex=1, las=3)
#mtext('Confidence Update', side=3, line=0.25, cex=1, font=2)
#mtext('Confidence Update, Peptide X in Exp. A', side=3, outer=T, line=0, cex=1, font=2)

# legend(x=21.1, y=0.35, xjust=0, yjust=0,
#        c('Inferred RT', 'Aligned\nCanonical RT', 'Null RT', 'PSM 1', 'PSM 2'), 
#        pch=c(NA, NA, NA, 16, 17),
#        col=c('blue', 'blue', 'red', 'black', 'black'), 
#        lty=c(1, 2, 1, NA, NA), lwd=c(2, 2, 2, NA, NA),
#        bty='n', cex=1.1, pt.cex=2, x.intersp=1, y.intersp=1.25,
#        adj=c(0, 0.5), inset=c(0.1,0))
legend(#x=21.1, y=0.35,
       'topleft',
       xjust=0, yjust=0,
       c('Aligned Canonical RT', 'Inferred RT', 'Null RT', 'PSM 1', 'PSM 2'), 
       pch=c(NA, 22, 22, 25, 24),
       col=c('blue', 'blue', 'red', 'black', 'black'), 
       pt.bg=c(NA, rgb(0,0,1,0.3), rgb(1,0,0,0.3), 'black', 'black'),
       pt.lwd=c(NA, 2, 2, NA, NA),
       pt.cex=c(NA, 3, 3, 1.5, 1.5),
       lty=c(2, 1, 1, NA, NA), lwd=c(2, NA, NA, NA, NA),
       bty='n', cex=1, x.intersp=1, 
       #y.intersp=1.25,
       #y.intersp=c(1, 1, 0.8, 0.8, 0.8),
       y.intersp=1.25,
       adj=c(0, 0.5), inset=c(0.03,-0.03))

par(#oma=c(0, 0, 1, 0),
    mar=c(0.75, 1, 2.5, 0.5),
    pty='m')

plot(0, 0, type='n',
     xlim=c(0, 2.75), ylim=c(0, 3.5),
     xlab=NA, ylab=NA,
     xaxs='i', yaxs='i',
     xaxt='n', yaxt='n')

barplot(-log10(c(5e-2, 0.3)), width=1, space=0.25, add=T,
        col=c(av[1], av[2]),
        xaxt='n', yaxt='n')

mtext('   PSM 1', side=3, line=0.1, cex=0.8, font=2)
points(0.25, 3.85, pch=25, cex=1.25, xpd=T,bg='black')
mtext('ID Confidence', 2, line=0.25, cex=1, at=-0.5, las=3)

par(#oma=c(1,0,0,0),
    mar=c(2.5, 1, 1, 0.5),
    pty='m')

plot(0, 0, type='n',
     xlim=c(0, 2.75), ylim=c(0, 3.5),
     xlab=NA, ylab=NA,
     xaxs='i', yaxs='i',
     xaxt='n', yaxt='n')

barplot(-log10(c(5e-2, 1e-3)), width=1, space=0.25, add=T,
        col=c(av[1], av[2]),
        xaxt='n', yaxt='n')

par(xpd=T)
legend(x=-0.5, y=0, c('Spectra', 'DART-ID'), 
       pch=22, pt.cex=2, pt.bg=c(av[1], av[2]), cex=0.9, col='black',
       y.intersp=1, bty='n')

mtext('   PSM 2', side=3, line=0.1, cex=0.8, font=2)
points(0.25, 3.8, pch=24, cex=1.25, xpd=T, bg='black')

dev.off()
