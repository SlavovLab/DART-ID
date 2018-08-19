library(tidyverse)
library(rmutil)
source('Rscripts/lib.R')

## Peptide Update ------

pdf(file='manuscript/Figs/confidence_update_v10.pdf', width=3.5, height=2.5)

layout(rbind(c(1, 2),
             c(1, 3)),
       widths=c(3.5, 1))

par(oma=c(2,0,0,1),
    mar=c(0.5,3,0.5,0.75),
    xaxs='i', yaxs='i', pty='m', cex.axis=1)

plot(0, 0, type='n',
     xlim=c(20, 25), ylim=c(-0.15, 1.95),
     xaxt='n', yaxt='n',
     #xlab="Retention Time (min)", ylab="Density",
     #main='Confidence Update',
     xlab=NA, ylab=NA)

mu <- 25.25
exp <- 0.95

denx <- seq(20, 25, length.out=1000)
dist_sd <- 0.275
#deny <- dnorm(denx, mean=mu*exp, sd=dist_sd)
deny <- dlaplace(denx, m=mu*exp, s=dist_sd)
polygon(c(denx, 25), c(deny, 0), 
        col=rgb(0, 0, 1, 0.3))
lines(denx[400:1000], deny[400:1000], col='blue', lwd=2)

#segments(x0=mu*exp, x1=mu*exp,
#         y0=0, y1=dnorm(mu*exp, mean=mu*exp, sd=dist_sd),
#         col='blue', lwd=2, lty=2)
abline(h=0, col='black', lwd=2)

# null distribution
den_null <- (dnorm(denx, mean=22, sd=2) * 1)
polygon(c(0, denx,25), c(0, den_null,0), 
        col=rgb(1, 0, 0, 0.3), border=NA)
lines(denx, den_null, col='red', lwd=2)

# points
x1 <- mu*exp-2
x2 <- mu*exp+0.3
points(c(x1, x2), rep(0,2), pch=c(25, 24), bg='white', col='black', lwd=2, cex=2)

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
par(lheight=0.85)
# legend(#x=21.1, y=0.35,
#        'topleft',
#        xjust=0, yjust=0,
#        c('Aligned\nCanonical RT', 
#          'Conditional\nLikelihood\n(PSM+)', 
#          'Conditional\nLikelihood\n(PSM-)'), 
#        pch=c(NA, 22, 22),
#        col=c('blue', 'blue', 'red'), 
#        pt.bg=c(NA, rgb(0,0,1,0.3), rgb(1,0,0,0.3)),
#        pt.lwd=c(NA, 2, 2),
#        pt.cex=c(NA, 4, 4),
#        lty=c(2, 1, 1), lwd=c(2, NA, NA),
#        bty='n', cex=1, x.intersp=1, 
#        #y.intersp=1.25,
#        #y.intersp=c(1, 1, 0.8, 0.8, 0.8),
#        y.intersp=c(0, 0.8, 1.1),
#        adj=c(0, 0.5), inset=c(0.03,-0.03))
#legend(20.25, 1.6, c('Aligned\nReference RT'), xjust=0, yjust=0.5,
#       pch=NA, col='blue', lty=2, lwd=2, 
#       bty='n', cex=1, inset=c(0,0))
text(20.4, 1.7, 'Conditional RT\nLikelihood', adj=c(0, 0.5), cex=1)
legend(20.25, 1.6, c('ID correct', 'ID incorrect'), col=c('blue', 'red'), xjust=0, yjust=1,
       pch=22, pt.bg=c(rgb(0,0,1,0.3), rgb(1,0,0,0.3)), pt.lwd=2, pt.cex=3,
       lty=1, lwd=NA, bty='n', cex=1, y.intersp=1.5, inset=c(0,0), adj=c(0, 0.5))

par(#oma=c(0, 0, 1, 0),
    mar=c(0.25, 1, 1.5, 0.5),
    pty='m')

plot(0, 0, type='n',
     xlim=c(0, 2.75), ylim=c(0, 3.5),
     xlab=NA, ylab=NA,
     xaxs='i', yaxs='i',
     xaxt='n', yaxt='n')

barplot(-log10(c(5e-2, 0.3)), width=1, space=0.25, add=T,
        col=c('white', 'black'),
        xaxt='n', yaxt='n')
#arrows(x0=0.6, x1=2, y0=2, y1=1, col='black', length=0.075, code=2)

mtext('     PSM 1', side=3, line=0.1, cex=0.8, font=2)
points(0.25, 4, pch=25, cex=1.25, xpd=T,bg='white', col='black', lwd=1.5)
mtext('ID Confidence', 2, line=0.35, cex=1, at=-0.5, las=3)

par(#oma=c(1,0,0,0),
    mar=c(0.5, 1, 1.25, 0.5),
    pty='m')

plot(0, 0, type='n',
     xlim=c(0, 2.75), ylim=c(0, 3.5),
     xlab=NA, ylab=NA,
     xaxs='i', yaxs='i',
     xaxt='n', yaxt='n')

barplot(-log10(c(5e-2, 3e-3)), width=1, space=0.25, add=T,
        col=c('white', 'black'),
        xaxt='n', yaxt='n')
#arrows(x0=0.6, x1=1.2, y0=1.5, y1=2.5, col='black', length=0.075, code=2)

legend(x=-0.75, y=0, c('Spectra', 'DART-ID'), 
       pch=22, pt.cex=2, pt.bg=c('white', 'black'), cex=0.9, col='black',
       y.intersp=1, bty='n', xpd=NA)

mtext('     PSM 2', side=3, line=0.1, cex=0.8, font=2)
points(0.25, 3.85, pch=24, cex=1.25, xpd=T, bg='white', col='black', lwd=1.5)

dev.off()
