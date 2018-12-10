library(tidyverse)
library(rmutil)
source('Rscripts/lib.R')

# no ggplot ---------------------------------------------------------------

pdf(file='manuscript/Figs/mixture_model_v5.pdf', width=3.15, height=4)

layout(c(1, 2))

par(oma=c(2.5,4,1.5,0.75), cex.axis=0.85, xaxs='i', yaxs='i')

x <- seq(0,60,by=0.1)
#y1 <- dlnorm(x, meanlog=4.663, sdlog=0.5089)
null <- dnorm(x, mean=38, sd=17)
#y <- dnorm(x, mean=12, sd=1.78)
y <- dlaplace(x, m=12, s=1.78)
pep <- 0.5

par(mar=c(0.5,0,0,0))

plot(0, 0, type='n', xlab=NA, ylab=NA,
     xlim=c(0, 60), ylim=c(0, 0.30),
     xaxt='n', yaxt='n')

lines(x, y*(1-pep), type='l', col='blueviolet', lwd=2)
polygon(c(x, max(x), min(x)), c(y*(1-pep), 0, 0), 
        border=NA, col=rgb(138,43,226,alpha=150, max=255))
lines(x, null*pep, type='l', col='steelblue', lwd=2)
polygon(c(x, max(x), min(x)), c(null*pep, 0, 0), 
        border=NA, col=rgb(70,130,180,alpha=150,max=255))

points(rep(24, 2), c(0.07, 0.20), pch=15, cex=3,
       col=c(rgb(70,130,180,alpha=150,max=255), rgb(138,43,226,alpha=150, max=255)))

text(20, 0.27, 'Alignment Process', cex=1, adj=c(0, 0.5))

axis(1, at=seq(0, 60, by=10), tck=-0.02, mgp=c(0, 0.2, 0), labels=NA)
axis(2, at=seq(0, 0.3, by=0.1), tck=-0.02, las=1, mgp=c(0, 0.4, 0))

par(mar=c(0,0,0.5,0))

plot(0, 0, type='n', xlab=NA, ylab=NA,
     xlim=c(0, 60), ylim=c(0, 0.30),
     xaxt='n', yaxt='n')
#lines(x, (y*(1-pep)) + (null*pep), col='blue', lwd=2)
lines(x, y, col='blue', lwd=2)
polygon(c(x, max(x), min(x)), c(y, 0, 0), 
        border=NA, col=rgb(0,0,1,0.3))
lines(x, null, col='red', lwd=2)
polygon(c(x, max(x), min(x)), c(null, 0, 0), border=NA, col=rgb(1,0,0,0.3))

points(rep(24, 2), c(0.10, 0.17), pch=15, cex=3,
       col=c(rgb(255,0,0,alpha=150,max=255), rgb(0,0,255,alpha=150, max=255)))

text(25, 0.27, 'Update Process', cex=1, adj=c(0, 0.5))

axis(1, at=seq(0, 60, by=10), tck=-0.02, mgp=c(0, 0.2, 0))
axis(2, at=seq(0, 0.3, by=0.1), tck=-0.02, las=1, mgp=c(0, 0.4, 0))


mtext('Retention time (min)', side=1, outer=T, cex=1, line=1.25)
mtext('Density', side=2, outer=T, cex=1, line=2)

mtext('Example Peptide, PEP = 0.5', side=3, outer=T, cex=1, font=2, line=0.25)

dev.off()