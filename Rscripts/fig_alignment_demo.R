pdf(file='manuscript/Figs/Fig_1B_v1.pdf', width=3, height=3)
par(mgp=c(1.2, 0.5, 0),
    mar=c(2.5,2,1.5,0.2),
    pin=c(2, 2),
    cex.lab=0.8, cex.axis=0.8, cex.main=1,
    xaxs='r', yaxs='r',
    pty='s')

set.seed(15)

exps <- c(1.06, 1.01, 0.95)

#mus <- runif(20, min=20,max=30)
mus <- c(20.30265, 20.63603, 20.92263, 22.13038, 22.43004, 23.09638, 23.49773,
         23.96849, 25.01534, 25.65994, 25.95323, 26.40581, 26.94700, 27.44276,
         27.61737, 28.42095, 28.68403, 29.08916, 29.23397, 29.79338)
points_per_exp <- 1
musx <- rep(rep(mus, each=points_per_exp),length(exps))
rts <- musx * rep(exps, each=points_per_exp*length(mus))
error <- rnorm(length(rts), 0, 0.25)
rts <- rts + error
rts[49] <- rts[49] - 0.3

cols <- rep(rgb(0,0,0,0.3), length(musx))
cols[-(1:(2*length(mus)))] <- 'blue'

plot(musx, rts, pch='x', cex=1, col=cols,
     xlim=c(22,30), ylim=c(22,32),
     xaxt='n',yaxt='n',
     xlab="Canonical RT (min)", ylab="Observed RT (min)")
title("Global RT Alignment", line=0.5)

#abline(a=0, b=1, col='black')
abline(a=0, b=exps[3], col='blue')
for(exp in exps) {
  abline(a=0, b=exp, col=rgb(0,0,0,0.3))
}
# lines for mus
abline(v=mus[-9], col=rgb(0,0,0,0.1), lty=2)


# highlight one of the peptides
#rect(xleft=mus[9]-0.4, xright=mus[9]+0.4, ybottom=rts[49]-1, ytop=rts[9]+1,
#     lty=1, lwd=1)
abline(v=mus[9], col=rgb(0,0,0,0.4), lty=2, lwd=1.5)
#text(x=9, y=0, labels="Peptide", cex=0.8, adj=c(0,0))
#text(x=mus, y=rep(0,3), labels=c('I', 'J', 'K'), cex=0.8, font=2, adj=c(0.5, 0))
#text(x=mus[13]+0.6, y=23.8, labels=expression(i), cex=0.8, adj=c(0,0.5))
text(x=mus[9]+0.3, y=rts[49], labels="Inferred RT Distribution\nPeptide X, Experiment A", 
     cex=0.7, adj=c(0,1))

# draw the sideways distribution
polygon(mus[9] - (dnorm(seq(22, 26, length.out=100), mean=(mus[9] * exps[3]), sd=0.5)),
        seq(22, 26, length.out=100), 
        col=rgb(0, 0, 1, 0.3), border=NA)


axis(side=1, tck=-0.02, padj=-0.6)
axis(side=2, tck=-0.02, padj=0.2)

legend('topleft',
       c("Transformed Canonical RT", "Other Experiments", "Experiment A"), 
       pch=c(NA,'x','x'),
       col=c(rgb(0,0,0,0.3), rgb(0,0,0,0.5), "blue"), lty=c(1, 1, 1), lwd=2,
       bty='n', pt.cex=1, cex=0.75, inset=c(0.03,0), y.intersp=1.1)

dev.off()