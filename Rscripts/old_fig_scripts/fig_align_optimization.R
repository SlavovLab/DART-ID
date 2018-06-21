## Optimization Cases ---------

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

pdf(file='manuscript/Figs/Fig_1B2_v1.pdf', width=4, height=3)

m <- rbind(c(1, 1, 1, 1, 1, 1, 4, 4, 4, 4), 
           c(2, 2, 2, 3, 3, 3, 5, 5, 5, 5))
layout(m)

par(mgp=c(0.2, 0.5, 0),
    mar=c(1,2,1.5,0.5),
    #pin=c(2,2),
    #asp=1,
    #pty='s',
    cex.lab=1, cex.axis=1, cex.main=1,
    xaxs='r', yaxs='r')

plot(0, 0,
     xlim=c(24.9, 25.8), ylim=c(23.5, 27.75),
     xlab=NA, ylab=NA,
     xaxt='n', yaxt='n')

#abline(a=0, b=exps[1:3], col='blue')
for(exp in exps) {
  abline(a=0, b=exp, col='blue')
}

x <- rep(mus[9],3)
y <- rts[c(9,29,49)]
y[1] <- y[1] + 0.7
y[3] <- y[3] + 0.7
offset <- 0.68

# segments(x0=mus[9],x1=mus[9],y0=0,y1=27.6,
#          lty=1, lwd=1.5, col=rgb(0,0,0,0.5))
# segments(x0=mus[9]+offset,x1=mus[9]+offset,y0=0,y1=27.6,
#          lty=1, lwd=1.5, col=rgb(0,0,0,0.5))
abline(v=mus[9], lty=1, lwd=1, col=rgb(0,0,0,0.5))
abline(v=mus[9]+offset, lty=1, lwd=1, col=rgb(0,0,0,0.5))

segments(x0=x+offset, x1=x, y0=y, y1=y, col=rgb(0,0,0,0.3), lty=2)

segments(x0=mus[9], x1=mus[9], y0=y, y1=mus[9]*exps, 
         col='red', lwd=2.5)
print(sum((y - mus[9]*exps)^2))
segments(x0=mus[9]+offset, x1=mus[9]+offset, y0=y, y1=(mus[9]+offset)*exps, 
         col='red', lwd=2.5)
print(sum((y - (mus[9]+offset)*exps)^2))

#points(x, y, pch='x', col=c(rep(rgb(0,0,0,0.9),2), 'blue'))
#points(x+offset, y, pch='x', col=c(rep(rgb(0,0,0,0.9),2), 'blue'))
points(x, y, pch='x', cex=1.5, col='blue')
points(x+offset, y, pch='x', cex=1.5, col='blue')

par(xpd=T)
mtext('Shifting Canonical RT', side=3, cex=0.8, line=0.25, adj=0.5)
mtext('Observed RT (min)', side=2, cex=0.8, line=0.5, at=22.5)

#text(x=26.1, y=27, labels="Shifting\nCanonical RT", cex=0.8, adj=c(0, 0.5))
#legend(x=26.025, y=26, c("Residual"), col='red', lty=1, lwd=3, cex=0.75,
#       bg='white', bty='n', inset=c(0, 0))

##
par(mar=c(1.5, 2, 1, 0), xpd=F)
plot(0, 0,
     xlim=c(24.9, 25.8), ylim=c(23.5, 27.75),
     xlab=NA, ylab=NA,
     xaxt='n', yaxt='n')

for(exp in exps) {
  abline(a=0, b=exp, col='blue')
}

x <- rep(c(25.08, 25.3, 25.5), each=3)
y <- rep(exps, 3) * x
y[5:6] <- y[5:6] + 0.4

segments(x0=x, x1=x, y0=y, y1=x*rep(exps, 3), 
         col='red', lwd=2.5)
print(sum((y - x*rep(exps,3))^2))

mtext('Shifting Experiment Transform', side=3, cex=0.8, line=0.35, at=25.95)
points(x, y, pch='x', cex=1.5, col='blue')

##
par(mar=c(1.5, 1.25, 1, 0.5))
plot(0, 0,
     xlim=c(24.9, 25.8), ylim=c(23.5, 27.75),
     xlab=NA, ylab=NA,
     xaxt='n', yaxt='n')

exps_new <- exps
exps_new[2:3] <- exps_new[2:3] * 1.007

for(exp in exps_new) {
  abline(a=0, b=exp, col='blue')
}

y <- rep(exps, 3) * x
y[5:6] <- y[5:6] + 0.4

segments(x0=x, x1=x, y0=y, y1=x*rep(exps_new, 3), 
         col='red', lwd=2.5)
print(sum((y - x*rep(exps_new,3))^2))

points(x, y, pch='x', cex=1.5, col='blue')

##
par(mar=c(1,2,1.5,1))
plot(0, 0, type='n',
     xlab=NA, ylab=NA,
     xlim=c(0, 2.75), ylim=c(0, 1),
     xaxs='i', yaxs='i',
     xaxt='n', yaxt='n')

barplot(c(0.784, 0.374), width=1, space=0.25, add=T,
        xaxt='n', yaxt='n')

axis(side=1, at=c(0.75, 2), tck=-0.02, padj=-0.2, labels=c('a', 'b'))

mtext('Total Residual RT (min)', side=2, cex=0.8, line=0.5, at=-0.1)

##
par(mar=c(1.5,2,1,1))
plot(0, 0, type='n',
     xlab=NA, ylab=NA,
     xlim=c(0, 2.75), ylim=c(0, 1),
     xaxs='i', yaxs='i',
     xaxt='n', yaxt='n')

barplot(c(0.32, 0.223), width=1, space=0.25, add=T,
        xaxt='n', yaxt='n')

axis(side=1, at=c(0.75, 2), tck=-0.02, padj=-0.2, labels=c('a', 'b'))

dev.off()