# no ggplot ---------------------------------------------------------------

pdf(file='manuscript/Figs/mixture_model_v2.pdf', width=4, height=4)

layout(c(1, 2))

par(oma=c(2,4,0.5,0.5), cex.axis=0.85, xaxs='i', yaxs='i')

x <- seq(0,60,by=0.1)
#y1 <- dlnorm(x, meanlog=4.663, sdlog=0.5089)
null <- dnorm(x, mean=38, sd=17)
y <- dnorm(x, mean=25, sd=1.78)
pep <- 0.5

par(mar=c(0.25,0,0,0))

plot(0, 0, type='n', xlab=NA, ylab=NA,
     xlim=c(0, 60), ylim=c(0, 0.18),
     xaxt='n', yaxt='n')

lines(x, y*(1-pep), type='l', col='blue', lwd=2)
polygon(c(x, max(x), min(x)), c(y*(1-pep), 0, 0), border=NA, col=rgb(0,0,1,0.3))
lines(x, null*pep, type='l', col='purple', lwd=2)
polygon(c(x, max(x), min(x)), c(null*pep, 0, 0), border=NA, col=rgb(1,0,1,0.3))

legend('topleft', parse(text=c('d[ijk]', 'd[k0]')), pch=15,
       col=c(rgb(0,0,1,0.5), rgb(1,0,1,0.5)), bty='n',
       cex=1.25, pt.cex=2.5, y.intersp=1, x.intersp=1.1)

axis(1, at=seq(0, 60, by=10), tck=-0.02, mgp=c(0, 0.2, 0))
axis(2, at=seq(0, 0.2, by=0.05), tck=-0.02, las=1, mgp=c(0, 0.4, 0))



mtext('Retention time (min)', side=1, outer=T, cex=1, line=1)
mtext('Density', side=2, outer=T, cex=1, line=2)

dev.off()