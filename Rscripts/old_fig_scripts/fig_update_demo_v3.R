source('Rscripts/lib.R')

## lets give this another go, and add the alignment step in as well --------

pdf(file='manuscript/Figs/Fig1B.pdf', width=3.25,height=2.5)
par(mgp=c(1.5, 0.5, 0),
    mar=c(2.5,2.5,0.2,0.2),
    cex.lab=0.8,
    cex.axis=0.8)

mu1 <- 39.2
mu2 <- 34.8

set.seed(3)

exps <- 1:10
# generated with rnorm(10, 1, 0.1)
exp_transforms <- c(1.02, 1.03, 0.98, 0.97, 1.01, 0.99, 1.01, 1.02, 0.96, 0.99)
exp_transforms <- exp_transforms + rnorm(10, 0, 0.0005)
# generated with rnorm(10, 1, 0.1)
#error <- c(0.9637004,1.2088793,1.0443407,0.9719422,0.9451481,1.0277092,0.9766545,1.0837655,1.0597011,0.8805878)
error1 <- rnorm(10, 1, 0.015)
error2 <- rnorm(10, 1, 0.015)
#exp_transforms <- exp_transforms * error

ymin = 33
ymax = 43.5

plot(exps, mu1*exp_transforms*error1, pch=16, col='black',
     xlim=c(1, length(exps)), ylim=c(ymin, ymax), xaxt='n',
     xlab="Experiment", ylab="Retention Time (min)", cex=0.75)
points(exps, mu1*exp_transforms*error1-(mu1*exp_transforms-mu1), pch=16, col='blue')
segments(x0=exps, x1=exps, y0=mu1*exp_transforms*error1, y1=mu1*exp_transforms*error1-(mu1*exp_transforms-mu1), col='black', lwd=2)
points(exps+0.2, mu1*exp_transforms*error1-(mu1*exp_transforms-mu1), pch=16, col='blue')
segments(x0=exps+0.2, x1=exps+0.2, y0=mu1*exp_transforms*error1-(mu1*exp_transforms-mu1), y1=mu1, col='red', lwd=2)
abline(a=mu1, b=0, col='blue', lty=2)

points(exps, mu2*exp_transforms*error2, pch=16, col='black')
points(exps, mu2*exp_transforms*error2-(mu2*exp_transforms-mu2), pch=16, col='blue')
segments(x0=exps, x1=exps, y0=mu2*exp_transforms*error2, y1=mu2*exp_transforms*error2-(mu2*exp_transforms-mu2), col='black', lwd=2)
points(exps+0.2, mu2*exp_transforms*error2-(mu2*exp_transforms-mu2), pch=16, col='blue')
segments(x0=exps+0.2, x1=exps+0.2, y0=mu2*exp_transforms*error2-(mu2*exp_transforms-mu2), y1=mu2, col='red', lwd=2)
abline(a=mu2, b=0, col='blue', lty=2)

axis(side=1, at=1:length(exps), cex.axis=1, tck=-0.02)
#axis(side=2, at=seq(ceiling(ymin), floor(ymax), 2), 
#     labels=c(T, F, T, F, T, F, T, F, T, F, T), tck=-0.02)

text(2.25, 40.5, labels="Peptide A", cex=0.8)
text(2.25, 36.2, labels="Peptide B", cex=0.8)

legend("top", 
       c("Observed RT", "Adjusted RT","", "RT Adjustment", "Canonical RT", "Residual (Error)"), 
       pch=c(16,16,16,NA,NA,NA), 
       lty=c(NA,NA,NA,1,2,1),
       lwd=c(NA,NA,NA,2,1,2), 
       pt.cex=c(1,1,0,NA,NA,NA),
       col=c('black', 'blue', 'transparent', 'black', 'blue', 'red'),
       bty='n', cex=0.75, ncol=2)

#par(resetPar())
#par(mgp=c(1.5, 0.5, 0))
dev.off()

## second panel ------

pdf(file='manuscript/Figs/Fig1C.pdf', width=2.75, height=2.5)
par(mgp=c(1.5, 0.5, 0),
    mar=c(2.5,2.5,0.3,0.2),
    cex.lab=0.8,
    cex.axis=0.8)

muij1 <- mu1*exp_transforms*error1-(mu1*exp_transforms-mu1)
res1 <- muij1 - mu1
muij2 <- mu2*exp_transforms*error2-(mu2*exp_transforms-mu2)
res2 <- muij2 - mu2
sigma1 <- sd(muij1) * 0.6
sigma2 <- sd(muij2) * 0.6
x = seq(ymin,ymax-2,length.out=1000)

plot(muij1,rep(0,length(exps)), pch=16, col='blue', 
     xlim=c(ymin, ymax-2), ylim=c(0, 2),
     xlab="Retention Time (min)", ylab="Density")
#abline(v=mu1, col='blue', lty=2)
segments(x0=mu1, x1=mu1, y0=0, y1=dnorm(mu1, mean=mu1, sd=sigma1), col='blue', lty=2)
lines(x[400:1000], dnorm(x[400:1000], mean=mu1, sd=sigma1))

points(muij2,rep(0,length(exps)), pch=16, col='blue')
#abline(v=mu2, col='blue', lty=2)
segments(x0=mu2, x1=mu2, y0=0, y1=dnorm(mu2, mean=mu2, sd=sigma2), col='blue', lty=2)
lines(x[1:400], dnorm(x[1:400], mean=mu2, sd=sigma2))

# null distribution
lines(x, dlnorm(x, meanlog=3.595, sdlog=2.85) * 20, col='red')

text(34.8, 1.75, labels="Peptide B")
text(39.2, 1.75, labels="Peptide A")

# pull out cases
case1=muij2[4]
segments(x0=case1, x1=case1, y0=0, y1=1, col='black', lwd=1.5)
segments(x0=case1, x1=36, y0=1, y1=1, col='black', lwd=1.5)
case2=muij2[9]
segments(x0=case2, x1=case2, y0=0, y1=0.6, col='black', lwd=1.5)
segments(x0=case2, x1=36, y0=0.6, y1=0.6, col='black', lwd=1.5)

text(37.2, 1, labels="Case (ii) ", cex=0.8)
text(37.2, 0.6, labels="Case (i)", cex=0.8)

dev.off()

