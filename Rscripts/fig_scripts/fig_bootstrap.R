library(tidyverse)
library(rmutil)
source('Rscripts/lib.R')


# generate data -----------------------------------------------------------

set.seed(6)

n <- 4

ref_rts <- rlaplace(n, m=30, s=0.1)
# row 1 = shift, row 2 = scale
exp_weights <- rbind(runif(n, -3, 3),runif(n, 0.75, 1.25))
exp_sigmas <- abs(rnorm(n, 0.15, 0.1))
peps <- runif(n, 0.01, 0.5)
abs_rts <- (ref_rts * exp_weights[2,]) + exp_weights[1,]

null_means <- rnorm(n, 30, 1)
null_sds <- rnorm(n, 10, 0.5)

# panel a

pdf(file='manuscript/Figs/bootstrap_demo_v3.pdf', width=7, height=5)

par(mar=c(0,0,0,0),
    oma=c(0,3.5,0,1.5),
    cex.axis=1.2)

layout(rbind(
  c(1,3),
  c(2,3),
  c(4,5),
  c(4,5)))

ip_space <- 3.5
xmin <- 29.7
xmax <- 30.4
ymin <- -0.5
ptcex <- 1.75
ptlwd <- 2.5

x <- seq(28, 32, length.out=1000)
x_null <- seq(0, 60, length.out=1000)
#mu <- median(ref_rts)
y <- sapply(1:n, function(i) {
  dlaplace(x, m=ref_rts[i], s=exp_sigmas[i]) * (1-peps[i])
})
y_null <- sapply(1:n, function(i) {
  dnorm(x_null, mean=null_means[i], sd=null_sds[i]) * (peps[i])
})

num_cols <- 20
colfunc <- colorRampPalette(c('blue', 'red'))
cols <- colfunc(num_cols)
col_lim <- seq(0, 0.3, length.out=num_cols)

# panel b1

par(mar=c(1, 2, ip_space, 6))

plot(0, 0, type='n', xlim=c(xmin, xmax), ylim=c(ymin, max(y)+1), 
     xaxt='n', yaxt='n', xlab=NA, ylab=NA)
abline(h=0, col=rgb(0,0,0,0.5))
#abline(v=ref_rts, lwd=2, lty=2, col=rgb(0,0,0,0.3))
for(i in 1:n) {
  lines(x, y[,i], col=paste0(cols[findInterval(exp_sigmas[i], col_lim)], 'AA'), lwd=2)
  polygon(x, y[,i], col=paste0(cols[findInterval(exp_sigmas[i], col_lim)], '22'), border=NA)
}
points(ref_rts, rep(0, n), pch=4, cex=ptcex, lwd=ptlwd, col=rgb(0,0,1,1))

text(30.15, 10, expression(delta[' '*ik]*' = '*1), adj=c(0, 0.5), cex=1.5)

axis(1, at=seq(29, 31, by=0.1), mgp=c(0, 0.5, 0), tck=-0.05)
axis(2, at=seq(0, max(y)+2, by=4), las=1, mgp=c(0, 0.6, 0), tck=-0.05)

mtext('Density', side=2, line=3, at=-2, cex=1)
mtext('Conditional densities', side=3, line=0.25, cex=1, font=1)
mtext(expression(''%*% (1-lambda[' '*ik])), 
      side=4, las=1, cex=1, xpd=T, line=0.2)

# panel b2

par(mar=c(ip_space, 2, 1, 6))

plot(0, 0, type='n', xlim=c(5, 55), ylim=c(-0.0025, max(y_null)+0.012), 
     xaxt='n', yaxt='n', xlab=NA, ylab=NA)

abline(h=0, col=rgb(0,0,0,0.5))

for(i in 1:n) {
  lines(x_null, y_null[,i], col=paste0(cols[findInterval(exp_sigmas[i], col_lim)], 'AA'),
        lwd=2)
  polygon(x_null, y_null[,i], col=paste0(cols[findInterval(exp_sigmas[i], col_lim)], '22'),
          border=NA)
}
#points(ref_rts, rep(0, n), pch=4, cex=ptcex, lwd=ptlwd, col=rgb(0,0,1,0.5))

text(37, 0.018, expression(delta[' '*ik]*' = '*0), adj=c(0, 0), cex=1.5)

axis(1, at=seq(10, 60, by=10), mgp=c(0, 0.5, 0), tck=-0.05)
axis(2, at=seq(0, max(y_null)+0.02, by=0.01), las=1, mgp=c(0, 0.5, 0), tck=-0.05)

mtext('Aligned RT (min)', side=1, line=2, cex=1)
mtext(expression(''%*%(lambda[' '*ik])), 
      side=4, las=1, cex=1, xpd=T, line=0.2)

mtext(expression(''%->%''), side=4, cex=2, xpd=T, las=1, line=1, at=0.035)
mtext('Combine', side=4, cex=1, xpd=T, las=1, line=4.5, at=0.036)

# panel c

par(mar=c(ip_space, 6, ip_space, 0.5))

plot(0, 0, type='n', xlim=c(29.6, xmax), ylim=c(ymin, 5.5), 
     xaxt='n', yaxt='n', xlab=NA, ylab=NA)

# regenerate both on the same x scale
x <- seq(0, 60, length.out=30000)
#mu <- median(ref_rts)
y <- sapply(1:n, function(i) {
  dlaplace(x, m=ref_rts[i], s=exp_sigmas[i]) * (1-peps[i])
})
y_null <- sapply(1:n, function(i) {
  dnorm(x, mean=null_means[i], sd=null_sds[i]) * (peps[i])
})

y_all <- (apply(y, 1, sum) / n) + (apply(y_null, 1, sum) / n)

# only take between 28 and 32
x <- x[(30000/60*28):(30000/60*32)]
y_all <- y_all[(30000/60*28):(30000/60*32)]

lines(x, y_all, col=cols[floor(num_cols/3)], lwd=2)
polygon(x, y_all, col=paste0(cols[floor(num_cols/3)], '44'), border=NULL)
#abline(v=mu, lwd=2, lty=2, col=rgb(0,0,0))

# sample random points from x
#x_sample <- rlaplace(n, m=ref_rts, s=exp_sigmas)
x_sample <- c(29.72, 29.93, 30.05, 30.13)

points(x_sample, rep(0, n), pch=4, cex=ptcex, lwd=2,
       col=rgb(1,0,0,1))
abline(v=median(x_sample), col=rgb(1,0,0), lty=2, lwd=2)
points(median(x_sample), 0, pch=1, cex=ptcex, lwd=2,
       col=rgb(1,0,0))

axis(1, at=seq(29, 31, by=0.1), mgp=c(0, 0.6, 0), tck=-0.02)
axis(2, at=seq(0, max(y)+2, by=1), las=1, mgp=c(0, 0.6, 0), tck=-0.02)

mtext('Aligned RT (min)', side=1, line=2, cex=1)
#mtext('Density', side=2, line=1.5, cex=0.85)
mtext('Sample from distribution', side=3, line=0.25, cex=1, font=1)

# panel d

par(mar=c(ip_space, 1, ip_space, 0.5))

# sample over k iterations
k <- 8
x_bootstrap <- sapply(1:k, function(i) {
  median(rlaplace(n, m=ref_rts, s=exp_sigmas))
})

plot(0, 0, type='n', xlim=c(29.85, 30.3), ylim=c(ymin, max(y_all)+1), 
     xaxt='n', yaxt='n', xlab=NA, ylab=NA)

abline(h=0)
points(x_bootstrap, rep(0, k), col=rgb(1,0,0), pch=1, cex=ptcex, lwd=ptlwd)

y_bootstrap <- sapply(1:k, function(i) {
  dlaplace(x, m=x_bootstrap[i], s=0.1)
})

for(i in 1:k) {
  y <- y_bootstrap[,i]
  #abline(v=x_bootstrap[i], col=rgb(1,0,0,0.25), lwd=1, lty=2)
  lines(x, y, col=rgb(1,0,0,0.5), lwd=2)
  polygon(x, y, col=rgb(1,0,0,0.05), border=NA)
}

axis(1, at=seq(29, 31, by=0.1), labels=seq(29, 31, by=0.1)+5,
     mgp=c(0, 0.5, 0), tck=-0.02)
axis(2, at=seq(0, max(y_all)+4, by=1), las=1, mgp=c(0, 0.4, 0), tck=-0.02)

mtext('Observed RT (min)', side=1, line=2, cex=1, at=30.35)
mtext('Density', side=2, line=1.5, cex=1)
mtext(expression('Combine '*mu[' '*i]*' samples'), side=3, line=0., cex=1, font=2)

# panel e

par(mar=c(ip_space, 0.5, ip_space, 1))

y_bootstrap_all <- apply(y_bootstrap, 1, sum) / k

plot(0, 0, type='n', xlim=c(xmin, xmax), ylim=c(ymin, max(y_all)+1), 
     xaxt='n', yaxt='n', xlab=NA, ylab=NA)

abline(h=0)
lines(x, y_bootstrap_all, col=rgb(1,0,0), lwd=2)
polygon(x, y_bootstrap_all, border=NA, col=rgb(1,0,0,0.25))

#abline(v=30.12, lwd=2, lty=2)
points(30.12, approx(x, y_bootstrap_all, 30.12)$y, pch=1, cex=ptcex*2, lwd=ptlwd*1.5)

#text(29.71, 4.8, expression('Evaluate '*P(rho[' '*i*' '==' '*alpha*', '*k*' = '*A]*' | '*delta[' '*i*' '==' '*alpha*', '*k*' = '*A] == 1)), cex=1.5, adj=c(0, 0.5))
text(29.83, 4.8, expression('Evaluate '*P(rho[' '*ik]*' | '*delta[' '*ik] == 1)), cex=1.5, adj=c(0, 0.5))

axis(1, at=seq(29, 31, by=0.2), labels=seq(29, 31, by=0.2)+5,
     mgp=c(0, 0.5, 0), tck=-0.02)
axis(2, at=seq(0, max(y_all)+4, by=1), labels=NA, las=1, mgp=c(0, 0.4, 0), tck=-0.02)

mtext('Posterior predictive distribution', side=3, line=0.25, cex=1)

dev.off()
