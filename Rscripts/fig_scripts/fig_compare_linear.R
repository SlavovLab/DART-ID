library(tidyverse)
source('Rscripts/lib.R')

# load data ---------------------------------------------------------------

ev_linear <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_linear_20180711_2/ev_updated.txt')
ev_twopiece <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180711_4/ev_updated.txt')

ev_linear_params <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_linear_20180711_2/exp_params.txt')
ev_twopiece_params <- read_tsv('/gd/bayesian_RT/Alignments/SQC_varied_20180711_4/exp_params.txt')

# compare selected alignment ----------------------------------------------

exp <- 4
ev_linear_f <- ev_linear %>% filter(exp_id == exp)
ev_twopiece_f <- ev_twopiece %>% filter(exp_id == exp)

pdf(file='manuscript/Figs/linear_fit_mu_v_rt_v2.pdf', width=3.5, height=6)

layout(rbind(c(1), c(2)))

par(oma=c(2.25, 1.25, 1, 0),
    mar=c(0.5, 0, 0.5, 0),
    pty='s')

plot(ev_linear_f$mu, ev_linear_f$`Retention time`, pch=16, cex=0.75,
     xlim=c(10, 60), ylim=c(10, 60),
     col=rgb(0,0,0,0.5), xlab=NA, ylab=NA, xaxt='n', yaxt='n')
abline(a=ev_linear_params$beta_0[exp], b=ev_linear_params$beta_1[exp],
       col='green')

text(10, 58, 'Linear Model', adj=c(0, 0.5), cex=1, font=2)

axis(1, at=seq(10, 60, by=10), labels=NA, tck=-0.02)
axis(2, at=seq(10, 60, by=10), tck=-0.02,
     mgp=c(0, 0.5, 0), las=1)

plot(ev_twopiece_f$mu, ev_twopiece_f$`Retention time`, pch=16, cex=0.75,
     xlim=c(10, 60), ylim=c(10, 60),
     col=rgb(0,0,0,0.5), xlab=NA, ylab=NA, xaxt='n', yaxt='n')
abline(v=ev_twopiece_params$split_point[exp], col="blue", lty=2)
segments(x0=0, 
         y0=ev_twopiece_params$beta_0[exp], 
         x1=ev_twopiece_params$split_point[exp], 
         y1=ev_twopiece_params$beta_0[exp] + 
           (ev_twopiece_params$beta_1[exp]*ev_twopiece_params$split_point[exp]), 
         col="green", lwd=1.5)
segments(x0=ev_twopiece_params$split_point[exp], 
         y0=ev_twopiece_params$beta_0[exp] + 
           (ev_twopiece_params$beta_1[exp]*ev_twopiece_params$split_point[exp]), 
         x1=400, 
         y1=ev_twopiece_params$beta_0[exp] + 
           (ev_twopiece_params$beta_1[exp]*ev_twopiece_params$split_point[exp]) + 
           (ev_twopiece_params$beta_2[exp]*(400-ev_twopiece_params$split_point[exp])), 
         col="red", lwd=1.5)

text(10, 58, 'Segmented Model', adj=c(0, 0.5), cex=1, font=2)

axis(1, at=seq(10, 60, by=10), tck=-0.02, mgp=c(0, 0.3, 0))
axis(2, at=seq(10, 60, by=10), tck=-0.02,
     mgp=c(0, 0.5, 0), las=1)

mtext('Canonical RT (min)', side=1, outer=T, cex=1, line=1)
mtext('Observed RT (min)', side=2, outer=T, cex=1, line=0.1)

dev.off()


# compare ECDF of residuals -----------------------------------------------

ecdf_linear <- ecdf(ev_linear %>%
  filter(!is.na(muij)) %>%
  mutate(dRT=`Retention time`-muij) %>% pull(dRT))
ecdf_twopiece <- ecdf(ev_twopiece %>%
  filter(!is.na(muij)) %>%
  mutate(dRT=`Retention time`-muij) %>% pull(dRT))

pdf(file='manuscript/Figs/linear_fit_ecdf.pdf', width=3.5, height=3)

par(mar=c(2.5, 2.75, 1.5, 0.75), cex.axis=0.85)

x <- seq(-1, 1, by=0.01)
plot(x, ecdf_linear(x), type='l', col='red', lwd=2,
     xlab=NA, ylab=NA, xaxt='n', yaxt='n',
     ylim=c(0, 1))
lines(x, ecdf_twopiece(x), type='l', col='blue', lwd=2, lty=1)

axis(1, at=seq(-1, 1, by=0.25), tck=-0.02, mgp=c(0, 0.3, 0))
axis(2, at=seq(0, 1, by=0.2), tck=-0.02, mgp=c(0, 0.5, 0), las=1)

legend('topleft', c('Linear', 'Segmented'), col=c('red', 'blue'), lty=1, lwd=2,
       bty='n', inset=c(0, 0), y.intersp=1, cex=0.85)

mtext('Residual RT (min)', side=1, line=1.5, cex=1)
mtext('Density', side=2, line=1.75, cex=1)
mtext('Fit Residual ECDF', side=3, line=0.25, font=2, cex=1)

dev.off()


# compare sigmaijs --------------------------------------------------------

pdf(file='manuscript/Figs/linear_fit_sigmas.pdf', widt=3.5, height=3)

par(oma=c(2, 2, 1, 2), pty='s',
    mar=c(0,0,0,0))

layout(rbind(c(2, 4),
             c(1, 3)),
       widths=c(4, 1),
       heights=c(1, 3))

set.seed(1)
n <- sample.int(nrow(ev_linear), 1e4)

par(mar=c(0.5, 0.5, 0.5, 0), cex.axis=0.75)

plot(ev_linear$sigmaij[n], ev_twopiece$sigmaij[n], 
     pch=16, col=rgb(0,0,0,0.3), cex=0.5,
     xlim=c(0, 1), ylim=c(0, 1),
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
abline(a=0, b=1, col='red')

axis(1, at=seq(0, 1, by=0.2), mgp=c(0, 0.1, 0), tck=-0.02)
axis(2, at=seq(0, 1, by=0.2), mgp=c(0, 0.4, 0), tck=-0.02, las=1)

mtext(bquote(sigma[.('ij')]*.(' - Linear Model')), side=1, line=1.5, cex=0.85)
mtext(bquote(sigma[.('ij')]*.(' - Segmented Model')), side=2, line=1.5, cex=0.85)

par(pty='m', mar=c(0, 1.9, 0, 1.5))
sigmas <- ev_linear$sigmaij[n]
hist(sigmas[sigmas < 1], breaks=seq(0, 1, by=0.025), freq=F,
     xlab=NA, ylab=NA, main=NA, xaxt='n', yaxt='n', ylim=c(0, 10))
axis(1, at=seq(0, 1, by=0.1), labels=NA, tck=-0.03)
axis(2, at=seq(0, 10, by=2), labels=NA, tck=-0.03)

mtext('Freq', side=2, line=1.75, cex=0.85)

par(mar=c(0.5, 0, 0.4, 0))
sigmas <- ev_twopiece$sigmaij[n]
s <- hist(sigmas[sigmas < 1], breaks=seq(0, 1, by=0.025), plot=F)
barplot(s$density, width=1, space=0, col='white', horiz=T,
        xlim=c(-0.5, 10), 
        xlab=NA, ylab=NA, main=NA, xaxt='n', yaxt='n')
axis(1, at=seq(0, 10, by=2), labels=NA, tck=-0.03)
axis(2, at=seq(0, 40, by=4), labels=NA, tck=-0.03)

mtext('Freq', side=1, line=1.25, cex=0.85)

dev.off()

# compare confidence shifts -----------------------------------------------

dpep_linear <- ev_linear %>%
  filter(PEP < 0.01) %>%
  filter(!is.na(pep_new)) %>%
  mutate(dPEP=log10(pep_new/PEP)) %>%
  pull(dPEP)
dpep_twopiece <- ev_twopiece %>%
  filter(PEP < 0.01) %>%
  filter(!is.na(pep_new)) %>%
  mutate(dPEP=log10(pep_new/PEP)) %>%
  pull(dPEP)

#plot(dpep_linear, dpep_twopiece, pch=16, col=rgb(0,0,0,0.1),
#     xlim=c(-5, 5), ylim=c(-5, 5))
#abline(a=0, b=1, col='red')
#abline(h=0, col='black', lty=2)
#abline(v=0, col='black', lty=2)

plot(density(dpep_linear, adjust=3), col='black')
lines(density(dpep_twopiece, adjust=3), col='red')

# compare residual histograms -------------------------------------------------------

ev_linear %>%
  filter(!is.na(muij)) %>%
  mutate(dRT=`Retention time`-muij) %>%
ggplot() +
  geom_histogram(aes(dRT), bins=50) +
  scale_x_continuous(limits=c(-2, 2))

ev_twopiece %>%
  filter(!is.na(muij)) %>%
  mutate(dRT=`Retention time`-muij) %>%
  ggplot() +
  geom_histogram(aes(dRT), bins=50) +
  scale_x_continuous(limits=c(-2, 2))
