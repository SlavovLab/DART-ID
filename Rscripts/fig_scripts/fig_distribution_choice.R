library(tidyverse)
library(rmutil)
source('Rscripts/lib.R')

#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180724_3/ev_updated.txt")
#ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180813_with_PI/ev_updated.txt')
ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180815_2/ev_updated.txt')

ev.f <- ev %>%
  filter(!is.na(pep_new)) %>%
  mutate(dRT=`Retention time`-muij)

# laplace dist -----------------------------------------------------------------------

pdf(file='manuscript/Figs/dist_laplace_v2.pdf', width=3.5, height=3.5)

par(mar=c(2.5,3,2,1), cex.axis=1)

hist(ev.f$dRT[abs(ev.f$dRT) <= 2], breaks=seq(-2, 2, length.out=100), freq=F,
     xlab=NA, ylab=NA, main=NA, xaxt='n', yaxt='n')

x <- seq(-2, 2, length.out=1000)
y <- dlaplace(x, m=median(ev.f$dRT), s=mad(ev.f$dRT))
lines(x, y, col='red', lwd=2)

axis(1, at=seq(-2, 2, by=0.5), mgp=c(0, 0.25, 0), tck=-0.02)
axis(2, at=seq(0, 4, by=1), mgp=c(0, 0.4, 0), tck=-0.02, las=1)

mtext('Residual RT (min)', side=1, cex=1, line=1.3)
mtext('Density', side=2, cex=1, line=1.25)
mtext('Laplace Distribution, All PSMs', side=3, cex=1, font=2, line=0.25)

dev.off()


# normal dist -------------------------------------------------------------

rts <- ev %>% 
  #filter(!is.na(pep_new)) %>% 
  filter(`Retention time` < 60) %>% 
  pull(`Retention time`)

pdf(file='manuscript/Figs/dist_normal.pdf', width=3.5, height=3.5)

par(mar=c(2.5,4,2,1), cex.axis=1)

hist(rts, breaks=seq(0, 60, length.out=50), freq=F,
  xlab=NA, ylab=NA, main=NA, xaxt='n', yaxt='n')
axis(1, at=seq(0, 60, by=10), mgp=c(0, 0.25, 0), tck=-0.02)
axis(2, at=seq(0, 0.025, by=0.005), mgp=c(0, 0.4, 0), tck=-0.02, las=1)

mtext('Observed RT (min)', side=1, cex=1, line=1.3)
mtext('Density', side=2, cex=1, line=2.75)
mtext('Normal Distribution, All PSMs', side=3, cex=1, font=2, line=0.25)

dev.off()