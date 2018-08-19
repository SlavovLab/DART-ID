library(tidyverse)
source('Rscripts/lib.R')

#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180621_2/ev_updated.txt")
#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180724_3/ev_updated.txt")
ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180813_with_PI/ev_updated.txt')

## Add percolator data --------

source('Rscripts/add_percolator.R')

# compare # of false positives --------------------------------------------

ev.f <- ev %>%
  mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .)))

k <- 500
peps <- logseq(1e-8, 1, k)

fp_mq <- vector(length=k)
fp_dart <- vector(length=k)
fp_perc <- vector(length=k)

for(i in 1:length(peps)) {
  cat('\r', i,'/',length(peps), '       ')
  fp_mq[i] <- ceil(sum(ev.f$PEP[ev.f$PEP < peps[i]]))
  fp_dart[i] <- ceil(sum(ev.f$PEP[ev.f$pep_updated < peps[i]]))
  fp_perc[i] <- ceil(sum(ev.f$pep_perc_updated[ev.f$pep_perc_updated < peps[i]]))
}

# vs. each other ----------------------------------------------------------

#plot(fp_mq, fp_dart, pch=16)
#abline(a=0, b=1, col='red', lwd=2)

# plot --------------------------------------------------------------------

pdf(file='manuscript/Figs/fdr_continuity_v1.pdf', width=1.75, height=1.5)

par(mar=c(1.75, 0.5, 0.25, 0.5), cex.axis=0.65,
    oma=c(0, 1.25, 0, 0))

lwidth <- 2

plot(0, 0, type='n',
     xlim=c(-4, 0.03), ylim=c(0, 1.3e5),
     xaxt='n', yaxt='n', xaxs='i', yaxs='r', xlab=NA, ylab=NA)
     #xlab='log10 PEP Thresold', ylab='log2 false positives',
     #main='Number of false positives')
lines(log10(peps), (fp_dart), lwd=lwidth, col=av[2])
lines(log10(peps), (fp_perc), lwd=lwidth, col=av[3])
lines(log10(peps), (fp_mq),   lwd=lwidth, col=av[1])

par(lheight=0.8)
legend('topleft', c(paste0('Spectra\nN = ', max(fp_mq)),
                    paste0('Percolator\nN = ', max(fp_perc)),
                    paste0('DART-ID\nN = ', max(fp_dart))), 
       col=c(av[1], av[3], av[2]), 
       lwd=2, seg.len=0.75, y.intersp=1.75, x.intersp=0.5,
       cex=0.7, bty='n', inset=c(-0.03, -0.09))

rng <- seq(-10, 0, 1)
axis(1, at=rng, labels=fancy_scientific(10^rng),
     tck=-0.02, mgp=c(0, 0, 0))
#axis(2, at=seq(0, 15, by=5),
#     tck=-0.02, mgp=c(0, 0.25, 0), las=1)
axis(2, at=seq(0, 1.3e5, by=2e4), labels=seq(0, 13, by=2),
     tck=-0.02, mgp=c(0, 0.25, 0), las=1)

mtext('PEP Threshold', side=1, line=0.8, cex=0.75)
mtext('# false positives * 10,000', side=2, line=0.5, cex=0.75, outer=T)


dev.off()