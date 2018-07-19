library(tidyverse)
source('Rscripts/lib.R')

ev <- read_tsv("/gd/Slavov_Lab/Albert/RTLib_Alignments/SQC_20180621_2/ev_updated.txt")

## --------

ev.f <- ev %>%
  select(c('Sequence', 'PEP', 'pep_new')) %>%
  filter(!is.na(pep_new)) %>%
  filter(PEP > 0 & pep_new > 0 & PEP > 1e-5 & pep_new > 1e-5) %>%
  mutate_at(c('PEP', 'pep_new'), funs(ifelse(. > 1, 1, .))) %>%
  mutate(pep_log=log10(PEP),
         pep_new_log=log10(pep_new))

## --------


nbins <- 50
x.bin <- seq(-5, 0, length=nbins)
y.bin <- seq(-5, 0, length=nbins)

freq <- as.data.frame(table(findInterval(ev.f$pep_log, x.bin),
                            findInterval(ev.f$pep_new_log, y.bin)))
freq[,1] <- as.numeric(freq[,1])
freq[,2] <- as.numeric(freq[,2])

freq2D <- diag(nbins)*0
freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]

colfunc <- colorRampPalette(c('white', 'red'))

## -------

pdf(file='manuscript/Figs/pep_scatter_v4.pdf', width=3.5, height=3)

par(mar=c(2.5,2,1.25,1),
    pty='s', las=1,
    cex.axis=0.85, cex.lab=1, cex.main=1)

cols <- colfunc(40)

# Normal
image(x.bin, y.bin, freq2D, col=cols,
      xlab=NA, ylab=NA,
      xaxs='i', yaxs='i',
      xaxt='n', yaxt='n', useRaster=F)

abline(a=0, b=1, col='black')
#abline(h=-2, col='black', lty=2)
#abline(v=-2, col='black', lty=2)
segments(x0=-2, x1=-2, y0=-5, y1=-2, col='black', lty=2)
segments(x0=-2, x1=0, y0=-2, y1=-2, col='black', lty=2)

text(x=-3.5, y=-1.9, labels='New PSMs selected\nat 0.01 PEP', adj=c(0, 0), cex=0.85)

rng <- seq(-5, 0, 1)
axis(1, tck=-0.02,  
     at=rng, labels=fancy_scientific(10^rng),
     mgp=c(0, 0.2, 0))
axis(2, tck=-0.02, 
     at=rng, labels=fancy_scientific(10^rng),
     mgp=c(0, 0.4, 0), las=1)

mtext('Spectra', 1, line=1, cex=1)
mtext('DART-ID', 2, line=1.85, cex=1, las=3)
mtext('Error Probability (PEP)', 3, line=0.1, cex=1, font=2)

scale = (length(cols)-1)/(-1-(-3))
for (i in 1:(length(cols)-1)) {
  y = (i-1)/scale + -3
  rect(-4.5, y, -4, y+1/scale, col=cols[i], border=NA, add=T)
}
text(x=-4.25, y=-0.85, labels='Density', adj=c(0.5, 0), cex=0.85)

# Log
#image(x.bin, y.bin, log(freq2D), col=r)

dev.off()

