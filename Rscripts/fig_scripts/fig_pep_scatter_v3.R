library(tidyverse)
source('Rscripts/lib.R')

#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180621_2/ev_updated.txt")
#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180724_3/ev_updated.txt")
ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180803_5exp_parametric_mixture_v2/ev_updated.txt")

## --------

conf_limit <- 1e-8

ev.f <- ev %>%
  select(c('Sequence', 'PEP', 'pep_new')) %>%
  filter(!is.na(pep_new)) %>%
  filter(PEP > 0 & pep_new > 0 & PEP > conf_limit & pep_new > conf_limit) %>%
  mutate_at(c('PEP', 'pep_new'), funs(ifelse(. > 1, 1, .))) %>%
  mutate(pep_log=log10(PEP),
         pep_new_log=log10(pep_new))

## --------


nbins <- 80
x.bin <- seq(log10(conf_limit), 0, length=nbins)
y.bin <- seq(log10(conf_limit), 0, length=nbins)

freq <- as.data.frame(table(findInterval(ev.f$pep_log, x.bin),
                            findInterval(ev.f$pep_new_log, y.bin)))
freq[,1] <- as.numeric(freq[,1])
freq[,2] <- as.numeric(freq[,2])

freq2D <- diag(nbins)*0
freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]

colfunc <- colorRampPalette(c('white', 'blue'))

## -------

pdf(file='manuscript/Figs/pep_scatter_v7.pdf', width=3.5, height=3)

layout(t(c(1, 2)), widths=c(5, 1))

par(mar=c(2.25,2.75,1.25,0.75),
    pty='s', las=1,
    cex.axis=0.85, cex.lab=1, cex.main=1)

cols <- colfunc(20)

# Normal
image(x.bin, y.bin, freq2D, col=cols,
      xlab=NA, ylab=NA,
      xaxs='i', yaxs='i',
      xaxt='n', yaxt='n', useRaster=F)

abline(a=0, b=1, col='black')
abline(h=-2, col='black', lty=2, lwd=1)
abline(v=-2, col='black', lty=2, lwd=1)
#segments(x0=-2, x1=-2, y0=log10(conf_limit), y1=-2, col='black', lty=2)
#segments(x0=-2, x1=0, y0=-2, y1=-2, col='black', lty=2)

rect(xleft=-2, xright=0, ybottom=log10(conf_limit), ytop=-2,
     border=NA, col=rgb(1,0,0,0.1))
rect(xleft=log10(conf_limit), xright=-2, ybottom=-2, ytop=0,
     border=NA, col=rgb(0,0,1,0.1))

text(-7, -1, 'Downgraded', cex=1, adj=c(0, 0.5))
text(-1, -4.5, 'Upgraded', cex=1, adj=c(0, 0), srt=270)


#text(x=-3.5, y=-1.9, labels='New PSMs selected\nat 0.01 PEP', adj=c(0, 0), cex=0.85)

#rng <- seq(-5, 0, 1)
rng <- seq(-10, 0, 2)
axis(1, tck=-0.02,  
     at=rng, labels=fancy_scientific(10^rng),
     mgp=c(0, 0.2, 0))
axis(2, tck=-0.02, 
     at=rng, labels=fancy_scientific(10^rng),
     mgp=c(0, 0.4, 0), las=1)

mtext('Spectra', 1, line=1.15, cex=1)
mtext('DART-ID', 2, line=1.85, cex=1, las=3)
#mtext('Error Probability (PEP)', 3, line=0.1, cex=1, font=2)
mtext('PEP - Parametric Bootstrap (Mixture)', 3, line=0.1, cex=1, font=2)


#x0 <- -4.5
#x1 <- -4
# x0 <- par('usr')[1] + 1.25
# x1 <- x0 + 0.75
# y0 <- -7
# y1 <- -4
# scale = (length(cols)-1)/(y1-(y0))
# for (i in 1:(length(cols)-1)) {
#   y = (i-1)/scale + y0
#   rect(x0, y, x1, y+1/scale, col=cols[i], border=NA)
# }
# text(x=mean(c(x0, x1)), y=y1+0.25, labels='Density', adj=c(0.5, 0), cex=0.85)

# Log
#image(x.bin, y.bin, log(freq2D), col=r)

par(mar=c(2.5, 0.5, 2, 1.25), pty='m')
image(matrix(seq(-1, 1, length.out=nbins), ncol=nbins), col=colfunc(nbins),
      xlab=NA, ylab=NA, xaxt='n', yaxt='n')
#axis(4, at=seq(0, 1, length.out=6), labels=seq(0.5, 1, length.out=6), 
#     tck=-0.1, las=1, mgp=c(0, 0.3, 0))
#text(1+0.5, 0.5, 'Density', adj=c(0.5, 0), srt=270, xpd=NA)
mtext('Density', side=3, line=0.1)

dev.off()

