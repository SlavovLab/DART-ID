library(tidyverse)
source('Rscripts/lib.R')

ev <- read_tsv("/gd/Slavov_Lab/Albert/RTLib_Alignments/SQC_20180621_2/ev_updated.txt")

## Add percolator data --------

source('Rscripts/add_percolator.R')

## qvalues -----
ev <- ev %>%
  # ceil PEPs to 1
  mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .))) %>%
  # calculate q-values
  mutate(qval=(cumsum(PEP[order(PEP)]) / 
                 seq(1, nrow(ev)))[order(order(PEP))],
         qval_updated=(cumsum(pep_updated[order(pep_updated)]) / 
                         seq(1, nrow(ev)))[order(order(pep_updated))])


## Fold change increase in IDs ------

x <- logseq(1e-3, 1, 100)

# frame to hold the results
df <- data.frame()
method.names <- c('Spectra', 'DART-ID', 'Percolator')
counter <- 1
for(i in x) {
  cat('\r', counter,'/',length(x), '       ')
  flush.console()
  counter <- counter + 1
  
  ratios <- c(
    1,
    sum(ev$pep_updated < i) /      sum(ev$PEP < i),
    sum(ev$pep_perc_updated < i) / sum(ev$PEP < i)
  )
  ident <- c(
    sum(ev$qval < i, na.rm=T) /              nrow(ev),
    sum(ev$qval_updated < i, na.rm=T) /      nrow(ev),
    sum(ev$fdr_perc < i, na.rm=T) / sum(!is.na(ev$fdr_perc))
  )
  
  df <- rbind(df, data.frame(
    x=as.numeric(i),
    ratio=as.numeric(ratios),
    ident=as.numeric(ident),
    Method=as.character(method.names)
  ))
}
df$Method <- factor(df$Method, levels=c('Spectra', 'Percolator', 'DART-ID'))

## ---------

pdf(file='manuscript/Figs/pep_fold_change_v3.pdf', width=2.32, height=2.25)

par(mar=c(2,2,1,1),
    pty='s', las=1,
    cex.axis=0.6, cex.lab=0.75, cex.main=0.75)

plot(0, 0, type='n',
     xlim=c(-3, 0), ylim=c(-25, 200),
     xlab=NA, ylab=NA,
     xaxs='i', yaxs='i',
     xaxt='n', yaxt='n')

abline(v=-2, col='black', lty=2)

methods <- c('Spectra', 'Percolator', 'DART-ID')
cols <- c(av[1], av[3], av[2])
for(i in 1:length(methods)) {
  df_a <- df %>% filter(Method == methods[i])
  lines(log10(df_a$x), (df_a$ratio-1)*100, col=cols[i], lwd=2)
}

rng <- seq(-3, 0, 1)
axis(1, tck=-0.02, padj=-2, 
     at=rng, labels=fancy_scientific(10^rng), 
     mgp=c(3, 1, 0))
axis(2, tck=-0.02, padj=0.5,
     at=c(-30,seq(0, 200, 50)), 
     mgp=c(0, 0.3, 0))

legend('topright', c('Spectra', 'Percolator', 'DART-ID'),
       lwd=2, lty=1, col=c(av[1], av[3], av[2]),
       bty='n', cex=0.7, y.intersp=1, inset=c(0.03, 0))

mtext('PEP Threshold', 1, line=0.75, cex=0.75)
mtext('% Increase', 2, line=1.25, cex=0.75, las=3)
mtext('Increase in Confident PSMs', 3, line=0.1, cex=0.75, font=2)

dev.off()

## num psms --------

pdf(file='manuscript/Figs/pep_num_psms_v3.pdf', width=2.32, height=2.25)

par(mar=c(2,2,1,1),
    pty='s', las=1,
    cex.axis=0.6, cex.lab=0.75, cex.main=0.75)

plot(0, 0, type='n',
     xlim=c(-3, 0), ylim=c(0, 1),
     xlab=NA, ylab=NA,
     xaxs='i', yaxs='i',
     xaxt='n', yaxt='n')

abline(v=-2, col='black', lty=2)

methods <- c('Spectra', 'Percolator', 'DART-ID')
cols <- c(av[1], av[3], av[2])
for(i in 1:length(methods)) {
  df_a <- df %>% filter(Method == methods[i])
  lines(log10(df_a$x), df_a$ident, col=cols[i], lwd=2)
}

rng <- seq(-3, 0, 1)
axis(1, tck=-0.02, padj=-2, 
     at=rng, 
     #labels=fancy_scientific(10^rng), 
     labels=c('0.1%', '1%', '10%', '100%'),
     mgp=c(3, 0.65, 0))
axis(2, tck=-0.02, padj=0.5,
     at=seq(0, 1, by=0.2), 
     mgp=c(0, 0.3, 0))

legend('bottomright', c('Spectra', 'Percolator', 'DART-ID'),
       lwd=2, lty=1, col=c(av[1], av[3], av[2]),
       bty='n', cex=0.7, y.intersp=1, inset=c(0.03, 0))

mtext('FDR Threshold', 1, line=0.75, cex=0.75)
mtext('Fraction of all PSMs', 2, line=1.1, cex=0.75, las=3)
mtext('Confident PSMs Selected', 3, line=0.1, cex=0.75, font=2)

dev.off()
