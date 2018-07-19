library(tidyverse)
library(pracma)
source('Rscripts/lib.R')

ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180621_2/ev_updated.txt")

## Add percolator data --------

source('Rscripts/add_percolator.R')


# calculate FDR -----------------------------------------------------------

ev <- ev %>%
  # ceil PEPs to 1
  mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .))) %>%
  # calculate q-values
  mutate(qval=(cumsum(PEP[order(PEP)]) / 
                 seq(1, nrow(ev)))[order(order(PEP))],
         qval_updated=(cumsum(pep_updated[order(pep_updated)]) / 
                         seq(1, nrow(ev)))[order(order(pep_updated))])

# flag peptides that don't have a single confident ID across all sets
new_peptides <- ev %>%
  #filter(!is.na(pep_new)) %>%
  group_by(`Modified sequence`) %>%
  summarise(min_pep=min(qval),
            min_pep_new=min(qval_updated)) %>%
  select(c('Modified sequence', 'min_pep', 'min_pep_new')) %>%
  filter(min_pep > 0.01 & min_pep_new < 0.01) %>%
  pull(`Modified sequence`)

ev <- ev %>%
  mutate(new_id=`Modified sequence` %in% new_peptides)


## Fold change increase in IDs ------

x <- logseq(1e-3, 1, 100)

# frame to hold the results
df <- data.frame()
method.names <- c('Spectra', 'DART-ID', 'DART-ID (conf only)', 'Percolator')
counter <- 1
for(i in x) {
  cat('\r', counter,'/',length(x), '       ')
  flush.console()
  counter <- counter + 1
  
  ratios <- c(
    1,
    #sum(ev$pep_updated < i) /      sum(ev$PEP < i),
    #sum(ev$pep_perc_updated < i) / sum(ev$PEP < i)
    sum(ev$qval_updated < i) /      sum(ev$qval < i),
    sum(ev$qval_updated < i & !ev$new_id) / sum(ev$qval < i & !ev$new_id),
    sum(ev$fdr_perc[!is.na(ev$fdr_perc)] < i) / sum(ev$qval[!is.na(ev$fdr_perc)] < i)
  )
  ident <- c(
    sum(ev$qval < i, na.rm=T) /              nrow(ev),
    sum(ev$qval_updated < i, na.rm=T) /      nrow(ev),
    sum(ev$qval_updated < i & !ev$new_id, na.rm=T) /      nrow(ev),
    sum(ev$fdr_perc < i, na.rm=T) / sum(!is.na(ev$fdr_perc))
  )
  
  df <- rbind(df, data.frame(
    x=as.numeric(i),
    ratio=as.numeric(ratios),
    ident=as.numeric(ident),
    Method=as.character(method.names)
  ))
}
df$Method <- factor(df$Method, levels=c('Spectra', 'Percolator', 
                                        'DART-ID', 'DART-ID (conf only)'))


# fold change of IDs at FDR thresh ----------------------------------------

pdf(file='manuscript/Figs/pep_fold_change_v5.pdf', width=2.5, height=3)

par(mar=c(2,2.15,1,0.25),
    #pty='s', 
    las=1,
    cex.axis=0.75, cex.lab=1, cex.main=1)

plot(0, 0, type='n',
     xlim=c(-3, -0.5), ylim=c(-15, 130),
     xlab=NA, ylab=NA,
     xaxs='i', yaxs='i',
     xaxt='n', yaxt='n')

abline(v=-2, col='black', lty=2, lwd=2)

methods <- c('Spectra', 'Percolator', 'DART-ID', 'DART-ID (conf only)')
cols <- c(av[1], av[3], av[2], av[2])
ltys <- c(1, 1, 1, 3)
for(i in 1:length(methods)) {
  df_a <- df %>% filter(Method == methods[i])
  lines(log10(df_a$x), (df_a$ratio-1)*100, col=cols[i], lty=ltys[i], lwd=2)
}

rng <- seq(-3, 0, 1)
axis(1, tck=-0.02,
     at=rng, 
     #labels=fancy_scientific(10^rng), 
     labels=c('0.1%', '1%', '10%', '100%'),
     mgp=c(0, 0.05, 0))
axis(2, tck=-0.02,
     at=c(-25,seq(0, 125, 25)), 
     mgp=c(0, 0.3, 0))

legend('topright', c('Spectra', 'Percolator', 
                     'DART-ID', 'DART-ID\nPrev IDs only'),
       lwd=2, lty=ltys, col=cols, seg.len=1,
       bty='n', cex=0.75, y.intersp=c(0.6,0.6,0.6,0.75), inset=c(0, -0.03))

mtext('FDR Threshold', 1, line=1, cex=1)
mtext('% Increase', 2, line=1.25, cex=1, las=3)
mtext('Confident PSMs', 3, line=0.1, cex=0.9, font=2)

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
