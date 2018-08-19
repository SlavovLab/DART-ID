library(tidyverse)
library(pracma)
source('Rscripts/lib.R')

#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180621_2/ev_updated.txt")
#ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180724_3/ev_updated.txt")
#ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180813_with_PI/ev_updated.txt')
ev <- read_tsv('/gd/bayesian_RT/Alignments/SQC_20180815_2/ev_updated.txt')

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
            min_pep_new=min(qval_updated),
            n=n()) %>%
  dplyr::select(c('Modified sequence', 'min_pep', 'min_pep_new', 'n')) %>%
  filter(min_pep > 0.01 & min_pep_new < 0.01) %>%
  pull(`Modified sequence`)

ev$qval_prev <- ev$qval_updated
ev$qval_prev[ev$`Modified sequence` %in% new_peptides] <- ev$qval[ev$`Modified sequence` %in% new_peptides]

ev$fdr_perc[is.na(ev$fdr_perc)] <- ev$qval[is.na(ev$fdr_perc)]

# same but with percolator
new_peptides <- ev %>%
  #filter(!is.na(pep_new)) %>%
  group_by(`Modified sequence`) %>%
  summarise(min_pep=min(qval),
            min_pep_new=min(fdr_perc)) %>%
  dplyr::select(c('Modified sequence', 'min_pep', 'min_pep_new')) %>%
  filter(min_pep > 0.01 & min_pep_new < 0.01) %>%
  pull(`Modified sequence`)

ev$fdr_perc_prev <- ev$fdr_perc
ev$fdr_perc_prev[ev$`Modified sequence` %in% new_peptides] <- ev$fdr_perc[ev$`Modified sequence` %in% new_peptides]

## Fold change increase in IDs ------

x <- logseq(5e-4, 1, 100)

# frame to hold the results
df <- data.frame()
method.names <- c('Spectra', 'DART-ID', 'DART-ID (conf only)', 'Percolator', 'Percolator (conf only)')
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
    sum(ev$qval_prev < i) /         sum(ev$qval < i),
    sum(ev$fdr_perc < i) / sum(ev$qval < i),
    sum(ev$fdr_per_prev < i) / sum(ev$qval < i)
  )
  ident <- c(
    sum(ev$qval < i) /              nrow(ev),
    sum(ev$qval_updated < i) /      nrow(ev),
    sum(ev$qval_prev < i) /         nrow(ev),
    sum(ev$fdr_perc < i) / nrow(ev),
    sum(ev$fdr_perc_prev < i) / nrow(ev)
  )
  
  df <- rbind(df, data.frame(
    x=as.numeric(i),
    ratio=as.numeric(ratios),
    ident=as.numeric(ident),
    Method=as.character(method.names)
  ))
}
df$Method <- factor(df$Method, levels=c('Spectra', 'Percolator', 'Percolator (conf only)',
                                        'DART-ID', 'DART-ID (conf only)'))


# fold change of IDs at FDR thresh ----------------------------------------

pdf(file='manuscript/Figs/fdr_fold_change_v4.pdf', width=1.75, height=3)

layout(c(1, 2), heights=c(1.8, 2))

par(oma=c(0, 0, 1.5, 0),
    mar=c(0.1,2.15,0,0.25),
    #pty='s', 
    las=1,
    cex.axis=0.65, cex.lab=1, cex.main=1)

plot(0, 0, type='n',
     xlim=c(-3.1, -0.3), ylim=c(-15, 130),
     xlab=NA, ylab=NA, xaxs='i', yaxs='i', xaxt='n', yaxt='n')

abline(v=-2, col='black', lty=2, lwd=1)

methods <- c('Spectra', 'Percolator', 'Percolator (conf only)', 'DART-ID', 'DART-ID (conf only)')
cols <- c(cb[1], cb[3], paste0(cb[3], 'CC'), cb[2], paste0(cb[2], 'CC'))
ltys <- c(1, 1, 2, 1, 2)
for(i in 1:length(methods)) {
  df_a <- df %>% filter(Method == methods[i])
  lines(log10(df_a$x), (df_a$ratio-1)*100, col=cols[i], lty=ltys[i], lwd=2)
}

rng <- seq(-3, 0, 1)
axis(1, tck=-0.01,
     at=rng, 
     #labels=fancy_scientific(10^rng), 
     #labels=c('0.1%', '1%', '10%', '100%'),
     labels=NA,
     mgp=c(0, 0.05, 0))
axis(2, tck=-0.02,
     #labels=NA,
     at=c(-25,seq(0, 125, 25)), 
     mgp=c(0, 0.3, 0))

legend('topright', c('Spectra', 'Percolator', 'DART-ID'),
       lwd=2, lty=1, col=c(cols[1], cols[2], cols[4]), seg.len=0.8,
       bty='n', cex=0.65, x.intersp=0.6, y.intersp=1.2, inset=c(-0.01, -0.02))

#mtext('FDR Threshold', 1, line=1, cex=1)

mtext('% Increase', 2, line=1.4, cex=0.85, las=3)
mtext('      Increase in confident PSMs', 3, line=0.2, cex=0.7, font=2, outer=T)

par(mar=c(1.75, 2.15, 0.1, 0.25))

plot(0, 0, type='n',
     xlim=c(-3.1, -0.3), ylim=c(0.18, 1.02),
     xlab=NA, ylab=NA, xaxs='i', yaxs='i', xaxt='n', yaxt='n')

abline(v=-2, col='black', lty=2)

methods <- c('Spectra', 'Percolator', 'Percolator (conf only)', 'DART-ID', 'DART-ID (conf only)')
#cols <- c(av[1], av[3], av[2])
for(i in 1:length(methods)) {
  df_a <- df %>% filter(Method == methods[i])
  lines(log10(df_a$x), df_a$ident, col=cols[i], lty=ltys[i], lwd=2)
}

rng <- seq(-3, 0, 1)
axis(1, tck=-0.02,
     at=rng, 
     #labels=fancy_scientific(10^rng), 
     labels=c('0.1%', '1%', '10%', '100%'),
     mgp=c(0, -0.05, 0))
axis(2, tck=-0.02,
     at=seq(0, 1, by=0.2),
     labels=seq(0, 1, by=0.2)*100, 
     mgp=c(0, 0.3, 0))

mtext('FDR Threshold', 1, line=0.75, cex=0.75)
mtext('% of all PSMs', 2, line=1.4, cex=0.85, las=3)
#mtext('Confident PSMs Selected', 3, line=0.1, cex=0.75, font=2)

dev.off()
