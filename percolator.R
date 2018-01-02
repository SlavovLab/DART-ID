library(tidyverse)

#library(devtools)
#install_github("ypriverol/pIR")
#library(pIR)

#msms <- read_tsv('dat/msms.txt')
ev <- read_tsv('dat/evidence.txt')

a <- ev[1:1e3,]

# cap protein IDs at 5
max.prots <- 5
prots <- ev$`Leading proteins`
# only take first [max.prots] proteins.
# if less than the max, pad with NAs
prots <- lapply(prots, function(prot.group) {
  .prots <- unlist(strsplit(prot.group, ';'))
  if(length(.prots) >= max.prots) return(.prots[1:max.prots])
  else {
    .prots[(length(.prots)+1):max.prots] <- NA
    return(.prots)
  }
})
# transform into matrix
prots <- do.call(rbind, c(deparse.level=1, prots))

# only use the first one?
#prots.all <- prots
#prots <- prots[,1]

# name the columns
#colnames(prot.ids) <- paste0('proteinId', 1:max.prot.ids)

#seqs <- unique(a$Sequence)
#seqs <- unique(ev$Sequence)
#pi <- as.numeric(sapply(seqs, function(seq) {
#  print(match(seq, seqs))
#  pIIterative(sequence=seq, pkSetMethod='solomon')
#}))

H.mass <- 1.007825035

pin <- ev %>%
  mutate(
    # 1 for target PSMs, -1 for decoys
    label=ifelse(is.na(Reverse), 1, -1),
    # Append dummy alanines on both terminuses
    # because percolator will complain if they're not there
    # not completely certain about how they're used
    Sequence=paste0('A.', Sequence, '.A'),
    pI=0) %>%
  rename(dRT='Retention time calibration',
         #dM='Uncalibrated - Calibrated m/z [Da]',
         #dM='Mass Error [Da]',
         Scan.number='MS/MS Scan Number') %>%
  mutate(CalcMass=Mass,
         ExpMass=((`m/z` * Charge) - (Charge * H.mass)),
         dM=CalcMass-ExpMass, absdM=abs(dM),
         dRT=abs(dRT),
         Charge1=as.numeric(Charge==1), Charge2=as.numeric(Charge==2), 
         Charge3=as.numeric(Charge==3), Charge4=as.numeric(Charge==4), 
         Charge5=as.numeric(Charge==5)) %>%
  mutate(dMdRT=dM * dRT) %>%
  select(id, label, Scan.number, ExpMass, CalcMass, `Retention time`, dM, 
         Score, `Delta score`, `Uncalibrated - Calibrated m/z [Da]`, Length,
         `Retention length`, dRT, PIF, `Number of isotopic peaks`,
         Charge1, Charge2, Charge3, Charge4, Charge5,
         CalcMass, dM, absdM, Sequence)
  
pin <- cbind(pin, prots)

pin <- pin %>%
  filter(!is.na(dM) & !is.na(PIF))

#write.table(pin, 'dat/pin.txt', row.names=F)
headers <- c('PSMId', 'Label', 'ScanNr', 'ExpMass', 'CalcMass', 'RT', 'dM', 
             'Score', 'dScore', 'MassError', 'PepLen',
             'RetLen', 'dRT', 'PIF', 'NumIsoPeaks',
             'Charge1', 'Charge2', 'Charge3', 'Charge4', 'Charge5',
             'Mass', 'dM', 'absdM', 'Peptide', 'Proteins')
write(x=headers, file='dat/pin.txt', ncolumns=length(headers), sep='\t')

write_tsv(pin, 'dat/pin.txt', col_names=F, append=T, na='')

## percolator command:
# percolator dat/pin.txt --doc --post-processing-tdc --testFDR 0.05 
# --results-peptides dat/pout-peptides.txt --decoy-results-peptides dat/pout-decoy-peptides.txt 
# --results-psms dat/pout-psms.txt --decoy-results-psms dat/pout-decoy-psms.txt 
# --results-proteins dat/pout-proteins.txt --decoy-results-proteins dat/pout-decoy-proteins.txt 
# --tab-out dat/pout-features.txt

# append percolator results to original evidence file
#pout <- read_tsv('dat/percolator_181217/pout-psms.txt')
pout <- read_tsv('dat/percolator_181217_nodoc/pout-psms.txt')

ev.out <- ev %>%
  # alias the modified peptide ID as the peptide ID
  rename(`Sequence ID`=`Peptide ID`) %>%
  rename(`Peptide ID`=`Mod. peptide ID`) %>%
  # select only few fields
  select(Sequence, Proteins, `Leading razor protein`, `Raw file`, 
         `Retention time`, `Retention length`, `PIF`, `PEP`, `Intensity`,
         `Reporter intensity corrected 0`, `Reporter intensity corrected 1`, 
         `Reporter intensity corrected 2`, `Reporter intensity corrected 3`, 
         `Reporter intensity corrected 4`, `Reporter intensity corrected 5`, 
         `Reporter intensity corrected 6`, `Reporter intensity corrected 7`, 
         `Reporter intensity corrected 8`, `Reporter intensity corrected 9`,
         Reverse, `Peptide ID`, `Sequence ID`, Modifications, id)

ev.out$PEP.new <- NA
ev.out$PEP.new[pout$PSMId+1] <- pout$posterior_error_prob
ev.out$qval <- NA
ev.out$qval[pout$PSMId+1] <- pout$`q-value`

write_tsv(ev.out, path='dat/ev.perc_181217_nodoc.txt')

##

sum(!is.na(ev.out$PEP.new) & ev.out$PEP < 0.05)
sum(ev.out$PEP.new < 0.05, na.rm=T)

b <- ev.out %>% filter(!is.na(PEP.new)) %>% sample_n(1e4)
plot(b$PEP, b$PEP.new, log='xy')


b <- ev.out %>%
  filter(!is.na(PEP.new)) %>%
  filter(PEP > 0 & PEP.new > 0) %>%
  mutate_at('PEP', funs(ifelse(. > 1, 1, .))) %>%
  mutate(dPEP=log10(PEP/PEP.new)) %>%
  arrange(desc(PEP)) %>%
  mutate(bin=cut(PEP, breaks=c(0, 1e-3, 1e-2, 5e-2, 1e-1, 5e-1, 9e-1, 1)))

plot.scatter <- b %>% sample_n(1e4) %>%
ggplot(aes(x=PEP, y=PEP.new)) +
  geom_point(alpha=0.3) +
  geom_abline(slope=1, intercept=0, color='red') +
  scale_x_log10(limits=c(1e-10, 1)) + 
  scale_y_log10(limits=c(1e-10, 1))


library(ggridges)
dpep.eq <- parse(text=paste0('log[10](frac(\'Spectral PEP\',\'Updated PEP\'))'))

plot.ridges <- ggplot(b) +
  geom_density_ridges(aes(x=dPEP, y=bin, group=bin), 
                    rel_min_height=0.01) +
  geom_vline(xintercept=0, color='red', linetype='longdash') +
  #annotate(geom='text', label='Decreased Confidence', x=-1.5, y=5.5, size=5) +
  #annotate(geom='text', label='Increased Confidence', x=1.5, y=5.5, size=5) +
  scale_x_continuous(limits=c(-2, 1), expand=c(0.01, 0)) +
  scale_y_discrete(expand=(c(0.01, 0)), position='right') +
  #scale_y_continuous(breaks=seq(1,length(levels(ev.f$bin))), 
  #                   labels=levels(ev.f$bin), 
  #                   expand=c(0.01, 0), position='right',
  #                   sec.axis=sec_axis(~., labels=NULL, name='Density')) +
  labs(x=dpep.eq, y=NULL) +
  #labs(x='dPEP', y=NULL) +
  #theme_bert()
  theme_ridges() +
  theme(axis.title.x = element_text(hjust=0.5),
        axis.title.y = element_text(hjust=0.5))

pdf(file='percolator_181217_nodoc.pdf', width=7, height=5)
grid.arrange(ggplotGrob(plot.scatter), ggplotGrob(plot.ridges), nrow=1, ncol=2)
dev.off()

