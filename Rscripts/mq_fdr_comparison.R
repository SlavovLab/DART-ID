library(tidyverse)
library(pracma)
library(ggridges)
source('Rscripts/lib.R')

ev <- read_tsv("/gd/bayesian_RT/Alignments/SQC_20180803_5exp_parametric_mixture_v2/ev_updated.txt")

ev.f <- ev %>%
  mutate_at(c('PEP', 'pep_updated'), funs(ifelse(. > 1, 1, .)))


# compare # of false positives --------------------------------------------

k <- 500
peps <- logseq(1e-8, 1, k)

fp_mq <- vector(length=k)
fp_dart <- vector(length=k)

for(i in 1:length(peps)) {
  fp_mq[i] <- ceil(sum(ev.f$PEP[ev.f$PEP < peps[i]]))
  fp_dart[i] <- ceil(sum(ev.f$PEP[ev.f$pep_updated < peps[i]]))
}

# plot --------------------------------------------------------------------

plot(fp_mq, fp_dart, pch=16)
abline(a=0, b=1, col='red', lwd=2)


# with pep on x-axis ------------------------------------------------------

plot(log10(peps), log2(fp_mq), type='l', lwd=3, col='red',
     xlab='log10 PEP Thresold', ylab='log2 false positives',
     main='Number of false positives')
lines(log10(peps), log2(fp_dart), lwd=3, col='blue')

legend('topleft', c('MaxQuant', 'DART-ID'), col=c('red', 'blue'), lwd=4,
       cex=1)


# compare # of false positives in pep bins --------------------------------

ev.f$pep_bin <- cut(ev.f$PEP, breaks=seq(0, 1, by=0.05))
ev.f$dpep <- log2(ev.f$PEP / ev.f$pep_updated)

ev.f %>% 
  filter(dpep > -1 & dpep < 6) %>% 
  filter(!is.na(pep_new)) %>%
  filter(PEP < 0.95) %>%
ggplot() +
  geom_density_ridges(aes(x=dpep, y=pep_bin, group=pep_bin), stat='binline', bins=50) +
  geom_vline(xintercept=0, color='red') +
  labs(x='Change in PEP (positive = upgrade)') +
  theme_ridges()
