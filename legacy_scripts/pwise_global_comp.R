library(tidyverse)
library(RColorBrewer)
source('lib.R')

ev.pwise.in <- read_tsv('dat/evidence+dRT.elite.txt')
#ev.pwise <- read_tsv('dat/ev+dRT.elite.txt')
ev.pwise <- ev.pwise.in %>%
  mutate(dRT=RT.corrected-RT.lib) %>%
  select(c('Raw.file', 'Sequence', 'PEP', 'PEP.new', 'Retention.time', 'dRT', 'Best.MS.MS')) %>%
  rename(`Retention time`='Retention.time',
         `Raw file`='Raw.file',
         `Best MS/MS`='Best.MS.MS') %>%
  mutate(pep_exp=paste(Sequence, `Raw file`, sep='|'),
         method='Pairwise')
  
ev.global.in <- read_tsv('dat/ev.adj.Fit3c.RTL.txt')
#ev.global.in <- read_tsv('dat/ev.adj.Fit3.txt')
ev.global <- ev.global.in %>%
  filter(!is.na(muijs)) %>%
  mutate(dRT=`Retention time`-muijs) %>%
  select(c('Raw file', 'Sequence', 'PEP', 'PEP.new', 'Retention time', 'dRT', 'Best MS/MS')) %>%
  mutate(pep_exp=paste(Sequence, `Raw file`, sep='|'),
         method='Global')

#pep_exp_pairs <- intersect(ev.pwise$pep_exp, ev.global$pep_exp)
#ev.global <- ev.global[ev.global$pep_exp %in% pep_exp_pairs,]
#ev.pwise <- ev.pwise[ev.pwise$pep_exp %in% pep_exp_pairs,]

msms_pairs <- intersect(ev.pwise$`Best MS/MS`, ev.global$`Best MS/MS`)
ev.global <- ev.global[ev.global$`Best MS/MS` %in% msms_pairs,]
ev.pwise <- ev.pwise[ev.pwise$`Best MS/MS` %in% msms_pairs,]

ev.global <- ev.global[order(ev.global$`Best MS/MS`),]
ev.pwise <- ev.pwise[order(ev.pwise$`Best MS/MS`),]

ev.tot <- rbind(ev.global, ev.pwise) %>%
  filter(!is.na(pep_bin)) %>%
  mutate(pep_bin=cut(PEP, breaks=c(0, 1e-3, 1e-2, 5e-2, 1e-1, 5e-1, 1)),
         dRT=abs(dRT))

ggplot(ev.tot, aes(dRT, color=method)) +
  geom_path(stat='density') +
  facet_grid(pep_bin~.) +
  scale_x_continuous(limits=c(0, 5)) +
  labs(x='dRT (mins)', color='Method')

## ----------

df.comp <- data.frame(
  Sequence=as.character(ev.global$Sequence),
  Raw.file=as.character(ev.global$`Raw file`))
df.comp$dRT.global <- ev.global$dRT
df.comp$dRT.pwise <- ev.pwise$dRT

# make sure 0s are not 0
df.comp$dRT.pwise[df.comp$dRT.pwise==0] <- .Machine$double.xmin
df.comp$dRT.global[df.comp$dRT.global==0] <- .Machine$double.xmin

df.comp %>%
  #sample_n(1e4) %>%
  filter(Raw.file %in% unique(ev.global$`Raw file`)[1:20]) %>%
ggplot(aes(x=dRT.pwise, y=dRT.global)) +
  #geom_point(alpha=0.1) +
  #stat_density2d(geom='polygon', n=10, aes(fill=..level..)) +
  stat_bin2d(bins=16, drop=TRUE, geom='tile', aes(fill=..count..)) +
  geom_abline(intercept=0, slope=1, color='red') +
  scale_x_log10(limits=c(1e-1, 20), expand=c(0,0)) +
  scale_y_log10(limits=c(1e-1, 20), expand=c(0,0)) +
  #scale_fill_gradient(low='red', high='yellow') +
  scale_fill_gradientn(colors=brewer.pal(5, 'Blues')) +
  labs(x='dRT - Pairwise', y='dRT - Global') +
  theme_bert()



