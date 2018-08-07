library(tidyverse)
library(reshape2)

## PEP vs. PEP.new fractions ----

ev.f <- ev %>%
  filter(!is.na(pep_new)) %>%
  select(c('Raw file', 'Sequence', 'PEP', 'pep_new', 'Retention time', 'muij')) %>%
  mutate(PEP=ifelse(PEP > 1, 1, PEP), 
         PEP=ifelse(PEP == 0, .Machine$double.xmin, PEP)) %>%
  mutate(dPEP=log2(PEP/pep_new),
         dRT=log10(abs(`Retention time`-muij))) %>%
  mutate_at(vars(dPEP, dRT), funs(ifelse(is.infinite(.), NA, .))) %>%
  #mutate(bin=cut(PEP, breaks=c(0, 5e-2, 1e-1, 5e-1, 9e-1, 1)))
  mutate(bin=cut(PEP, breaks=seq(0, 1, by=0.1)))


dmat <- matrix(nrow=10, ncol=10)
for(i in 1:10) {
  b <- levels(ev.f$bin)[i]
  ev.a <- ev.f %>% filter(bin == b)
  bin.new <- cut(ev.a$pep_new, breaks=seq(0, 1, by=0.1))
  for(j in 1:10) {
    b2 <- levels(bin.new)[j]
    dmat[i, j] <- nrow(ev.a %>% filter(bin.new==b2)) / nrow(ev.a)
  }
}

dmat <- melt(dmat, varnames=c('y', 'x'))

labels <- rep('=100', 10)
labels <- paste(labels, ' (n=', ev.f %>% group_by(bin) %>% summarise(n=length(bin)) %>% pull(n), ')',
                sep='')

## ------

#pep.bins <-
ggplot(dmat, aes(x=x, y=y, fill=value*100)) +
  geom_tile(color='black') +
  geom_text(aes(label=round(value*100, digits=1)), size=4) +
  scale_fill_gradientn(colors=c('white','red'), values=c(0, 1.25), guide=F) +
  scale_x_continuous(expand=c(0,0), breaks=seq(1,10), labels=levels(ev.f$bin)) +
  scale_y_continuous(expand=c(0,0), breaks=seq(1,10), labels=levels(ev.f$bin),
                     sec.axis=sec_axis(trans=~., breaks=seq(1,10), labels=labels)) +
  labs(x='Updated PEP', y='Spectral PEP', fill='Percent') +
  theme_bert() +
  theme(
    axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1),
    axis.text.y=element_text(size=12),
    panel.border = element_rect(size=0.25, color='black', fill=NA)
  )

ridges <- 
  ggplot(ev.f) +
  #geom_density(aes(dPEP))
  geom_density_ridges(aes(x=dRT, y=bin, group=bin), 
                      rel_min_height=0.01, scale=1.5) +
  #stat='binline', bins=60) +
  geom_vline(xintercept=0, color='red') +
  scale_x_continuous(limits=c(-2.5, 2.5), breaks=seq(-2, 2, by=1)) +
  scale_y_discrete(expand=c(0.05,0)) +
  #scale_fill_gradientn(colors=rev(viridis(n=10)), guide=F) +
  labs(x='log2(dRT) (min)', y='Spectral PEP') +
  theme_bert() +
  #theme_ridges() +
  theme(
    plot.margin = unit(c(0.25, 0.25, 1.25, 0.25), 'cm'),
    axis.title.x = element_text(size=16, hjust=0.5, vjust=0.5),
    axis.title.y = element_text(size=16, hjust=0.5, vjust=0.5),
    panel.grid = element_line(color='grey50', size=0.25)
  )

pep.bins <- ggplotGrob(pep.bins)
ridges <- ggplotGrob(ridges)

grid.arrange(pep.bins, ridges, nrow=1, ncol=2, widths=c(2,1))
