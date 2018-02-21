## Figure 2 -------

library(tidyverse)
library(pracma)
library(gridExtra)
library(grid)
library(gtable)
library(RColorBrewer)
source('lib.R')

ev <- read_tsv('dat/ev.adj.Fit3.txt')

## PEP vs. PEP.new scatterplot -----

ev.f <- ev %>%
  select(c('Sequence', 'PEP', 'PEP.new')) %>%
  filter(!is.na(PEP.new)) %>%
  filter(PEP > 0 & PEP.new > 0) %>%
  mutate_at('PEP', funs(ifelse(. > 1, 1, .))) %>%
  mutate_at('PEP.new', funs(ifelse(. > 1, 1, .))) %>%
  mutate(dPEP=PEP-PEP.new) %>%
  arrange(desc(PEP)) %>%
  mutate(bin=cut(PEP, breaks=c(0, 1e-3, 1e-2, 5e-2, 1e-1, 5e-1, 7.5e-1, 1)))

ev.f %>%
  sample_n(1e5) %>%
  ggplot(aes(x=PEP, y=PEP.new)) +
  geom_point(alpha=0.1) +
  scale_x_log10(limits=c(1e-10, 1)) +
  scale_y_log10(limits=c(1e-10, 1))


## Fig 1A - PEP vs. PEP.new Scatter/Density -----

pdf(file='manuscript/Figs/Fig_1A.pdf', width=7, height=7, onefile=F)

p <- ggplot(ev.f, aes(x=PEP, y=PEP.new)) +
  #ggplot(ev.f, aes(x=PEP, y=PEP.new)) +
  stat_bin2d(bins=100, drop=TRUE, geom='tile', aes(fill=..count..)) +
  #stat_density2d(aes(alpha=..level.., fill=..level..), bins=9, geom='polygon') +
  geom_abline(slope=1, intercept=0, color='red') +
  scale_fill_gradientn(colors=bHeatmap, 
                       values=c(0, 0.05, 0.1, 0.2, 0.5, 1)
                       #labels=c(0, 500, 1000, 1500, 2000)) +
  )+
  scale_alpha_continuous(guide=FALSE) +
  scale_x_log10(limits=c(1e-10, 1), 
                breaks=logseq(1e-10, 1, 6), labels=fancy_scientific) +
  #breaks=logseq(1e-10, 1, 6)) +
  scale_y_log10(limits=c(1e-10, 1), 
                breaks=logseq(1e-10, 1, 6), labels=fancy_scientific) +
  #breaks=logseq(1e-10, 1, 6)) +
  theme_bert() +
  theme(legend.position=c(0.9, 0.3),
        legend.key.height = unit(0.05, 'npc'),
        legend.text = element_text(colour = "black", size = 12),
        plot.margin=unit(c(0,0,0.5,0.5), 'cm')) +
  labs(fill='Count', x='Spectral PEP', y='Updated PEP')

density.top <- ggplot(ev.f, aes(PEP)) + 
  #ggplot(ev.f, aes(PEP)) + 
  geom_histogram(bins=30, aes(y=..density..), color='black', fill='white') +
  stat_density(adjust=8, geom='path', color='red') + 
  scale_x_log10(limits=c(1e-10, 1)) +
  labs(y='Density') +
  theme_bert() +
  theme(axis.text=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x=element_line(color='ivory3', size=0.5),
        panel.grid.minor.x=element_line(color='ivory2', size=0.5),
        plot.margin=unit(c(0.5,0.5,0,0.5), 'cm'))
density.right <- ggplot(ev.f, aes(PEP.new)) + 
  geom_histogram(bins=30, aes(y=..density..), color='black', fill='white') +
  stat_density(adjust=5, geom='path', color='red') + 
  scale_x_log10(limits=c(1e-10, 1)) + 
  labs(y='Density') +
  coord_flip() +
  theme_bert() +
  theme(axis.text=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y=element_line(color='ivory3', size=0.5),
        panel.grid.minor.y=element_line(color='ivory2', size=0.5),
        plot.margin=unit(c(0.5,0.5,0.5,0), 'cm'))

bTitle <- textGrob(label='[Method Name] PEP Shift', 
                   x=unit(1, 'cm'),
                   just=c('left', 'centre'),
                   gp=gpar(fontsize=24))

main.grob <- ggplotGrob(p)
den.top.grob <- ggplotGrob(density.top)
den.right.grob <- ggplotGrob(density.right)

den.top.grob$widths <- main.grob$widths
den.right.grob$heights <- main.grob$heights

grid.arrange(#bTitle,        nullGrob(),
  den.top.grob,  nullGrob(),
  main.grob,     den.right.grob, 
  #nrow=3, ncol=2, widths=c(5, 1), heights=c(1, 2, 7))
  nrow=2, ncol=2, widths=c(5, 1), heights=c(1, 5))

dev.off()

## PEP vs. PEP.new fractions ----

ev.f <- ev %>%
  filter(!is.na(PEP.new)) %>%
  select(c('Raw file', 'Sequence', 'PEP', 'PEP.new', 'Retention time', 'muijs')) %>%
  mutate(PEP=ifelse(PEP > 1, 1, PEP), 
         PEP=ifelse(PEP == 0, .Machine$double.xmin, PEP)) %>%
  mutate(dPEP=log2(PEP/PEP.new),
         dRT=log10(abs(`Retention time`-muijs))) %>%
  mutate_at(vars(dPEP, dRT), funs(ifelse(is.infinite(.), NA, .))) %>%
  #mutate(bin=cut(PEP, breaks=c(0, 5e-2, 1e-1, 5e-1, 9e-1, 1)))
  mutate(bin=cut(PEP, breaks=seq(0, 1, by=0.1)))


dmat <- matrix(nrow=10, ncol=10)
for(i in 1:10) {
  b <- levels(ev.f$bin)[i]
  ev.a <- ev.f %>% filter(bin == b)
  bin.new <- cut(ev.a$PEP.new, breaks=seq(0, 1, by=0.1))
  for(j in 1:10) {
    b2 <- levels(bin.new)[j]
    dmat[i, j] <- nrow(ev.a %>% filter(bin.new==b2)) / nrow(ev.a)
  }
}

dmat <- melt(dmat, varnames=c('y', 'x'))

labels <- rep('=100', 10)
labels <- paste(labels, ' (n=', ev.f %>% group_by(bin) %>% summarise(n=length(bin)) %>% pull(n), ')',
                sep='')

ggplot(dmat, aes(x=x, y=y, fill=value*100)) +
  geom_tile(color='black') +
  geom_text(aes(label=round(value*100, digits=1))) +
  scale_fill_gradient(low='white', high='red') +
  scale_x_continuous(expand=c(0,0), breaks=seq(1,10), labels=levels(ev.f$bin)) +
  scale_y_continuous(expand=c(0,0), breaks=seq(1,10), labels=levels(ev.f$bin),
                     sec.axis=sec_axis(trans=~., breaks=seq(1,10), labels=labels)) +
  labs(x='PEP.new', y='PEP', fill='Percent') +
  theme_bert() +
  theme(
    axis.text.x=element_text(size=10, angle=45, hjust=1, vjust=1),
    axis.text.y=element_text(size=10)
  )


## PEP Fold Change Line Plot -------

source('lib.R')

ev <- read_tsv('dat/ev.adj.Fit3c.txt')
ev.perc <- read_tsv('dat/ev.perc_181217_noPEP.txt')
ev.perc.nodoc <- read_tsv('dat/ev.perc_181217_nodoc.txt')

## -----

exps <- list(RTLib=ev,
             Percolator=ev.perc)
#Percolator_NoDoc=ev.perc.nodoc)

df <- fold.change.comp(exps, num.steps=100)
df <- fold.change.comp(exps, begin=0, end=1, num.steps=30, log=F)

ggplot(df, aes(x=x, y=PEP, color=Method)) +
  geom_path() +
  scale_x_log10() +
  annotation_logticks(sides='b') +
  #scale_y_continuous(limits=c(0.95,2.25), breaks=c(1, 1.25, 1.5, 1.75, 2, 2.25)) +
  labs(x='PEP Threshold', y='Fold Change Increase in IDs',
       title=paste0('Fold Change Increase of PSM IDs\n',
                    '= #Adjusted PEPs / #Original PEPs above PEP Threshold'))

pdf('manuscript/Figs/Fig_1B.pdf')

ggplot(df, aes(x=x, y=PEP, color=Method)) +
  #geom_hline(yintercept=1, color='red', linetype='longdash') +
  geom_path() +
  #geom_point() +
  scale_x_log10(name='PEP Interval', breaks=logseq(1e-5, 1, 6)) +
  #labels=fancy_scientific) +
  annotation_logticks(sides='b') +
  #scale_x_continuous(breaks=seq(0, 1, by=0.2)) +
  scale_y_continuous(breaks=seq(0.4,2,by=0.2)) +
  labs(x='PEP Threshold', 
       y=parse(text='frac(\'# Updated PEP < Threshold\', \'# Spectral PEP < Threshold\')'),
       #y='asdf',
       #title=parse(text="frac(\'#Adjusted PEPs\', \'#Original PEPs above PEP Threshold\')")) +
       title=paste0("Fold Change Increase of Confident PSMs")) +
  theme_bert()

dev.off()
