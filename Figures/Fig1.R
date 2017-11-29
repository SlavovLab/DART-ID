## Figure 1 -------

library(tidyverse)
library(pracma)
library(gridExtra)
library(grid)
library(gtable)
library(RColorBrewer)
source('lib.R')

ev <- read_tsv('dat/ev.adj.Fit2.expcentric.txt')
#ev <- parse.ev.adj(ev)

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

## Fig 1A - PEP vs. PEP.new Scatter/Density -----

pdf(file='manuscript/Figs/Fig_1A.pdf', width=7, height=7, onefile=F)

p <- ggplot(ev.f, aes(x=PEP, y=PEP.new)) +
#ggplot(ev.f, aes(x=PEP, y=PEP.new)) +
  stat_bin2d(bins=100, drop=TRUE, geom='tile', aes(fill=..count..)) +
  #stat_density2d(aes(alpha=..level.., fill=..level..), bins=9, geom='polygon') +
  geom_abline(slope=1, intercept=0, color='red') +
  scale_fill_gradientn(colors=bHeatmap, 
                       values=c(0, 0.05, 0.1, 0.2, 0.5, 1),
                       labels=c(0, 500, 1000, 1500, 2000)) +
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

## -----

ev.f <- ev %>%
  select(c('Sequence', 'PEP', 'PEP.new')) %>%
  filter(!is.na(PEP.new)) %>%
  filter(PEP > 0 & PEP.new > 0) %>%
  mutate_at('PEP', funs(ifelse(. > 1, 1, .))) %>%
  mutate(dPEP=log10(PEP/PEP.new)) %>%
  arrange(desc(PEP)) %>%
  mutate(bin=cut(PEP, breaks=c(0, 5e-2, 1e-1, 5e-1, 9e-1, 1)))

## binned violin plot ----

p <- ggplot(ev.f) +
  geom_hline(yintercept=0, color='red', linetype='longdash') +
  geom_violin(aes(x=bin, y=dPEP, group=bin), fill='black', scale='width', trim=F)  +
  #scale_x_discrete(expand=c(0, 1)) +
  labs(x='Spectral PEP', y=parse(text=paste0('log[10](frac(\'Spectral PEP\',\'Updated PEP\'))'))) +
  #labs(x='PEP', y='eeeeeeeeeeeeeeeeee') +
  theme_bert() + 
  theme(axis.title.y=element_text(vjust=0.5, angle=0),
        axis.text.x=element_text(angle=45, hjust=1, margin=margin(t=4)),
        panel.grid.major.y=element_line(color='ivory3', size=0.25),
        panel.grid.minor.y=element_line(color='ivory2', size=0.25))

main.grob <- ggplotGrob(p)
inc.text <- textGrob(label='Increased Confidence', y=0.85)
dec.text <- textGrob(label='Decreased Confidence', y=0.15)

#rm(g)
g <- gtable_add_grob(main.grob, t=6, l=2, grobs=list(inc.text, dec.text), name=c(0, 1))

grid.newpage()
grid.draw(g)


## PEP Bin -- with ggridges ---------
library(ggridges)

dpep.eq <- parse(text=paste0('log[10](frac(\'Spectral PEP\',\'Updated PEP\'))'))

#pdf(file='manuscript/Figs/Fig_1B.pdf', width=7, height=7, onefile=F)

p <- ggplot(ev.f) +
  geom_density_ridges(aes(x=dPEP, y=bin, group=bin), 
                      rel_min_height=0.00) +
  geom_vline(xintercept=0, color='red', linetype='longdash') +
  #annotate(geom='text', label='Decreased Confidence', x=-1.5, y=5.5, size=5) +
  #annotate(geom='text', label='Increased Confidence', x=1.5, y=5.5, size=5) +
  scale_x_continuous(expand=c(0.01, 0)) +
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

main.grob <- ggplotGrob(p)

g <- gtable_add_grob(main.grob, t=6, l=5, name='yaxis.title.right', z=Inf,
                     textGrob(label='Spectral PEP Bin', hjust=0.5, y=0.7, gp=gpar(fontsize=14)))
g <- gtable_add_grob(g, t=6, l=3, name='yaxis.title.left', z=Inf,
                     textGrob(label='Density', rot=90, gp=gpar(fontsize=14), hjust=0.5, vjust=0.5))
g$widths[3] <- unit(0.4, 'in')
g$widths[5] <- unit(1.7, 'in')

grid.newpage()
grid.draw(g)

#dev.off()

## PEP Fold Change Line Plot -------

#t <- c(0, logseq(1e-10, 1, n=50))
t <- seq(0,1,by=0.05)
p <- c()
for(i in 1:length(t)-1) {
  # number of PSMs with PEP.new within interval / same with PEP
  p[i] <- sum(ev.f$PEP.new > t[i] & ev.f$PEP.new < t[i+1]) / 
    sum(ev.f$PEP > t[i] & ev.f$PEP < t[i+1])
  #p[i] <-sum(ev.f$PEP.new < t[i+1]) / sum(ev.f$PEP < t[i+1])
}

ggplot(data.frame(), aes(x=t[-1], y=p)) +
  geom_hline(yintercept=1, color='red', linetype='longdash') +
  geom_path() +
  geom_point() +
  #scale_x_log10(name='PEP Interval', breaks=logseq(1e-10, 1, 6)) +
  #              #labels=fancy_scientific) +
  scale_x_continuous(breaks=seq(0, 1, by=0.2)) +
  labs(x='PEP Threshold', 
       y=parse(text='frac(\'# Updated PEP < Threshold\', \'# Spectral PEP < Threshold\')')) +
  theme_bert()




## Alignment Ridges -----

ev.f <- ev %>%
  select(c('Sequence', 'Leading razor protein', 'Raw file', 'PEP', 'PEP.new', 'Peptide ID', 
           'Retention time', 'muijs', 'sigmas')) %>%
  filter(!is.na(PEP.new)) %>%
  #filter(PEP < 0.05) %>%
  filter(!grepl('REV*', `Leading razor protein`)) %>% # Remove Reverse matches
  filter(!grepl('CON*',`Leading razor protein`)) # Remove Contaminants

# clean up raw file names and extract an experiment ID from it
# need the experiment ID (19A, 30B, etc.), to match with sample metadata from the
# experiment description excel sheet/.csv
ev.f$file <- clean.file.name(ev.f$`Raw file`)
ev.f$exp <- str_extract(ev.f$file, '[1-9][0-9][A-Z](\\d)?')

psm.count <- ev.f %>%
  group_by(`Peptide ID`) %>%
  summarise(Sequence=unique(Sequence), 
            Protein=unique(`Leading razor protein`), 
            count=length(`Peptide ID`)) %>%
  arrange(desc(count))


library(ggridges)
library(rmutil)
## ------

i <- 5
pep.id <- psm.count$`Peptide ID`[i]
ev.p <- ev.f %>% 
  filter(`Peptide ID`==pep.id)

# take the top [num.files]
num.files <- 5

alignment.dat <- ev.p %>%
  group_by(exp) %>%
  summarise(count=length(exp), 
            mu=first(unique(muijs)), 
            sigma=first(unique(sigmas))) %>%
  arrange(desc(count)) %>%
  top_n(num.files, wt=count)

files <- alignment.dat$exp[1:num.files]

density.dat <- data.frame(x=numeric(), y=numeric(), exp=character())
range <- seq(0,250,by=0.5)
for(j in 1:num.files) {
  density.dat <- rbind(density.dat, data.frame(
    x=as.numeric(range),
    y=as.numeric(dlaplace(range, m=alignment.dat$mu[j], s=alignment.dat$sigma[j])),
    exp=as.character(files[j])
  ))
}

ev.p %>% 
  filter(exp %in% files) %>%
  filter(PEP.new < 0.05) %>%
  ggplot() +
  geom_density_ridges(aes(x=`Retention time`, y=exp, group=exp), 
                      rel_min_height=0.01) +
  geom_path()
  #geom_hline(yintercept=2, color='red') +
  #scale_x_continuous(limits=c(40, 150)) +
  labs(title=unique(ev.p$Sequence))

