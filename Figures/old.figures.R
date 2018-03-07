## Update Demo - BAD ------

x <- seq(0, 250, by=0.1)
y1 <- dlnorm(x, meanlog=4, sdlog=0.5)
y2 <- dnorm(x, mean=120, sd=2)

# color palette
pal <- brewer.pal(4, 'Set1')

#plot(x, y2, type='l', col='blue')
#lines(x, y1, type='l', col='red')

density.plot <- 
  ggplot(data.frame()) +
  geom_path(aes(x=x, y=y2, color='b'), size=1) +
  geom_path(aes(x=x, y=y1, color='a'), size=1) +
  geom_segment(aes(x=90, y=0, xend=90, yend=0.25, color='c'), linetype='dashed', size=0.75) +
  geom_segment(aes(x=119, y=0, xend=119, yend=0.25, color='c'), linetype='dashed', size=0.75) +
  #annotate(geom='text', x=90, y=-0.05, label='Observation 1', size=5) +
  scale_color_manual(values=c(pal[1], pal[2], 'black'),
                     labels=c('Null RT Density (All RTs)', 
                              'Predicted RT Density\nfor Top Match', 
                              'Observed RT')) +
  scale_x_continuous(limits=c(0, 160), expand=c(0,0)) +
  scale_y_continuous(limits=c(0, 0.25), expand=c(0.01,0)) +
  labs(x=NULL, y='Density', color='Distribution') +
  theme_bert() +
  theme(
    plot.margin = unit(c(0.25,0.25,0,0.25), 'cm'),
    #legend.position = c(0.1, 0.5),
    legend.position = c(0.23, 0.55),
    #legend.background = element_rect(fill='white', color=NULL, size=0.25),
    legend.background = element_blank(),
    legend.key.height = unit(1.1, 'cm'),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    #axis.line = element_line(size=0.25, color='black')
  )

den.plot <- ggplotGrob(density.plot)


bdf <- data.frame(
  a=factor(c('Spectra Only', 'Updated'), levels=c('Spectra Only', 'Updated')),
  b=as.numeric(c(1.9, 3))
)

bar.plot.1 <- ggplot(bdf) +
  #ggplot(bdf) +
  geom_bar(aes(x=a, y=b), stat='identity',
           fill='grey90', color='black') +
  labs(x=NULL, y='Error Prob. (PEP)\nfor Top Match') +
  scale_y_continuous(expand=c(0,0), limits=c(0,3.5)) +
  #annotate(geom='text', size=5, x=0.9, y=2.5, label='PSM B') +
  labs(title='PSM A') +
  theme_bert() %+replace% theme(
    plot.margin = unit(c(0.5,0.1,0.4,0.25), 'cm'),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank(),
    plot.title=element_text(size=16)
  )
bdf$b <- c(1.1, 0.5)
bar.plot.2 <- 
  ggplot(bdf) +
  geom_bar(aes(x=a, y=b), stat='identity',
           fill='grey90', color='black') +
  labs(x=NULL, y=' ') +
  scale_y_continuous(expand=c(0,0), limits=c(0,3.2)) +
  labs(title='PSM B') +
  theme_bert() %+replace% theme(
    plot.margin = unit(c(0.5,0.25,0.4,0.1), 'cm'),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank(),
    plot.title=element_text(size=16)
  )

bar.plot.1 <- ggplotGrob(bar.plot.1)
bar.plot.2 <- ggplotGrob(bar.plot.2)

#grid.arrange(den.plot.1, den.plot.2, bar.plot.1, bar.plot.2, 
#             ncol=2, nrow=2, widths=c(1,1), heights=c(3,2))

lay <- rbind(c(1,1,1,1,1,1),
             c(2,2,2,3,4,NA),
             c(5,5,5,6,6,6))
gs <- list(den.plot, 
           textGrob('Retention Time (mins)', hjust=0.6), 
           #textGrob('PSM A'), textGrob('PSM B'), 
           textGrob(''), textGrob(''), 
           bar.plot.1, bar.plot.2)
grid.arrange(grobs=gs, layout_matrix=lay, 
             widths=c(1,1,1,1,1,1), heights=c(1,0.1,1))


## dPEP Ridges ------

ev.f <- ev %>%
  filter(!is.na(PEP.new)) %>%
  select(c('Raw file', 'Sequence', 'PEP', 'PEP.new', 'Retention time', 'muijs')) %>%
  mutate(PEP=ifelse(PEP > 1, 1, PEP), 
         PEP=ifelse(PEP == 0, .Machine$double.xmin, PEP)) %>%
  mutate(dPEP=log2(PEP/PEP.new),
         dRT=log10(abs(`Retention time` - muijs))) %>%
  mutate_at(vars(dPEP, dRT), funs(ifelse(is.infinite(.), NA, .))) %>%
  #mutate(bin=cut(PEP, breaks=c(0, 5e-2, 1e-1, 5e-1, 9e-1, 1)))
  mutate(bin=cut(PEP, breaks=seq(0, 1, by=0.1)))
#mutate(bin=cut(PEP, breaks=seq(0, 1, by=0.2)))

## dPEP/dRT Ridges - Plotting ----

ggplot(ev.f) +
  #geom_density(aes(dPEP))
  geom_density_ridges(aes(x=dRT, y=bin, group=bin), rel_min_height=0.01, scale=1.2) +
  #stat='binline', bins=60) +
  geom_vline(xintercept=0, color='red') +
  scale_x_continuous(limits=c(-2.5, 2.5), breaks=seq(-2, 2, by=1)) +
  scale_y_discrete(expand=c(0.05,0)) +
  labs(x='dRT (min)', y='Spectral PEP') +
  theme_ridges() +
  theme(
    axis.title.x = element_text(size=16, hjust=0.5, vjust=0.5),
    axis.title.y = element_text(size=16, hjust=0.5, vjust=0.5)
  )

## -----

ev.f <- ev %>%
  select(c('Sequence', 'PEP', 'PEP.new')) %>%
  filter(!is.na(PEP.new)) %>%
  filter(PEP > 0 & PEP.new > 0) %>%
  mutate_at('PEP', funs(ifelse(. > 1, 1, .))) %>%
  mutate(dPEP=log10(PEP/PEP.new)) %>%
  arrange(desc(PEP)) %>%
  mutate(bin=cut(PEP, breaks=c(0, 5e-2, 1e-1, 5e-1, 9e-1, 1)))


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



## ------

library(ggridges)
library(rmutil)


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

## ------

pep.bins <-
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
