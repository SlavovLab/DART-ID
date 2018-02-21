## Figure 1 -------

library(tidyverse)
library(pracma)
library(gridExtra)
library(grid)
library(gtable)
library(RColorBrewer)
library(ggridges)
library(reshape2)
source('lib.R')

ev <- read_tsv('dat/ev.adj.Fit3c.RTL.txt')
ev <- read_tsv('dat/ev.adj.Fit3.txt')
ev <- read_tsv('dat/ev.adj.Fit2.txt')

## Update Demo -------

x <- seq(0, 250, by=0.1)
y1 <- dlnorm(x, meanlog=4, sdlog=0.5)
y2 <- dnorm(x, mean=120, sd=2)

#plot(x, y2, type='l', col='blue')
#lines(x, y1, type='l', col='red')

density.plot <- 
ggplot(data.frame()) +
  geom_path(aes(x=x, y=y1, color='Null RT Density (All RTs)')) +
  geom_path(aes(x=x, y=y2, color='Predicted RT Density\nfor Top Match')) +
  #geom_vline(xintercept=80, color='black', linetype='longdash') +
  scale_color_manual(values=c('red', 'blue', 'black')) +
  scale_x_continuous(limits=c(0, 170), expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0.003)) +
  labs(x=NULL, y='Density', color='Distribution') +
  theme_bert() +
  theme(
    #plot.margin = unit(c(0.5,0.5,0.5,0.5), 'cm'),
    legend.position = c(0.1, 0.5),
    legend.background = element_rect(fill='white', color=NULL, size=0.25),
    legend.key.height = unit(1.25, 'cm'),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

den.plot.1 <- 
density.plot + 
  #geom_vline(xintercept=80, color='black', linetype='longdash') +
  geom_segment(aes(x=80, y=0, xend=80, yend=0.2, color='Observed RT'), linetype='longdash') +
  #annotate('text', x=80, y=-0.015, label='Observed RT', size=5) +
  scale_color_manual(values=c('red', 'black', 'blue'), guide=F)

den.plot.2 <- 
density.plot + 
  #geom_vline(xintercept=119, color='black', linetype='longdash') +
  geom_segment(aes(x=119, y=0, xend=119, yend=0.2, color='Observed RT'), linetype='longdash') +
  #annotate('text', x=119, y=-0.015, label='Observed RT', size=5) +
  scale_color_manual(values=c('red', 'black', 'blue')) +
  labs(y=' ') 
den.plot.1 <- ggplotGrob(den.plot.1)
den.plot.2 <- ggplotGrob(den.plot.2)


bdf <- data.frame(
  a=factor(c('Spectra Only', 'Updated'), levels=c('Spectra Only', 'Updated')),
  b=as.numeric(c(1.5, 3))
)

bar.plot.1 <- ggplot(bdf) +
#ggplot(bdf) +
  geom_bar(aes(x=a, y=b), stat='identity',
           fill='grey90', color='black') +
  labs(x=NULL, y='Error Prob. (PEP)\nfor Top Match') +
  scale_y_continuous(expand=c(0,0), limits=c(0,3.2)) +
  theme_bert() %+replace% theme(
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank()
  )
bdf$b <- c(1.5, 0.5)
bar.plot.2 <- 
ggplot(bdf) +
  geom_bar(aes(x=a, y=b), stat='identity',
           fill='grey90', color='black') +
  labs(x=NULL, y=' ') +
  scale_y_continuous(expand=c(0,0), limits=c(0,3.2)) +
  theme_bert() %+replace% theme(
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank()
  )

bar.plot.1 <- ggplotGrob(bar.plot.1)
bar.plot.2 <- ggplotGrob(bar.plot.2)

grid.arrange(den.plot.1, den.plot.2, bar.plot.1, bar.plot.2, 
             ncol=2, nrow=2, widths=c(1,1), heights=c(3,2))



## dPEP Ridges ------

ev.f <- ev %>%
  filter(!is.na(PEP.new)) %>%
  select(c('Raw file', 'Sequence', 'PEP', 'PEP.new', 'Retention time', 'muijs')) %>%
  mutate(PEP=ifelse(PEP > 1, 1, PEP), 
         PEP=ifelse(PEP == 0, .Machine$double.xmin, PEP)) %>%
  mutate(dPEP=log2(PEP/PEP.new),
         dRT=log2(abs(`Retention time`-muijs))) %>%
  mutate_at(vars(dPEP, dRT), funs(ifelse(is.infinite(.), NA, .))) %>%
  #mutate(bin=cut(PEP, breaks=c(0, 5e-2, 1e-1, 5e-1, 9e-1, 1)))
  mutate(bin=cut(PEP, breaks=seq(0, 1, by=0.1)))

## dPEP Ridges - Plotting ----

ggplot(ev.f) +
  #geom_density(aes(dPEP))
  geom_density_ridges(aes(x=dRT, y=bin, group=bin), rel_min_height=0.01) +
  geom_vline(xintercept=0) +
  scale_x_continuous(limits=c(-3, 3))
  #theme_ridges()

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

