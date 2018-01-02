# load different fits

#load('dat/params.Fit2.RData')
#pars2 <- pars
load('dat/params.Fit3.RData')
pars3 <- pars
load('dat/params.Fit3b.RData')
pars3b <- pars
load('dat/params.Fit3c.RData')
pars3c <- pars

evidence <- read_tsv('dat/evidence_elite.txt')

load('dat/alignment_data.RData')

get.mean.res <- function(exp) {
  ## get experiment indices
  exp_indices <- data$muij_to_exp==exp
  observed <- data$retention_times[exp_indices[data$muij_map]]
  predicted <- (muijs[data$muij_map])[exp_indices[data$muij_map]]
  
  residual <- observed - predicted
  
  mean(abs(residual))
}


muijs <- pars3[grep('muij', names(pars3))]
mean.res.3 <- sapply(1:data$num_experiments, get.mean.res)

muijs <- pars3b[grep('muij', names(pars3b))]
mean.res.3b <- sapply(1:data$num_experiments, get.mean.res)

pdf(file='Figures/Fit3_Fit3b_Residuals.pdf', width=7, height=5)
par(mfrow=c(1, 2))

plot(mean.res.3, mean.res.3b, main='mean(abs(Residuals)) by Exp.\nFit3 vs. Fit3b')
abline(a=0, b=1, col='red')

plot(mean.res.3, mean.res.3b, xlim=c(1, 4), ylim=c(1, 4))
abline(a=0, b=1, col='red')
dev.off()

muijs <- pars3c[grep('muij', names(pars3c))]
mean.res.3c <- sapply(1:data$num_experiments, get.mean.res)

pdf(file='Figures/Fit3_Fit3c_Residuals.pdf', width=7, height=5)
par(mfrow=c(1, 2))

plot(mean.res.3, mean.res.3c, main='mean(abs(Residuals)) by Exp.\nFit3 vs. Fit3c')
abline(a=0, b=1, col='red')

plot(mean.res.3, mean.res.3c, xlim=c(1, 4), ylim=c(1, 4))
abline(a=0, b=1, col='red')
dev.off()

## ridges -----

library(tidyverse)
library(ggridges)

ev3 <- read_tsv('dat/ev.adj.Fit3.txt')
ev3b <- read_tsv('dat/ev.adj.Fit3b.txt')
ev3c <- read_tsv('dat/ev.adj.Fit3c.txt')

dpep.eq <- parse(text=paste0('log[10](frac(\'Spectral PEP\',\'Updated PEP\'))'))

breaks <- c(0, 5e-2, 1e-1, 5e-1, 1)

ev.f <- ev3 %>%
  select(c('Sequence', 'PEP', 'PEP.new', 'Retention time', 'muijs')) %>%
  filter(!is.na(PEP.new)) %>%
  filter(PEP > 0 & PEP.new > 0) %>%
  mutate_at('PEP', funs(ifelse(. > 1, 1, .))) %>%
  mutate(dPEP=log10(PEP/PEP.new)) %>%
  arrange(desc(PEP)) %>%
  mutate(bin=cut(PEP, breaks=10))

pdf('Figures/Fit3c_Ridges.pdf', width=7, height=7)

ggplot(ev.f) +
  geom_density_ridges(aes(x=dPEP, y=bin, group=bin), 
                      rel_min_height=0.01) +
  geom_vline(xintercept=0, color='red', linetype='longdash') +
  annotate(geom='text', label='Decreased Confidence', x=-1.5, y=11.5, size=5) +
  annotate(geom='text', label='Increased Confidence', x=1.5, y=11.5, size=5) +
  scale_x_continuous(limits=c(-2.5,2.5), expand=c(0.01, 0)) +
  scale_y_discrete(expand=(c(0.01, 0)), position='right') +
  # scale_y_continuous(breaks=seq(1,length(levels(ev.f$bin))),
  #                   labels=levels(ev.f$bin),
  #                   expand=c(0.01, 0), position='right',
  #                   sec.axis=sec_axis(~., labels=NULL, name='Density')) +
  labs(x=dpep.eq, y='Spectral PEP Bin', title='Change in PEP by Spectral PEP - Fit 3c') +
  #labs(x='dPEP', y=NULL) +
  theme_ridges() +
  theme(axis.title.x = element_text(hjust=0.5),
        axis.title.y = element_text(hjust=0.5))

dev.off()

pdf('Figures/Fit3_Ridges_dRT.pdf', width=7, height=7)

ggplot(ev.f) +
  geom_density_ridges(aes(x=log2(abs(muijs-`Retention time`)), y=bin, group=bin), 
                      rel_min_height=0.01) +
  geom_vline(xintercept=0, color='red', linetype='longdash') +
  #annotate(geom='text', label='Decreased Confidence', x=-1.5, y=11.5, size=5) +
  #annotate(geom='text', label='Increased Confidence', x=1.5, y=11.5, size=5) +
  scale_x_continuous(expand=c(0.01, 0)) +
  scale_y_discrete(expand=(c(0.01, 0)), position='right') +
  # scale_y_continuous(breaks=seq(1,length(levels(ev.f$bin))),
  #                   labels=levels(ev.f$bin),
  #                   expand=c(0.01, 0), position='right',
  #                   sec.axis=sec_axis(~., labels=NULL, name='Density')) +
  labs(x='log2(abs(dRT))', y='Spectral PEP Bin', title='log2(abs(dRT)) by Spectral PEP - Fit 3') +
  #labs(x='dPEP', y=NULL) +
  theme_ridges() +
  theme(axis.title.x = element_text(hjust=0.5),
        axis.title.y = element_text(hjust=0.5))

dev.off()


## validation ----




