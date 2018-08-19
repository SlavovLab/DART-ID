library(tidyverse)
library(pracma)
library(ggridges)
library(RColorBrewer)
source('Rscripts/lib.R')

#ev <- read_tsv("/gd/Slavov_Lab/Albert/RTLib_Alignments/SQC_1000cell_20180607_6/ev_updated.txt")
ev_a <- read_tsv('/gd/bayesian_RT/Alignments/FP_validation_set1_20180706_2/ev_updated.txt')
ev_b <- read_tsv('/gd/bayesian_RT/Alignments/FP_validation_set2_20180706_1/ev_updated.txt')
ev_all <- read_tsv('/gd/bayesian_RT/Alignments/FP_validation_allsets_20180706_1/ev_updated.txt')

conf_thresh <- 0.01

# find peptides in ev_a which are boosted but have no PSMs with PEP < conf_thresh
peps_a <- ev_a %>%
  filter(!is.na(pep_new)) %>%
  group_by(`Modified sequence`) %>%
  summarise(n=n(),
            rt=median(`Retention time`),
            rt_sd=sd(`Retention time`),
            spec=sum(PEP < conf_thresh),
            new=sum(PEP > conf_thresh & pep_new < conf_thresh)) %>%
  filter(spec == 0) %>%
  filter(new > 0) %>%
  arrange(`Modified sequence`)

# find peptides in ev_b which have confident PSMs
peps_b <- ev_b %>%
  group_by(`Modified sequence`) %>%
  summarise(n=n(),
            rt=median(`Retention time`),
            rt_sd=sd(`Retention time`),
            spec=sum(PEP < conf_thresh)) %>%
  filter(spec > 0) %>%
  arrange(`Modified sequence`)

common_peps <- intersect(peps_a$`Modified sequence`, peps_b$`Modified sequence`)

# ev_c <- data.frame(
#   `Modified sequence`=common_peps,
#   rt_a=as.numeric(peps_a$rt[peps_a$`Modified sequence` %in% common_peps]),
#   rt_b=as.numeric(peps_b$rt[peps_b$`Modified sequence` %in% common_peps]),
#   rt_sd_a=as.numeric(peps_a$rt_sd[peps_a$`Modified sequence` %in% common_peps]),
#   rt_sd_b=as.numeric(peps_b$rt_sd[peps_b$`Modified sequence` %in% common_peps])
# )

ev_c <- ev_a %>% 
  filter(`Modified sequence` %in% common_peps) %>%
  dplyr::select(c('Modified sequence', 'Retention time'))

ev_c$rt_conf <- peps_b$rt[match(ev_c$`Modified sequence`, peps_b$`Modified sequence`)]

## ----------

pdf(file='manuscript/Figs/rt_validation_v6.pdf', width=3.5, height=1.75)

#layout(rbind(c(1,1),
#             c(2,3)))

par(pty='s', mar=c(1.75,1.25,0.25,0.5), cex.axis=0.85)

#plot(ev_c$rt_b, ev_c$rt_a, pch=16, cex=0.75,
plot(ev_c$rt_conf, ev_c$`Retention time`, pch=16, cex=0.4,
     xlim=c(12, 52), ylim=c(12, 52),
     xaxt='n', yaxt='n', xlab=NA, ylab=NA)
#errorbar(ev_c$rt_b, ev_c$rt_a, ev_c$rt_sd_b, ev_c$rt_sd_a, 
#         col='red', add=T)
abline(a=0, b=1, col='red')

text(#11.5, 48, adj=c(0, 0.5)
     52, 16, adj=c(1, 0.5),
     labels=paste0('n = ', nrow(ev_c), ' PSMs\n',
                   length(common_peps), ' shared peptides'),
     cex=0.7)

axis(1, at=seq(10, 50, by=10), tck=-0.02, mgp=c(0, 0.01, 0))
axis(2, at=seq(10, 50, by=10), tck=-0.02, mgp=c(0, 0.3, 0), las=1)

mtext('Set B: Median RT (min)', 
      side=1, line=0.85, cex=0.8)
mtext('Set A: RT (min)', 
      side=2, line=1.25, cex=0.8)
#mtext('PSM Validation with Technical Replicates', side=3, line=0.25)

# legend('bottomright', c('Median RT', 'Std RT'), 
#        pch=c(16, NA), lty=c(NA, 1), lwd=c(NA, 2), col=c('black', 'red'),
#        cex=0.9, pt.cex=1.25,
#        #bty='n', 
#        y.intersp=1.25, inset=c(0, 0))

dev.off()

## ----------

#par(cex.axis=0.85, pty='m', 
#    mar=c(2, 2.5, 2.5, 0.5))

# find peptides in ev_a which are boosted but have no PSMs with PEP < conf_thresh
peps_a <- ev_a %>%
  filter(!is.na(pep_new)) %>%
  group_by(`Modified sequence`) %>%
  summarise(n=n(),
            rt=median(`Retention time`),
            rt_sd=sd(`Retention time`),
            spec=sum(PEP < conf_thresh),
            new=sum(PEP > conf_thresh & pep_new < conf_thresh)) %>%
  filter(spec == 0) %>%
  filter(new > 0) %>%
  arrange(`Modified sequence`)

# find peptides in ev_b which have confident PSMs
peps_b <- ev_b %>%
  group_by(`Modified sequence`) %>%
  summarise(n=n(),
            rt=median(`Retention time`),
            rt_sd=sd(`Retention time`),
            spec=sum(PEP < conf_thresh)) %>%
  filter(spec > 0) %>%
  arrange(`Modified sequence`)

common_peps <- intersect(peps_a$`Modified sequence`, peps_b$`Modified sequence`)

ev_c <- ev_a %>% 
  filter(`Modified sequence` %in% common_peps) %>%
  dplyr::select(c('Modified sequence', 'Retention time'))
ev_c$rt_conf <- peps_b$rt[match(ev_c$`Modified sequence`, peps_b$`Modified sequence`)]
ev_c$set <- 'new'

# confident
common_peps <- intersect(
  unique(ev_a %>% filter(!is.na(pep_new)) %>% filter(PEP < conf_thresh) %>% pull(`Modified sequence`)), 
  unique(ev_b %>% filter(!is.na(pep_new)) %>% filter(PEP < conf_thresh) %>% pull(`Modified sequence`)))
rts_a <- ev_a %>% 
  filter(!is.na(pep_new)) %>%
  filter(PEP < conf_thresh) %>% 
  filter(`Modified sequence` %in% common_peps) %>%
  arrange(`Modified sequence`) %>%
  dplyr::select(c('Modified sequence', 'Retention time'))
rts_b <- ev_b %>% 
  filter(!is.na(pep_new)) %>%
  filter(PEP < conf_thresh) %>% 
  filter(`Modified sequence` %in% common_peps) %>%
  group_by(`Modified sequence`) %>%
  summarise(rt=median(`Retention time`)) %>%
  arrange(`Modified sequence`)
rts_a$rt_conf <- rts_b$rt[match(rts_a$`Modified sequence`, rts_b$`Modified sequence`)]
rts_a$set <- 'confident'
ev_c <- rbind(ev_c, rts_a)

# scramble
rts <- ev_all %>%
  filter(!is.na(pep_new)) %>%
  filter(PEP < conf_thresh) %>%
  dplyr::select(c('Modified sequence', 'Retention time', 'muij'))
rts$rt_conf <- sample(rts$muij, size=nrow(rts))
rts <- rts[,-3]
rts$set <- 'null'
ev_c <- rbind(ev_c, rts)

ev_c <- ev_c %>%
  mutate(drt=`Retention time`-rt_conf)


# horizontal boxplot ------------------------------------------------------

boxs <- list(
  null=log10(abs(ev_c$drt[ev_c$set=='null'])),
  new=log10(abs(ev_c$drt[ev_c$set=='new'])),
  conf=log10(abs(ev_c$drt[ev_c$set=='confident']))
)


# a -----------------------------------------------------------------------

pdf(file='manuscript/Figs/rt_validation_boxplot.pdf', width=3.5, height=1.75)
par(mar=c(2.5,5.5,0.5,2), cex.axis=0.85)

boxplot(boxs, horizontal=T, xlab=NA, ylab=NA,
        #col=c(av[4], av[2], av[1]),
        xaxt='n', yaxt='n',
        outpch=16, outcol=rgb(0,0,0,0.1), outcex=0.5)
axis(1, at=seq(-16, 10, by=1), tck=-0.02, mgp=c(0, 0.2, 0))
axis(2, at=1:3, labels=c('Random', 'Set A vs. Set B', 'Confident'), las=1, tck=-0.02,
     mgp=c(0, 0.5, 0))
mtext('log10 | Residual RT |  (min)', side=1, cex=1, line=1.25)
#mtext('')

dev.off()

# density on one line -----------------------------------------------------

cols <- brewer.pal(3, 'Set1')[c(2, 3, 1)]

p <- ggplot(ev_c) +
  #geom_density(aes(drt, fill=set), alpha=0.5) +
  stat_density(aes(drt, color=set), geom='line') +
  stat_density(aes(drt, fill=set), geom='area', alpha=0.3) +
  geom_abline(slope=0, intercept=0) +
  #geom_histogram(aes(drt, y=..density.., fill=set), bins=80) +
  scale_x_continuous(limits=c(-5, 5), expand=c(0,0)) +
  scale_y_continuous(limits=c(0, 7), expand=c(0,0)) +
  scale_color_manual(values=cols, guide=F) +
  scale_fill_manual(values=cols, guide=F) +
  labs(x='Residual RT (min)', y='Density') +
  theme_bert()


# no ggplot ---------------------------------------------------------------

cols <- brewer.pal(3, 'Set1')[c(2, 3, 1)]

plot(0, 0, type='n',
     xlim=c(-2, 2), ylim=c(0, 3.5), xlab=NA, ylab=NA, 
     xaxt='n', yaxt='n', xaxs='i', yaxs='i')
#x <- seq(-10, 10, by=0.01)
den1 <- density(ev_c$drt[ev_c$set=='null'], n=1000, adjust=2)
den2 <- density(ev_c$drt[ev_c$set=='confident'], n=1000, adjust=3)
den3 <- density(ev_c$drt[ev_c$set=='new'], n=1000, adjust=2)
lines(den1$x, den1$y, col=cols[3])
lines(den2$x, den2$y, col=cols[1])
lines(den3$x, den3$y, col=cols[2])
polygon(c(-5, den1$x, 5), c(0, den1$y, 0), border=NA, col=paste0(cols[3], '66'))
polygon(c(-5, den2$x, 5), c(0, den2$y, 0), border=NA, col=paste0(cols[1], '66'))
polygon(c(-5, den3$x, 5), c(0, den3$y, 0), border=NA, col=paste0(cols[2], '66'))

# violin ------------------------------------------------------------------

ggplot(ev_c) +
  #geom_density_ridges(aes(drt, set, group=set), bandwidth=0.02) +
  geom_histogram(aes(drt, y=..density.., fill=set), bins=80) +
  facet_grid(set~.) +
  labs(x='Residual RT', y='Density') +
  scale_x_continuous(limits=c(-2, 2))

#dev.off()

## -----------

# cor mat between confident RTs in the technical replicate sets
common_peps <- ev_all %>%
  filter(PEP < 0.01) %>%
  group_by(`Modified sequence`) %>%
  summarise(n=length(unique(`Raw file`))) %>%
  filter(n==11) %>%
  pull(`Modified sequence`)

exps <- unique(ev_all$`Raw file`)
cormat <- zeros(length(exps), length(exps))

for(i in 1:length(exps)) {
  rt_i <- ev_all %>% 
    filter(`Raw file`==exps[i]) %>% 
    filter(`Modified sequence` %in% common_peps) %>%
    group_by(`Modified sequence`) %>%
    summarise(rt=median(`Retention time`)) %>%
    arrange(`Modified sequence`) %>%
    pull(rt)
  
  for(j in 1:length(exps)) {
    if(i == j) {
      cormat[i, j] = NA
      next
    }
    rt_j <- ev_all %>% 
      filter(`Raw file`==exps[j]) %>% 
      filter(`Modified sequence` %in% common_peps) %>%
      group_by(`Modified sequence`) %>%
      summarise(rt=median(`Retention time`)) %>%
      arrange(`Modified sequence`) %>%
      pull(rt)
    
    cormat[i, j] <- cor(rt_i, rt_j)
  }
}


## ------------

pdf(file='manuscript/Figs/tech_replicate_rts.pdf', width=3.5, height=3)

layout(rbind(c(1, 2)), widths=c(4, 1))

par(mar=c(2,2,1,0.15), cex.axis=0.75)

ncols <- 50
colfunc <- colorRampPalette(c('white', 'red'))
diag(cormat) <- 1-(1e-5)

image(-log10(1-cormat)[nrow(cormat):1,], zlim=c(2, 5), col=colfunc(ncols),
      xaxt='n', yaxt='n', xlab=NA, ylab=NA)

alabels <- LETTERS[seq(1, length(exps))]
#alabels[seq(2,length(alabels),2)] <- NA
axis(1, tck=-0.02, at=seq(0,1, length.out=length(exps)),
     labels=alabels, mgp=c(0,0.2,0))
axis(2, tck=-0.02, at=seq(0,1, length.out=length(exps)),
     labels=rev(alabels), mgp=c(0,0.5,0), las=1)

mtext('Experiment', side=1, line=1.05)
mtext('Experiment', side=2, line=1.15)
mtext('RT Correlation in Technical Replicates', side=3, line=0.15, font=2)

par(mar=c(2,0.15,1,2.5), cex.axis=0.75)

image(matrix(seq(2, 5, length.out=ncols), ncol=ncols), col=colfunc(ncols),
      xlab=NA, ylab=NA, xaxt='n', yaxt='n', useRaster=T)
#mtext('Correlation', side=3, line=0.15)
axis(4, tck=-0.02, at=seq(0, 1, length.out=4), 
     labels=1-(10^(-seq(2, 5, length.out=4))), mgp=c(0, 0.2, 0), las=1)

dev.off()

## plot the downgraded ones? -----------

conf_thresh <- 1e-2

peps <- ev_all %>%
  filter(!is.na(pep_new)) %>%
  mutate(confident=(PEP < conf_thresh),
         downgraded=(PEP < conf_thresh & pep_new > conf_thresh),
         upgraded=(PEP > conf_thresh & pep_new < conf_thresh)) %>%
  group_by(`Modified sequence`) %>%
  summarise(n=n(),
            confident=sum(confident),
            downgraded=sum(downgraded),
            upgraded=sum(upgraded)) %>%
  filter(confident > 0) %>%
  filter(confident != downgraded & confident != upgraded) %>%
  #filter(downgraded > 0 | upgraded > 0) %>%
  filter(downgraded > 0) %>%
  arrange(`Modified sequence`)

#sapply(peps$`Modified sequence`, function(s) {
#  ev_f <- ev_all %>% filter(`Modified sequence`==s)
#  list(ev_f %>% filter(PEP < conf_thresh & pep_new < conf_thresh) %>% pull(`Retention time`))
#})

## --------

ev_f <- ev_all %>% 
  filter(`Modified sequence` %in% peps$`Modified sequence`) %>%
  filter(PEP < conf_thresh & pep_new < conf_thresh) %>%
  dplyr::select(c('Modified sequence', 'Retention time', 'PEP', 'pep_new', 'muij')) %>%
  mutate(sequence=as.factor(`Modified sequence`))
  
stripchart(`Retention time`~sequence, data=ev_f, col=rainbow(nrow(peps)), 
           pch=1, vertical=T, xaxt='n', 
           xlab='Peptide Sequence', ylab='Retention time (min)',
           main='Downgraded PSMs\no = Confident, x = Downgraded')

ev_f <- ev_all %>% 
  filter(`Modified sequence` %in% peps$`Modified sequence`) %>%
  filter(PEP < conf_thresh & pep_new > conf_thresh) %>%
  dplyr::select(c('Modified sequence', 'Retention time', 'PEP', 'pep_new', 'muij', 'sigmaij')) %>%
  mutate(sequence=as.factor(`Modified sequence`))

stripchart(`Retention time`~sequence, data=ev_f, col=rainbow(nrow(peps)), 
           pch=4, vertical=T, add=T, cex=1.25)

abline(v=seq(1,nrow(peps)), col=rainbow(nrow(peps), alpha=0.3))

axis(1, at=seq(0,nrow(peps), by=10), tck=-0.03)

## ------------
# 180413S_X_FP18E
ev4 <- ev_all %>% 
  filter(exp_id == 4) %>% 
  filter(!is.na(pep_new)) %>%
  dplyr::select(c('Retention time', 'muij', 'PEP', 'pep_new', 'sigmaij')) %>%
  mutate(sigma_ms=sigmaij*60*1000,
         drt=abs(`Retention time`-muij))

sigma_levels <- cut(ev4$sigma_ms, breaks=12)

plot(ev4$PEP, ev4$pep_new, log='xy', 
     xlim=c(1e-8, 1), ylim=c(1e-8, 1),
     col=heat.colors(12, alpha=0.8)[as.numeric(sigma_levels)],
     pch=16,
     xlab='Spectra PEP', ylab='DART-ID PEP')
abline(a=0, b=1, col='black')

# plot(ev4$PEP, ev4$pep_new, log='xy', 
#      xlim=c(1e-8, 1), ylim=c(0.4, 1),
#      col=heat.colors(12, alpha=0.8)[as.numeric(sigma_levels)],
#      pch=16,
#      xlab='Spectra PEP', ylab='DART-ID PEP')
# abline(a=0, b=1, col='black')

sum(log10(ev4$pep_new / ev4$PEP) > 0) / nrow(ev4)

## Function to plot color bar ----------
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  #dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}
color.bar(heat.colors(12), min=50, max=100)
mtext('Sigma (milliseconds)', side=2, line=2.5)

## ----------

plot(density(-log10(ev4$PEP), adjust=3), col='black',
     xlim=c(0, 20), xlab='Confidence: -log10(PEP)', lwd=2,
     main=NA)
lines(density(-log10(ev4$pep_new), adjust=3), col='red', lwd=2)
legend('topright', c('Spectra PEP', 'DART-ID PEP'), lty=1, lwd=2, col=c('black', 'red'))

## --------

hist(ev4$sigma_ms, breaks=20, xlab='Sigma (milliseconds)', main=NA)
