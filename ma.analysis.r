## analyze effectiveness of Multiple Alignment method
# init
library(readr)
library(ggplot2)
library(reshape2)
library(hexbin)

# load data
#ev <- read_tsv('dat/ev.adjusted.txt')
ev <- read_tsv('dat/ev.adj.elite.txt')
ev.f <- ev[, c('Sequence', 'Proteins', 'Raw file', 
               'Retention time','PEP', 'rt.minus', 'rt.plus',
               'muijs', 'sigmas', 'PEP.new')]

# set original PEP > 1 to PEP = 1
ev.f$PEP[ev.f$PEP > 1] = 1

# plot distributions of PEP and PEP.new
par(mfcol=c(2, 1))
with(ev.f, hist(PEP[PEP < 0.05], breaks=seq(0, 0.05, by=1e-3),
  main=paste('Original PEP\n', 
             'PEP > 0.05:', sum(ev$PEP > 0.05), '(',
             formatC(sum(ev$PEP > 0.05) / length(ev$PEP) * 100, digits=4),'%)',
             '\nPEP < 0.05:', sum(ev$PEP < 0.05), '(',
             formatC(sum(ev$PEP < 0.05) / length(ev$PEP) * 100, digits=4),'%)'),
  xlab='Original PEP'))

with(ev.f, hist(PEP.new[PEP.new < 0.05], breaks=seq(0, 0.05, by=1e-3),
  main=paste('Adjusted PEP\n', 
             'PEP > 0.05:', sum(ev.f$PEP.new > 0.05), '(',
             formatC(sum(ev.f$PEP.new > 0.05) / length(ev.f$PEP.new) * 100, digits=4),'%)',
             '\nPEP < 0.05:', sum(ev.f$PEP.new < 0.05), '(',
             formatC(sum(ev.f$PEP.new < 0.05) / length(ev.f$PEP.new) * 100, digits=4),'%)'),
  xlab='Adjusted PEP'))

## dPEP (PEP - PEP.new)
ev.f <- cbind(ev.f, ev.f$PEP.new - ev.f$PEP)
names(ev.f)[11] <- 'dPEP'

# change in PEP by experiment
dpep.exp <- aggregate(dPEP ~ `Raw file`, data=ev.f, FUN=mean)$dPEP
hist(dpep.exp, breaks=30)

# dPEP overall
hist(ev.f$dPEP, breaks=50, xlab='PEP.new - PEP', 
     main=paste('Change in PEP (PEP.new - PEP)\n', 
                '-1 = completely incorrect -> completely correct\n',
                '1 = completely correct -> completely incorrect') )

# for these complete reversals (abs(dPEP) == 1), what is the difference in retention time?
ev.f <- cbind(ev.f, ev.f$`Retention time` - ev.f$muijs)
names(ev.f)[12] <- 'dRT'
hist(ev.f$dRT, main='Difference in RT (dRT) for all PSMs') # for all

# can we correlate dRT and dPEP?
ndPEP <- (ev.f$dPEP + 1)/2 # normalize on a scale of 0 -> 1

# sample 10,000 at a time for sanity
c <- sample.int(nrow(ev.f), size=10000)
# run linear regression
mod.1 <- lm(ndPEP[c] ~ ev.f$dRT[c])

plot(resid(mod.1) ~ fitted(mod.1),
  xlab = "Fitted Values",
  ylab = "Residuals",
  main = "Residual Diagnostic Plot"
)
plot(ndPEP[c], ev.f$dRT[c])
ggplot() +
  geom_point(aes(x=ndPEP[c], y=ev.f$rt.plus[c]), alpha=0.1) + 
  scale_x_log10()

# no way this is even close to correlating...

## look at dPEP for CON_ and REV_ proteins
inds.CON = grep('CON_', ev.f$Proteins)
inds.REV = grep('REV_', ev$`Leading razor protein`)

hist(ev.f$dPEP[inds.CON], breaks=100,
     main='Change in PEP for CON_ Contaminants',
     xlab='Change in PEP (PEP.new - PEP)')

hist(ev.f$dPEP[inds.REV], breaks=100,
     main='Change in PEP for REV_ Reverse Matches',
     xlab='Change in PEP (PEP.new - PEP)')

## look at dPEP for certain ranges of original PEP
# originally confident PEPs
dPEP <- ev.f$PEP.new - ev.f$PEP

# PEP < 0.05
inds.l005 <- ev.f$PEP < 0.05
hist(ev.f$dPEP[inds.l005], breaks=100,
     main='Change in PEP for PEP < 0.05',
     xlab='Change in PEP (PEP.new - PEP)')

# PEP < 5e-3
hist(ev.f$dPEP[ev.f$PEP < 5e-3], breaks=100,
     main='Change in PEP for PEP < 5e-3',
     xlab='Change in PEP (PEP.new - PEP)')

# PEP < 1e-4
hist(ev.f$dPEP[ev.f$PEP < 1e-4], breaks=100,
     main='Change in PEP for PEP < 1e-4',
     xlab='Change in PEP (PEP.new - PEP)')

## scatterplot original vs. adjusted PEPs
# we have 1.5M data points, so we either need density, hex bins, or contours

# try some random samples
c <- sample.int(nrow(ev.f), size=50000)
ggplot(ev.f[c,], aes(x=PEP, y=PEP.new)) + 
  geom_point(alpha=0.1) + 
  geom_abline(color='red') +
  scale_y_log10(limits=c(1e-10, 1), breaks=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1)) +
  scale_x_log10(limits=c(1e-10, 1), breaks=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1)) + 
  labs(title=paste0('PEP.new vs. PEP - 10,000 PSMs\n',
                    'Elite Experiments Only'))


inds.trash = (ev.f$PEP >= 1 & ev.f$PEP.new >= 1) | 
  seq(1,nrow(ev.f)) %in% inds.CON |
  seq(1,nrow(ev.f)) %in% inds.REV
ggplot(ev.f[(seq(1,nrow(ev.f)) %in% c) & !inds.trash,], aes(x=PEP, y=PEP.new)) +
  stat_density2d(geom='raster', aes(fill=..density..), contour=FALSE) +
  geom_abline(color='white') +
  scale_y_log10(limits=c(1e-10, 1)) +
  scale_x_log10(limits=c(1e-10, 1)) + 
  labs(title='PEP vs. PEP.new - 10,000 PSMs')

ggplot(ev.f[(seq(1,nrow(ev.f)) %in% c),], aes(x=PEP, y=PEP.new)) +
  geom_point(alpha=0.1) + 
  geom_abline(color='red') +
  scale_y_log10(limits=c(1e-3, 1)) +
  scale_x_log10(limits=c(1e-1, 1)) + 
  labs(title='PEP.new vs. PEP - 10,000 PSMs')

# plot fractional increase of confident IDs as a function of p (alpha)
# and plot fractional increase of FDRs as a function of %

lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

# plot ratios of PEP.new / PEP over a log scale of p values (t)
t <- lseq(1e-3, 1e-1, length.out=100)
id.frac <- c(length=100)
id.frac.drt <- c(length=100)
counter <- 0
for (i in t) {
  counter <- counter + 1
  id.frac[counter] <- (sum(ev.f$PEP.new < i) / sum(ev.f$PEP < i))
  id.frac.drt[counter] <- (sum(ev.drt$PEP.new < i) / sum(ev.drt$PEP < i))
}

# calculate FDR and FDR.new
# where FDR for peptide i is sum(PEPi) / sum(i) - sum of PEPs over # of PSMs
ev.peps <- unique(ev.f$Sequence)
ev.FDR <- aggregate(PEP ~ Sequence, data=ev.f, FUN=function(x) {
  sum(x) / length(x)
})$PEP
ev.FDR.new <- aggregate(PEP.new ~ Sequence, data=ev.f, FUN=function(x) {
  sum(x) / length(x)
})$PEP.new

ev.FDR.drt <- aggregate(PEP ~ Sequence, data=ev.drt, FUN=function(x) {
  sum(x) / length(x)
})$PEP
ev.FDR.drt.new <- aggregate(PEP.new ~ Sequence, data=ev.drt, FUN=function(x) {
  sum(x) / length(x)
})$PEP.new

# plot ratios of FDR.new / FDR over a log scale of % values (tq)
tq <- lseq(1e-3, 1e-1, length.out=100)
fdr.frac <- c(length=100)
fdr.frac.drt <- c(length=100)
counter <- 0
for (i in tq) {
  counter <- counter + 1
  fdr.frac[counter] <- (sum(ev.FDR.new < i) / sum(ev.FDR < i))
  fdr.frac.drt[counter] <- (sum(ev.FDR.drt.new < i) / sum(ev.FDR.drt < i))
}

df <- data.frame(PEP=as.numeric(rbind(id.frac, id.frac.drt)), 
                 FDR=as.numeric(rbind(fdr.frac, fdr.frac.drt)),
                 Method=as.factor(c('MA-STAN', 'dRT - Ali')), stringsAsFactors=FALSE)
df <- melt(df)

ggplot(df, aes(x=c(t, t, t*100, t*100), y=value)) +
  facet_wrap(Method~variable, nrow=2, scales='free') +
  geom_smooth(aes(color=variable)) +
  scale_x_log10() + 
  labs(x='PEP Threshold | FDR Threshold (%)', 
       y='Fractional Increase',
       title='Fractional Increase in PEP and FDR \nwith STAN Alignment')



##########################

#ev2 <- read_tsv('dat/ev.adjusted.2.txt')
ev2 <- read_tsv('dat/ev.adjusted.txt')
ev2.f <- ev2[, c('Sequence', 'Proteins', 'Raw file', 
               'Retention time','PEP', 'rt.minus', 'rt.plus',
               'muijs', 'sigmas', 'PEP.new')]
ev2.f$PEP[ev2.f$PEP > 1] = 1

ev.dPEP <- ev.f$PEP.new - ev.f$PEP
ev2.dPEP <- ev2.f$PEP.new - ev2.f$PEP
# account for NAs
ev2.dPEP[is.na(ev2.f$PEP.new)] <- 0
ev2.PEP.new <- ev2.f$PEP.new
ev2.PEP.new[is.na(ev2.f$PEP.new)] <- ev2.f$PEP[is.na(ev2.f$PEP.new)]

ev3 <- read_tsv('dat/evidence+dRT.txt')
ev3.dPEP <- ev3$PEP.new - ev3$PEP

d <- data.frame(
  dPEP=as.numeric(rbind(ev.dPEP, ev2.dPEP)),
  PEP.new=as.numeric(rbind(ev.f$PEP.new, ev2.f$PEP.new)),
  Method=as.factor(c('Filter at 0.05','No Filter'))
)

ggplot(d) +
  geom_density(aes(x=dPEP, y=..density.., fill=Method), alpha=0.3) +
  xlim(c(-1, 1))

ggplot(d, aes(x=PEP.new, y=..density.., fill=Method)) +
  geom_density(alpha=0.3)

ggplot(d, aes(x=PEP.new, fill=Method)) +
  geom_histogram(binwidth=3e-2, position='dodge')

## plot fold changes between old method (parameters based on PEP < 0.05) 
# and new method (parameters with all PSMs)

t <- lseq(1e-3, 1e-1, length.out=100)
id.frac <- c(length=300)
counter <- 0
for (i in t) {
  counter <- counter + 1
  id.frac[counter] <- (sum(ev.a$PEP.new < i) / sum(ev.a$PEP < i))
  id.frac[100+counter] <- (sum(ev2.f$PEP.new < i) / sum(ev2.f$PEP < i))
  id.frac[200+counter] <- (sum(ev3$PEP.new < i) / sum(ev3$PEP < i))
}
  
df <- data.frame(
  PEP=as.numeric(id.frac),
  Method=as.character(rep(c('Filter at 0.05','No Filter','Ali RT'),each=100))
)

ggplot(df, aes(x=rep(t, 3), y=PEP, color=Method, fill=Method)) +
  geom_smooth(method='loess') +
  scale_x_log10() +
  labs(x='PEP Threshold', y='Fold Change Increase in IDs',
       title=paste0('Fold Change Increase of PSM IDs\n',
                    '= #Adjusted PEPs / #Original PEPs above Threshold'))

## how many increased PEP (dPEP < 0) come from PSMs that come from sequences only
# identified once? (in the method w/ all PSMs forming the library)

# number of PSMs per sequence
ev.psms <- aggregate(dPEP ~ Sequence, data=ev.f, FUN=length)
names(ev.psms)[2] <- 'Count'

# smallest PEP for that sequence
ev.min.peps <- aggregate(dPEP ~ Sequence, data=ev.f, FUN=min)
#names(ev.min.peps)[2] <- 'dPEP'

# get all sequences with 1 PSM
i <- ev.psms$Count == 1
ev.sole.psms <- cbind(ev.psms[i,],ev.min.peps[i,]$dPEP)
names(ev.sole.psms)[3] <- 'dPEP'




## look at original PEPs below 1e-5 that get their PEP boosted by our method
# first create a subset of data for which we have updated params on
ev.a <- subset(ev.f, !is.na(ev.f$PEP.new))

# normalized dPEP
ev.a$PEP[ev.a$PEP == 0] <- 1e-100
ev.a$dPEP[ev.a$dPEP == 0] <- 1e-100
dPEP.norm <- ev.a$dPEP / ev.a$PEP
# look at 2-fold increases
inds <- ev.a$PEP < 1e-5 & dPEP.norm > 2
ev.twofold <- ev.a[inds,]
# look at general increases
inds <- ev.a$PEP < 1e-5 & dPEP.norm > 0
ev.increase <- ev.a[inds,]


# combine into one data set
a <- as.data.frame(rbind(
  cbind(ev.twofold$dRT, ev.twofold$rt.plus, '2-fold Increase'),
  cbind(ev.increase$dRT, ev.increase$rt.plus, 'Any Increase'),
  cbind(ev.a$dRT[ev.a$PEP < 1e-5], ev.a$rt.plus[ev.a$PEP < 1e-5], 'All')
), stringsAsFactors = FALSE)

names(a) <- c('dRT', 'rt.plus', 'Class')
a$Class <- factor(a$Class, levels=c('All', 'Any Increase', '2-fold Increase'))
a$dRT <- as.numeric(a$dRT)
a$rt.plus <- as.numeric(a$rt.plus)

ggplot(a, aes(x=abs(dRT), color=Class)) +
  geom_density(size=1) +
  scale_x_continuous(limits=c(0,70), breaks=seq(0,70,by=10)) +
  labs(title=paste0('abs(dRT) for PSMs with PEP < 1e-5'), 
       color='Increase of PEP.new')

## calculate FDR

ev.PEP <- ev.f$PEP
ev.PEP <- sort(ev.PEP)
