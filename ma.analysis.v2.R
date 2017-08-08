## analyze effectiveness of Multiple Alignment method
# init
library(readr)
library(ggplot2)
library(reshape2)
library(hexbin)
library(pracma)
library(qvalue) ## usage here: http://svitsrv25.epfl.ch/R-doc/library/qvalue/html/qvalue.html
library(dplyr)

# load data
#ev <- read_tsv('dat/ev.adjusted.txt')
ev <- read_tsv('dat/ev.adj.elite.txt')
ev.f <- ev[, c('Sequence', 'Proteins', 'Raw file', 
               'Retention time','PEP', 'rt.minus', 'rt.plus',
               'muijs', 'sigmas', 'PEP.new')]

## calculate some important features, and fill in some missing data

# PEP > 1 defaults to PEP = 1
ev.f$PEP[ev.f$PEP > 1] <- 1

# PEP.updated - PEP.new except when its NA, and then defaults to PEP
ev.f$PEP.updated <- ev.f$PEP.new
ev.f$PEP.updated[is.na(ev.f$PEP.new)] <- ev.f$PEP[is.na(ev.f$PEP.new)]
# PEP == 0 defaults to PEP = min(PEP)
ev.f$PEP.updated[ev.f$PEP.updated <= 0] <- min(ev.f$PEP.updated)

ev.f$dPEP <- ev.f$PEP.new - ev.f$PEP
ev.f$dRT <- ev.f$muijs - ev.f$`Retention time`

# calculate q values
#qval <- qvalue(ev.f$PEP.updated, lambda=0.05, pi0.method='smoother')
qval <- qvalue(ev.f$PEP.updated)
ev.f$q <- qval$qvalues

inds <- seq(1,nrow(ev.f),by=100)

ggplot(ev.f[inds,], aes(x=inds)) +
  geom_path(aes(y=sort(p.adjust(PEP.updated, method='fdr'))), color='red') +
  geom_path(aes(y=sort(p.adjust(PEP, method='fdr'))), color='magenta') +
  geom_path(aes(y=sort(PEP)), color='blue') +
  geom_path(aes(y=sort(PEP.updated)), color='cyan')

df <- data.frame(
  x=as.numeric(inds),
  PEP=as.numeric(sort(ev.f$PEP[inds])),
  PEP.new=as.numeric(sort(ev.f$PEP.updated[inds])),
  FDR=as.numeric(sort(p.adjust(ev.f$PEP[inds], method='fdr'))),
  FDR.new=as.numeric(sort(p.adjust(ev.f$PEP.updated[inds], method='fdr')))
)

dfm <- melt(df, id.vars='x')
ggplot(dfm, aes(x=x, y=value, color=variable)) +
  geom_path()

ggplot(dfm %>% filter(variable %in% c('FDR', 'FDR.new')), aes(x=x, y=value*100, color=variable)) +
  geom_path(size=1) +
  scale_x_continuous(breaks=c()) +
  labs(title='FDR Comparison\nRTLib PEP Adjustment', 
       x='Index', y='FDR (%)')

ggplot(dfm, aes(x=x, y=value, color=variable)) +
  geom_path(size=1) +
  scale_x_continuous(breaks=c()) +
  labs(title='RTLib PEP Adjustment\nFDR, PEP Comparison', 
       x='Index', y='')

# fold change difference of FDR
t <- seq(0.01, 1, by=0.01)
y <- vector(length=length(t))
counter <- 0
for(tt in t) {
  counter <- counter + 1
  y[counter] <- sum(df$FDR.new < tt) / sum(df$FDR < tt)
}

ggplot(data.frame(), aes(x=t*100, y=y)) +
  geom_path(color='red') +
  geom_point(color='black', shape=4) +
  labs(x='FDR (%)', y='Fold Change',
       title='RTLib PEP Adjustment\nFold Change in # PSMs < FDR')
