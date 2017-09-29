## analyze effectiveness of Multiple Alignment method
# init
library(readr)
library(ggplot2)
library(reshape2)
library(hexbin)
library(pracma)
library(qvalue) ## usage here: http://svitsrv25.epfl.ch/R-doc/library/qvalue/html/qvalue.html
library(dplyr)
library(RColorBrewer)

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
qval <- qvalue(ev.f$PEP.updated, lambda=0.05, pi0.method='smoother')
#qval <- qvalue(sort(ev.f$PEP.updated))
ev.f$q[order(ev.f$PEP.updated)] <- cumsum(sort(ev.f$PEP.updated)) / seq(1,nrow(ev.f))
ev.f$q.old[order(ev.f$PEP)] <- cumsum(sort(ev.f$PEP)) / seq(1,nrow(ev.f))

inds <- seq(1,nrow(ev.f),by=100)

pal = brewer.pal(4, 'Set1')
ggplot(ev.f[inds,], aes(x=inds)) +
  geom_path(aes(y=sort(PEP), color='PEP')) +
  geom_path(aes(y=sort(PEP.updated), color='PEP.new')) +
  geom_path(aes(y=sort(q), color='FDR.new')) +
  geom_path(aes(y=sort(q.old), color='FDR')) +
  scale_color_manual(values=pal) +
  labs(title='RTLib: Updated PEP and FDR')

ggplot(ev.f[inds,], aes(x=inds)) +
  geom_path(aes(y=sort(qval$qvalues)))

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


# See the effects of changing PSM RT standard deviations (sigmas) -------

ev.f.0.5 <- read_tsv('dat/ev.adj.elite_sigmas0.5.txt')
ev.f.1 <- read_tsv('dat/ev.adj.elite.txt')
ev.f.1.5 <- read_tsv('dat/ev.adj.elite_sigmas1.5.txt')
ev.f.2 <- read_tsv('dat/ev.adj.elite_sigmas2.txt')
ev.f.3 <- read_tsv('dat/ev.adj.elite_sigmas3.txt')

# combine data sets, and only take what we need
cols <- c('Sequence', 'Raw file', 'Retention time', 'PEP', 'rt.minus', 'rt.plus', 'muijs',
          'sigmas', 'PEP.new')

ev.f.0.5 <- ev.f.0.5[,cols]
ev.f.1 <- ev.f.1[,cols]
ev.f.1.5 <- ev.f.1.5[,cols]
ev.f.2 <- ev.f.2[,cols]
ev.f.3 <- ev.f.3[,cols]

# factorize data sets so we can separate the data later
ev.f.0.5$sigma.mod <- 0.5
ev.f.1$sigma.mod <- 1
ev.f.1.5$sigma.mod <- 1.5
ev.f.2$sigma.mod <- 2
ev.f.3$sigma.mod <- 3

# sorted PEP.new
ev.f.0.5$PEP.sort <- sort(ev.f.0.5$PEP.new, na.last=FALSE)
ev.f.1$PEP.sort <- sort(ev.f.1$PEP.new, na.last=FALSE)
ev.f.1.5$PEP.sort <- sort(ev.f.1.5$PEP.new, na.last=FALSE)
ev.f.2$PEP.sort <- sort(ev.f.2$PEP.new, na.last=FALSE)
ev.f.3$PEP.sort <- sort(ev.f.3$PEP.new, na.last=FALSE)

# combine data sets
ev.f <- rbind(ev.f.0.5, ev.f.1, ev.f.1.5, ev.f.2, ev.f.3)
# factorize
ev.f$sigma.mod <- factor(ev.f$sigma.mod)
ev.f$index <- rep(seq(1,nrow(ev.f.1)), 5)

x = seq(1,nrow(ev.f),by=100)
ggplot(ev.f[x,], aes(x=index, y=PEP.sort, color=sigma.mod)) +
  geom_path() +
  scale_x_continuous(limits=c(1.2e6,1.4e6), breaks=c()) +
  labs(x='Index (truncated)', y='sort(PEP.new)', color='RT std (sigma) Weight',
       title='RTLib:\nEffect of RT standard deviation weight\non updated PEP (PEP.new)')


## fold change between methods -----
source('parse.ev.adj.R')
source('adjust.pep.ali.R')
ev <- parse.ev.adj('dat/evidence+dRT.elite.txt', type='Ali')
#ev.ali <- parse.ev.adj('dat/evidence+dRT.elite.txt', type='Ali')
ev.pep <- parse.ev.adj('dat/ev.adj.elite.txt')
#ev.exp <- parse.ev.adj('dat/ev+dRT.elite.txt', type='Ali')
ev.exp <- adjust.pep.ali(path.out=NULL)
#ev.exp.05 <- adjust.pep.ali(path.out=NULL)

acc <- 100
t <- logseq(1e-3, 1e-1, n=acc)
#id.frac <- c(length=acc*2)
id.frac <- c(length=acc*3)
counter <- 0
for(i in t) {
  counter <- counter + 1
  #id.frac[counter] <- sum(ev$PEP.updated < i) / sum(ev$PEP < i)
  #id.frac[acc+counter] <- sum(ev.pep$PEP.updated < i) / sum(ev.pep$PEP < i)
  #id.frac[(2*acc)+counter] <- sum(ev.exp.05$PEP.updated < i) / sum(ev.exp.05$PEP < i)
  id.frac[counter] <- sum(ev$PEP.new < i, na.rm=T) / sum(ev$PEP < i & !is.na(ev$PEP.new))
  id.frac[acc+counter] <- sum(ev.pep$PEP.new < i, na.rm=T) / sum(ev.pep$PEP < i & !is.na(ev.pep$PEP.new))
  id.frac[(2*acc)+counter] <- sum(ev.exp$PEP.new < i, na.rm=T) / sum(ev.exp$PEP < i & !is.na(ev.exp$PEP.new))
  cat(counter, '/', acc, '\r')
  flush.console()
}

df <- data.frame(
  PEP=as.numeric(id.frac),
  Method=as.factor(rep(c('Experiment-Centric (Ali)', 'Peptide-Centric (STAN)',
                         'Experiment-Centric (STAN)'), each=acc))
)

ggplot(df, aes(x=rep(t,3), y=PEP, color=Method, fill=Method)) +
  geom_point(size=1) + 
  geom_path() +#geom_smooth(method='loess') +
  scale_x_log10() +
  annotation_logticks(sides='b') +
  #scale_y_continuous(limits=c(0.95,2.25), breaks=c(1, 1.25, 1.5, 1.75, 2, 2.25)) +
  labs(x='PEP Threshold', y='Fold Change Increase in IDs',
       title=paste0('Fold Change Increase of PSM IDs\n',
                    '= #Adjusted PEPs / #Original PEPs above PEP Threshold'))


## look at subset of confident PSMs from spectra alone (low PEP) that
## get downgraded by our method
ev <- ev.new

# PEP from  < 1e-5 -> > 1e-2
ev.a <- subset(ev.new, PEP < 1e-5 & PEP.updated > 1e-2)

plot.rt.dist <- function(PSM) {
  sequence <- PSM$Sequence
  x <- ev[ev$Sequence == sequence,]$RT.new
  h <- hist(x, breaks=30,
            main='',
            xlab='Retention Time')
  xfit<-seq(min(x),max(x),length=40) 
  yfit<-dnorm(xfit,mean=PSM$muijs,sd=PSM$sigmas) 
  yfit <- yfit*diff(h$mids[1:2])*length(x)
  lines(xfit, yfit, col="red", lwd=2)
  lines(rep(PSM$`Retention time`,2), c(0, 100), col='blue', lty=2)
  mtext(paste0(sequence, ' [', PSM$ID, ']\n',
               'CON/REV: ', 
               !length(grep('CON|REV', PSM$Proteins)) == 0,
               '\n','Mean=',
               format(PSM$muijs, digits=3, scientific=FALSE),
               ', Sd=', 
               format(PSM$sigmas, digits=3, scientific=FALSE)), 
        line=0.5,
        cex=0.75)
}

i <- 4
plot.rt.dist(ev.a[i,])

par(mfrow=c(2,4), mar=c(6, 3, 6, 3))
for(i in sample.int(nrow(ev.a), size=8)) {
  plot.rt.dist(ev.a[i,])
}
# 215779, 98762, ##1200724, 436462, 195181

par(mfrow=c(1,1))
plot.rt.dist(ev.new[215779,])

par(mfrow=c(2,2), mar=c(4, 4, 4, 4))
for(i in c(215779, 98762, 436462, 195181)) {
  plot.rt.dist(ev.new[i,])
}

# investigate neighbors of 436462 (RT-wise)

psm <- ev.new[436462,]
ev.a <- subset(ev.new, `Raw file`==psm$`Raw file` & 
                 abs(RT.new - psm$RT.new) < 1 &
                 PEP < 1e-5)
ev.a <- ev.a[order(ev.a$`Retention time`),]
# 3, 4, 6, 8?, 9, 10, 11, 15, 16, 19, 20, 21, 22, 23
par(mfrow=c(1,1))
plot.rt.dist(ev.a[22,])

par(mfrow=c(2, 2))
for(i in c(6, 10, 16, 22)) {
  plot.rt.dist(ev.a[i,])
}

# KS-test between empirical and predicted distribution
# psms - subset of ev data
psm.ks.test <- function(psms) {
  apply(psms[,c('Sequence', 'muijs', 'sigmas')], MARGIN=1, FUN=function(psm) {
    x <- ev[ev$Sequence == psm[['Sequence']],]$RT.new
    suppressWarnings(
      a <- ks.test(x=x, y='pnorm', 
                   mean=as.numeric(psm[['muijs']]), 
                   sd=as.numeric(psm[['sigmas']]))
    )
    return(a[['statistic']])
  })
}

psm.ks.test(psms)

# Compare KS scores of confident PSMS (PEP < 1e-5) around
# potentially misaligned PSMs : 98762, 436462
# to global sample

ev.conf <- subset(ev, PEP < 1e-3)
ks.scores <- data.frame()
ks.scores <- rbind(ks.scores, cbind(
  as.numeric(psm.ks.test(ev.conf[sample.int(nrow(ev.conf), size=1000),])),
  'Global'
))
for(i in ev.a$ID) {
  cat(i, '\r')
  flush.console()
  psm <- ev.new[i,]
  ks.scores <- rbind(ks.scores, cbind(
    as.numeric(psm.ks.test(
      subset(ev.new, `Raw file`==psm$`Raw file` & 
               abs(RT.new - psm$RT.new) < 3 & # within 3 minutes
               PEP < 1e-3) # PEP < 1e-3
    )),
    as.character(i)
  ))
}
ks.score.freq <- as.data.frame(table(ks.scores$V2))
ks.scores$V1 <- as.numeric(as.character(ks.scores$V1))
ks.scores.a <- subset(ks.scores, V2=='Global' |
                        #V2 %in% as.character(ks.score.freq$Var1[ks.score.freq$Freq > 100]))
                        V2 %in% c('46585', '202339', '436462', '505784', '562059', '827651', '1006356',
                                  '1084307', '1200724'))
names(ks.scores.a) <- c('KS.Score', 'Peptide')
ggplot(ks.scores.a, 
       aes(x=Peptide, y=KS.Score)) + 
  geom_boxplot(aes(fill=Peptide)) +
  scale_fill_brewer(palette='Set3') +
  scale_color_brewer(palette='Set3') +
  labs(title=paste0('K-S Test between empirical and estimated RT distributions\n',
                    'PSMs taken from same experiment, and within 3 min RT\n',
                    'All PSMs have PEP < 1e-3'))
