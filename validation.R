# exploratory script for validating results of RTLib method
# will turn into generalized method later

library(tidyverse)
library(reshape2)
library(stringr)

source('lib.R')
ev <- parse.ev.adj('dat/ev.adj.Fit2.txt')
# assign protein IDs
ev$Protein_ID <- as.numeric(as.factor(ev$Proteins))

ev$file <- clean.file.name(ev$`Raw file`)
ev$exp <- str_extract(ev$file, '[1-9][0-9][A-Z](\\d)?')

# load excel description
desc <- process.desc()
desc.types <- dcast(desc, Exp~Type, value.var='Type')
desc.types$Exp <- as.character(desc.types$Exp)

# choose experiments with at least one J and one H channel
exps <- desc.types[desc.types$J > 1 & desc.types$H > 1,'Exp']
ev <- ev[ev$exp %in% exps,]

data.cols <- grep('intensity', colnames(ev))
# ignore empty and carrier channels - we won't use them
#data.cols <- data.cols[-c(8, 10)]

# remove peptides w/o a protein
ev <- ev[!is.na(ev$Proteins),]
# remove non-unique peptides
ev <- ev[!grepl(';',ev$Proteins),]

# remove REV proteins
ev <- ev[!grepl('REV*',ev$Proteins),]
# remove CON proteins?
ev <- ev[!grepl('CON*',ev$Proteins),]
# remove PEP > 0.05
#ev <- ev[ev$PEP < 0.05,]

# take subset of experiments
#ev <- ev[grepl('30[A-J]|29[A-C]|26[A-E]|25[A-C]|24[A-C]', ev$exp),]
ev <- ev[grepl('19A', ev$exp),]

# set 0 to NA
ev[,data.cols][ev[,data.cols]==0] = NA

# remove rows without quantitation
#ev <- data[!(apply(data, 1, FUN=sum, na.rm=TRUE)==0),]
ev <- ev[!(apply(ev[,data.cols], MARGIN=1, FUN=sum, na.rm=TRUE)==0),]

## remove carrier and empty channels
## these channels vary between experiments, so we're gonna have to
## loop thru these
for(i in exps) {
  if(sum(ev$exp==i) <= 0) next
  inds.remove <- desc[desc$Exp==i & (desc$Type=='MIXED' | is.na(desc$Type)),'ch']
  # make relative to the ev data frame
  inds.remove <- inds.remove + data.cols[1] - 1
  # set to NA
  ev[ev$exp==i,inds.remove] <- NA
}

## normalize data

# first normalize by column, by median
# this assumes the same amount of protein per channel
# (not by cell, but in aggregate - this erases any ionization or 
# collision cell or any other quantitative biases)
ev[,data.cols] <- ev[,data.cols] / apply(ev[,data.cols], MARGIN=2, FUN=median, na.rm=TRUE)

# now normalize across rows, to get the difference between the channels
# AKA the difference between different cells/conditions
ev[,data.cols] <- ev[,data.cols] / t(apply(ev[,data.cols], MARGIN=1, FUN=median, na.rm=TRUE))

# take median of J cells and median of H cells
for(i in exps) {
  cat('\r', match(i, exps), '/', length(exps), '  ')
  flush.console()
  #ev[ev$exp==i]
  inds <- ev$exp==i
  if(sum(inds) <= 0) next
  # only take single cell channels
  j.inds <- desc[desc$Exp==i & desc$Type=='J' & 
                   !is.na(desc$Type) & desc$Quantity < 10,'ch'] + data.cols[1] - 1
  h.inds <- desc[desc$Exp==i & desc$Type=='H' & 
                   !is.na(desc$Type) & desc$Quantity < 10,'ch'] + data.cols[1] - 1
  ev[inds,'J'] <- apply(ev[inds,j.inds], MARGIN=1, FUN=median, na.rm=TRUE)
  ev[inds,'H'] <- apply(ev[inds,h.inds], MARGIN=1, FUN=median, na.rm=TRUE)
  ev[inds, 'JHvar'] <- apply(ev[inds,data.cols], MARGIN=1, FUN=var, na.rm=TRUE)
}
# remove rows that don't have J and H to compare against
ev <- ev[apply(is.na(ev[,c('J','H')]), MARGIN=1, FUN=sum)==0,]

# count # of peptides per protein in our entire dataset
prot.map <- as.data.frame(table(ev$Proteins))
prot.map <- prot.map[order(prot.map$Freq, decreasing=TRUE),]

# only use proteins with more than 10 peptides attached to it
prots <- as.character(prot.map$Var1[prot.map$Freq >= 20])

cors <- as.data.frame(t(sapply(prots, FUN=function(prot) {
  cat('\r', match(prot, prots), '/', length(prots), '                           ')
  flush.console()
  
  # get PSMs for this protein
  
  # ORIGINAL 
  # -- originally good IDs (PEP < 0.05), ignoring result of RTLib method
  prot.data <- ev[ev$Proteins==prot & ev$PEP < 0.05,data.cols]
  
  # NEW
  # -- originally bad IDs, upgraded to good ID (PEP > 0.05 & PEP.new < 0.05)
  # -- set of NEW and ORIGINAL should be disjoint
  prot.data.new <- ev[ev$Proteins==prot & ev$PEP > 0.05 & ev$PEP.new < 0.05,data.cols]
  
  # TOTAL
  # -- good IDs after bayesian update, regardless of if it was upgraded beyond threshold or not
  # -- **NOT** the union of ORIGINAL and NEW (but does contain all of NEW), 
  # -- as some PSMs will be downgraded below PEP threshold via. the RTLib method
  # --
  # -- reminder: PEP.updated is PEP superceded by PEP.new, so that there are no NAs
  # -- that is, a PSM that was updated by the bayesian update will have its PEP replaced by PEP.new
  prot.data.tot <- ev[ev$Proteins==prot & ev$PEP.updated < 0.05,data.cols]
  
  # init output
  prot.cor <- NULL
  prot.cor.new <- NULL
  prot.cor.tot <- NULL
  
  # lets say... must have more than 10 observations to have a meaninful correlation matrix
  # not enough data can push extremely low/high correlations, especially when so many
  # of our columns are NA anyways
  #
  # the pairwise.complete.obs method of the cor() function can exacerbate 
  # this effect if there aren't enough observations
  cor.size.thresh <- 10 
  
  if(nrow(prot.data) > cor.size.thresh) {
    # pairwise complete to account to missing/NA values
    # results in matrix that is not positive semi-definite, and rows that may contain NA
    prot.cor <- cor(t(prot.data), use='pairwise.complete.obs', method='pearson')
    # set diagonal to NA
    diag(prot.cor) <- NA
  }
  # do the same for prot.data.new, prot.data.tot
  if(nrow(prot.data.new) > cor.size.thresh) {
    prot.cor.new <- cor(t(prot.data.new), use='pairwise.complete.obs', method='pearson')
    diag(prot.cor.new) <- NA
  }
  if(nrow(prot.data.tot) > cor.size.thresh) {
    prot.cor.tot <- cor(t(prot.data.tot), use='pairwise.complete.obs', method='pearson')
    diag(prot.cor.tot) <- NA
  }
  
  # sometimes either prot.data, new, or tot will have no rows
  # in that case just output NA and have the density() function
  # or whatever is evaluating this later exclude it from calculations
  c(ifelse(nrow(prot.data) > cor.size.thresh, median(prot.cor, na.rm=TRUE), NA),
    ifelse(nrow(prot.data.new) > cor.size.thresh, median(prot.cor.new, na.rm=TRUE), NA),
    ifelse(nrow(prot.data.tot) > cor.size.thresh, median(prot.cor.tot, na.rm=TRUE), NA))
})))
# clean up the resultant data frame
cors$Protein <- rownames(cors)
rownames(cors) <- NULL
colnames(cors) <- c('Original', 'New', 'Total', 'Protein')
cors <- cors[,c(4, 1, 2, 3)]

# get null distribution for correlation
# pairwise correlation between peptides from random proteins
cors$Null <- NA
set.seed(1)
# take n at a time
# n = median # of observations for this set of proteins
n = median(prot.map$Freq[prot.map$Var1 %in% prots])
for(i in 1:nrow(cors)) {
  #prot.inds <- which(ev$Proteins %in% prots)
  #prot.data <- ev[sample(prot.inds, size=n), data.cols]
  prot.data <- ev[sample(1:nrow(prot.data), size=n), data.cols]
  prot.cor <- cor(t(prot.data), use='pairwise.complete.obs', method='pearson')
  # set diagonal to NA
  diag(prot.cor) <- NA
  cors$Null[i] <- median(prot.cor, na.rm=TRUE)
}

cors.m <- melt(cors, id.var='Protein', variable.name='type', value.name='Correlation')
cors.m %>%
  filter(type != 'Total') %>%
ggplot(aes(x=Correlation, color=type)) +
  #stat_density(geom='line') +
  #stat_density() +
  geom_line(stat='density') +
  #geom_hline(yintercept=0, color='black') +
  scale_x_continuous(limits=c(0.25,1)) +
  scale_color_manual(labels=c(
    'Original (PEP < 0.05)', 'New (PEP > 0.05 & PEP.new < 0.05)', 'Null Distribution'
  ), values=c('blue', 'green', 'red')) +
  labs(title=paste0('RTLib: Peptide-Centric Fit for 19A'),
       color='Type',
       x='Peptide Correlation w.r.t. Protein')



# get null distribution for correlation
# pairwise correlation between peptides from random proteins
null.cors <- vector(length=nrow(cors))
set.seed(1)
par(mfrow=c(2,2))
for(i in 1:4) {
  # take 50 at a time
  # randomize this quantity????
  #n = 50
  # n = median # of observations for this set of proteins
  n = median(prot.map$Freq[prot.map$Var1 %in% prots])
  
  for(i in 1:length(cors)) {
    prot.inds <- which(ev$Proteins %in% prots)
    prot.data <- ev[sample(prot.inds, size=n), data.cols]
    prot.cor <- cor(t(prot.data), use='pairwise.complete.obs', method='pearson')
    # set diagonal to NA
    diag(prot.cor) <- NA
    null.cors[i] <- median(prot.cor, na.rm=TRUE)
  }
  
  #plot(density(null.cors), col='red')
  #lines(density(cors), col='blue')
  plot(density(cors), col='blue', 
       main='Correlation=Blue\nNull=Red')
  lines(density(null.cors), col='red')
  #lines(c(0,0),c(0,1e5),type='l')
}

par(mfrow=c(1,1))

ev.a <- ev[ev$Proteins=='sp|P13639|EF2_HUMAN',]
ev.cor <- cor(t(ev.a[,data.cols]), use='pairwise.complete.obs')
#ev.cor <- cor(t(ev.a[,data.cols]), use='everything')
# remove rows/cols with all NAs
inds.remove <- which(apply(ev.cor, MARGIN=1, FUN=sum, na.rm=TRUE)==0)

# remove rows that contain NA after this (unlucky pairwise comparisons), 
# not counting the rows that were all NAs in the first place
# keep track for indexing purposes
inds.remove <- sort(c(inds.remove, 
                     which(apply(is.na(ev.cor), MARGIN=1, FUN=sum) > length(inds.remove))))

ev.cor.dist <- as.dist((1-ev.cor[-inds.remove,-inds.remove])^2)
ev.hclust <- hclust(ev.cor.dist)

image(as.matrix(ev.cor.dist)[ev.hclust$order,ev.hclust$order])
  